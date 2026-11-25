#!/usr/bin/env python3
#
#   Copyright (C) 2026 John-Marc Chandonia
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2, or (at your option) any
#   later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301, USA.
#
#   This file is part of the PAST (PDB Archival Snapshot Toolkit) package.
#
"""Mock wwPDB file server for the PAST test harness (Python stdlib only).

usage: mockserver.py DOCROOT CONTROLDIR

Serves files below DOCROOT.  A directory is served from its index.html.
ETag is the MD5 of the file, Last-Modified its modification time.
CONTROLDIR/faults.json is read on every request and maps a URL path (or
"*") to faults:
    status            answer with this HTTP status
    truncate          send only this many bytes of the body (GET)
    etag              send this ETag instead, quoted (null: none)
    etag_raw          send this ETag header value exactly as given
    no_length         send no Content-Length (the connection is closed instead)
    no_last_modified  send no Last-Modified
    delay             sleep this many seconds first
The chosen port is written to CONTROLDIR/port; every request is appended
to CONTROLDIR/access.log as "METHOD path".
"""
import email.utils
import hashlib
import json
import os
import sys
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

DOCROOT = os.path.abspath(sys.argv[1])
CONTROL = os.path.abspath(sys.argv[2])


def faults_for(path):
    try:
        with open(os.path.join(CONTROL, "faults.json")) as fh:
            faults = json.load(fh)
    except (OSError, ValueError):
        faults = {}
    merged = dict(faults.get("*", {}))
    merged.update(faults.get(path, {}))
    return merged


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, *args):
        pass

    def answer(self, send_body):
        path = self.path.split("?")[0]
        with open(os.path.join(CONTROL, "access.log"), "a") as fh:
            fh.write("%s %s\n" % (self.command, path))
        fault = faults_for(path)
        if fault.get("delay"):
            time.sleep(fault["delay"])
        local = os.path.normpath(os.path.join(DOCROOT, path.lstrip("/")))
        if not (local == DOCROOT or local.startswith(DOCROOT + os.sep)):
            return self.error(403)
        if os.path.isdir(local):
            local = os.path.join(local, "index.html")
        if fault.get("status"):
            return self.error(fault["status"])
        if not os.path.isfile(local):
            return self.error(404)
        with open(local, "rb") as fh:
            data = fh.read()
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        no_length = bool(fault.get("no_length"))
        if no_length:
            self.send_header("Connection", "close")
        else:
            self.send_header("Content-Length", str(len(data)))
        if "etag_raw" in fault:
            self.send_header("ETag", fault["etag_raw"])
        elif "etag" in fault:
            if fault["etag"] is not None:
                self.send_header("ETag", '"%s"' % fault["etag"])
        else:
            self.send_header("ETag", '"%s"' % hashlib.md5(data).hexdigest())
        if not fault.get("no_last_modified"):
            self.send_header("Last-Modified", email.utils.formatdate(
                os.path.getmtime(local), usegmt=True))
        truncated = send_body and "truncate" in fault
        if truncated:
            self.send_header("Connection", "close")
        self.end_headers()
        if send_body:
            try:
                self.wfile.write(data[:fault["truncate"]] if truncated else data)
            except (BrokenPipeError, ConnectionResetError):
                self.close_connection = True
        if truncated or no_length:
            self.close_connection = True

    def error(self, status):
        body = b"error\n"
        self.send_response(status)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def do_GET(self):
        self.answer(True)

    def do_HEAD(self):
        self.answer(False)


def main():
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    server.daemon_threads = True
    with open(os.path.join(CONTROL, "port.tmp"), "w") as fh:
        fh.write(str(server.server_address[1]))
    os.rename(os.path.join(CONTROL, "port.tmp"), os.path.join(CONTROL, "port"))
    server.serve_forever()


if __name__ == "__main__":
    main()
