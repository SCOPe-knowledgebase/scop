/*
 * Software to build and maintain SCOPe, https://scop.berkeley.edu/
 *
 * Copyright (C) 2026 The Regents of the University of California
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * Version 2.1 of the License, or (at your option) any later version.
 */
package gov.lbl.scop.app;

import gov.lbl.scop.util.RAF;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Derive the RAF positions of a chain x domain FATCAT alignment directly
 * from the FATCAT .aln rows and the two RAF records.
 *
 * Every position is a 0-indexed index into the RAF bound to the chain
 * (astral_chain.raf_id) or to the chain of the domain.  The structural
 * alignment runs from the first to the last column FATCAT marks as
 * equivalent ("1" in the line between the two sequences).  FATCAT reads the
 * ATOM records of the first model and numbers the residues that have a CA
 * atom in file order, printing the author id of the first residue of each
 * 70-column row; so the k-th residue after that anchor is the k-th
 * following residue of the file, whose author id is known exactly and
 * looked up in the RAF.  No letters are compared to place residues.  Unaligned runs between the two ends, the residues the
 * RAF lists without coordinates inside them, and the flanks outside the
 * alignment are reported as typed gaps (the gap types of
 * scop_node_fatcat_chain_gap: 1 chain only, 2 domain only, 3 both; +4 for
 * a residue of unknown type, +8 for a residue without coordinates; chain
 * before domain in the odd/even pairs 4/6, 5/7, 8/10, 9/11).
 *
 * This replaces the RAF-versus-FATCAT sequence alignment that
 * ChainDomainFatcatPair.constructNewSeqs used until November 2026, which
 * rendered every insertion-coded residue as a slash and therefore could
 * not place chains such as 2qrtA and 3nwmA whose residues all carry one.
 */
public class ChainDomainRafPositions {

    public static class Gap {
        public int gapId;
        public int gapType;
        public Integer chainStart, chainEnd, domainStart, domainEnd;
        public int chainLength() { return chainStart==null ? 0 : chainEnd-chainStart+1; }
        public int domainLength() { return domainStart==null ? 0 : domainEnd-domainStart+1; }
    }

    public static class Result {
        public Integer chainStart, chainEnd, domainStart, domainEnd;
        public ArrayList<Gap> gaps = new ArrayList<Gap>();
        public int rows, alignedColumns, chainMismatches, domainMismatches;
        /** FATCAT printed an equivalence mark in a column where one side has
         *  no residue (its row was printed with spurious leading gaps); the
         *  extent then comes from the printed residues, as before 2026. */
        public boolean inconsistentMarks;
        /** Residues the file has but the RAF lists without coordinates. */
        public int chainSkipped, domainSkipped;
        public boolean shortMatch;
        public boolean hasAlignment() { return chainStart != null; }
    }

    /** One 70-column row of the .aln: residue id and text of the chain
     *  line, the equivalence marks, residue id and text of the domain line. */
    public static class Row {
        public String chainId, chainText, marks, domainId, domainText;
    }

    static final Pattern ROW_ID = Pattern.compile("^Chain \\d+:\\s*(-*\\d+[A-Za-z]*) ");
    static final Pattern RANGE = Pattern.compile("^(-?\\d+[A-Za-z]?)-(-?\\d+[A-Za-z]?)$");

    /** Collect the alignment rows of a FATCAT .aln file. */
    public static ArrayList<Row> parseAln(List<String> lines, boolean[] shortMatch) throws Exception {
        ArrayList<Row> rows = new ArrayList<Row>();
        for (int i = 0; i < lines.size(); i++) {
            String line = lines.get(i);
            if (line.startsWith("Short match")) { shortMatch[0] = true; continue; }
            if (!line.startsWith("Chain 1:")) continue;
            if (i + 2 >= lines.size()) throw new Exception("truncated alignment row at line " + (i + 1));
            String marks = lines.get(i + 1);
            String domainLine = lines.get(i + 2);
            Row r = new Row();
            Matcher m = ROW_ID.matcher(line);
            if (!m.find()) throw new Exception("no residue id in chain row: " + line);
            r.chainId = m.group(1);
            m = ROW_ID.matcher(domainLine);
            if (!m.find()) throw new Exception("no residue id in domain row: " + domainLine);
            r.domainId = m.group(1);
            r.chainText = line.length() > 14 ? line.substring(14) : "";
            r.domainText = domainLine.length() > 14 ? domainLine.substring(14) : "";
            r.marks = marks.length() > 14 ? marks.substring(14) : "";
            rows.add(r);
            i += 2;
        }
        return rows;
    }

    /** RAF index range covered by a domain description such as
     *  "1hqc A:243-318", "1abc A:" or "2xyz A:1-100,A:150-200"; the whole
     *  chain when the description has no usable range. */
    public static int[] domainRange(String domainRaf, String description) {
        int last = domainRaf.length() / 7 - 1;
        int lo = -1, hi = -1;
        boolean usable = description != null;
        if (usable) {
            String d = description.trim();
            int sp = d.indexOf(' ');
            if (sp >= 0) d = d.substring(sp + 1);
            for (String seg : d.split(",")) {
                int colon = seg.indexOf(':');
                String range = colon < 0 ? "" : seg.substring(colon + 1).trim();
                if (range.length() == 0) { usable = false; break; }
                Matcher m = RANGE.matcher(range);
                if (!m.matches()) { usable = false; break; }
                int a = RAF.indexOf(domainRaf, m.group(1), true);
                int b = RAF.indexOf(domainRaf, m.group(2), false);
                if (a < 0 || b < 0 || b < a) { usable = false; break; }
                if (lo < 0 || a < lo) lo = a;
                if (b > hi) hi = b;
            }
        }
        if (!usable || lo < 0) return new int[]{0, last};
        return new int[]{lo, hi};
    }

    private static boolean observed(String raf, int idx) { return raf.charAt(idx * 7 + 5) != '.'; }

    private static boolean unknownType(String raf, int idx) {
        char c = Character.toLowerCase(raf.charAt(idx * 7 + 5));
        return c == 'x' || c == 'b' || c == 'z';
    }

    static final java.util.Map<String, Character> THREE_TO_ONE = new java.util.HashMap<String, Character>();
    static {
        String[] t = {"ALA","A","ARG","R","ASN","N","ASP","D","CYS","C","GLN","Q","GLU","E","GLY","G","HIS","H",
                      "ILE","I","LEU","L","LYS","K","MET","M","PHE","F","PRO","P","SER","S","THR","T","TRP","W",
                      "TYR","Y","VAL","V"};
        for (int i = 0; i < t.length; i += 2) THREE_TO_ONE.put(t[i], t[i + 1].charAt(0));
    }

    /** The residues FATCAT sees in a PDB-format file, in file order: ATOM
     *  records of the first model, one entry per residue id (columns
     *  23-27), only residues that have a CA atom.  Each entry is
     *  {author residue id, one-letter code or X}. */
    public static ArrayList<String[]> fileResidues(java.nio.file.Path file) throws java.io.IOException {
        ArrayList<String[]> res = new ArrayList<String[]>();
        String lastKey = null; boolean lastHasCA = false; String lastName = null;
        try (java.io.BufferedReader br = java.nio.file.Files.newBufferedReader(file, java.nio.charset.StandardCharsets.ISO_8859_1)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("ENDMDL")) break;
                if (!line.startsWith("ATOM")) continue;
                if (line.length() < 27) continue;
                String key = line.substring(22, 27);
                if (!key.equals(lastKey)) {
                    if (lastKey != null && lastHasCA) res.add(new String[]{lastKey.trim(), lastName});
                    lastKey = key; lastHasCA = false;
                    Character one = THREE_TO_ONE.get(line.substring(17, 20).trim());
                    lastName = one == null ? "X" : String.valueOf(one);
                }
                if (line.length() >= 16 && line.substring(12, 16).trim().equals("CA")) lastHasCA = true;
            }
        }
        if (lastKey != null && lastHasCA) res.add(new String[]{lastKey.trim(), lastName});
        return res;
    }

    /** Map the residues of one row to RAF indices (-1 for a gap column).
     *  FATCAT drops residues whose name it does not know (UNK and the like)
     *  from the sequence but not from its residue-id array, so the ids it
     *  prints at the start of the second and later rows are stale after
     *  such a block (observed in 5oqj4, residues 89-114 UNK); only the first
     *  row's id is exact.  The residues of the whole alignment are therefore
     *  consecutive entries of the filtered residue list starting at the
     *  first row's residue; k is the position of the row's first residue in
     *  that list and the position after the row is returned in next[0].
     *  stats[0] counts residues whose letter in the file differs from the
     *  letter FATCAT printed (an alignment computed from a different file);
     *  stats[1] counts residues the RAF lists without coordinates although
     *  the file has them. */
    private static int[] mapRow(String raf, ArrayList<String[]> known, int k, String text,
                                int[] stats, int[] next) throws Exception {
        int[] idx = new int[text.length()];
        for (int c = 0; c < text.length(); c++) {
            char ch = text.charAt(c);
            if (ch == '-' || ch == ' ') { idx[c] = -1; continue; }
            if (k >= known.size()) throw new Exception("alignment runs past the end of the structure file");
            String[] r = known.get(k);
            int p = RAF.indexOf(raf, r[0], true);
            if (p < 0) throw new Exception("residue " + r[0] + " of the structure file is not in the RAF");
            if (r[1].charAt(0) != ch && r[1].charAt(0) != 'X') stats[0]++;
            if (!observed(raf, p)) stats[1]++;
            idx[c] = p;
            k++;
        }
        next[0] = k;
        return idx;
    }

    /** Position in the filtered list of the residue FATCAT printed on the
     *  first row; a residue of unknown type there maps to the next kept one. */
    private static int firstPosition(ArrayList<String[]> fileRes, ArrayList<String[]> known, String resId)
        throws Exception {
        int u = -1;
        for (int i = 0; i < fileRes.size(); i++) if (fileRes.get(i)[0].equals(resId)) { u = i; break; }
        if (u < 0) throw new Exception("residue " + resId + " is not in the structure file");
        int k = 0;
        for (int i = 0; i < u; i++) if (!fileRes.get(i)[1].equals("X")) k++;
        return k;
    }

    /** The residues FATCAT keeps: those with a standard residue name. */
    private static ArrayList<String[]> knownResidues(ArrayList<String[]> fileRes) {
        ArrayList<String[]> k = new ArrayList<String[]>();
        for (String[] r : fileRes) if (!r[1].equals("X")) k.add(r);
        return k;
    }

    private static int gapType(String chainRaf, int c, String domainRaf, int d) {
        boolean hasC = c >= 0, hasD = d >= 0;
        boolean cMissing = hasC && !observed(chainRaf, c), dMissing = hasD && !observed(domainRaf, d);
        boolean cUnknown = hasC && !cMissing && unknownType(chainRaf, c);
        boolean dUnknown = hasD && !dMissing && unknownType(domainRaf, d);
        if (hasC && !hasD) return cUnknown ? 4 : cMissing ? 8 : 1;
        if (!hasC && hasD) return dUnknown ? 6 : dMissing ? 10 : 2;
        if (cUnknown) return 5;
        if (dUnknown) return 7;
        if (cMissing) return 9;
        if (dMissing) return 11;
        return 3;
    }

    /** Flank columns: the residues of both sides outside the alignment,
     *  paired column by column against the alignment (right-aligned before
     *  it, left-aligned after it), the excess of the longer side one-sided. */
    private static void addFlank(ArrayList<int[]> stream, String chainRaf, int c0, int c1,
                                 String domainRaf, int d0, int d1, boolean before) {
        int nc = Math.max(0, c1 - c0 + 1), nd = Math.max(0, d1 - d0 + 1);
        int paired = Math.min(nc, nd);
        if (before) {
            for (int i = 0; i < nc - paired; i++) stream.add(new int[]{c0 + i, -1});
            for (int i = 0; i < nd - paired; i++) stream.add(new int[]{-1, d0 + i});
            for (int i = 0; i < paired; i++) stream.add(new int[]{c0 + (nc - paired) + i, d0 + (nd - paired) + i});
        } else {
            for (int i = 0; i < paired; i++) stream.add(new int[]{c0 + i, d0 + i});
            for (int i = paired; i < nc; i++) stream.add(new int[]{c0 + i, -1});
            for (int i = paired; i < nd; i++) stream.add(new int[]{-1, d0 + i});
        }
    }

    public static Result derive(String chainRaf, ArrayList<String[]> chainRes,
                                String domainRaf, ArrayList<String[]> domainRes,
                                int domainFirst, int domainLast, List<Row> rows) throws Exception {
        ArrayList<String[]> chainKnown = knownResidues(chainRes), domainKnown = knownResidues(domainRes);
        int[] kc = new int[1], kd = new int[1];
        if (rows.size() > 0) {
            kc[0] = firstPosition(chainRes, chainKnown, rows.get(0).chainId);
            kd[0] = firstPosition(domainRes, domainKnown, rows.get(0).domainId);
        }
        Result r = new Result();
        r.rows = rows.size();
        int[] cm = new int[2], dm = new int[2];
        ArrayList<int[]> cols = new ArrayList<int[]>();   // {chainIdx, domainIdx, aligned}
        for (int ri = 0; ri < rows.size(); ri++) {
            Row row = rows.get(ri);
            int[] ci = mapRow(chainRaf, chainKnown, kc[0], row.chainText, cm, kc);
            int[] di = mapRow(domainRaf, domainKnown, kd[0], row.domainText, dm, kd);
            int w = Math.max(ci.length, di.length);
            for (int c = 0; c < w; c++) {
                int a = c < ci.length ? ci[c] : -1;
                int b = c < di.length ? di[c] : -1;
                boolean marked = c < row.marks.length() && row.marks.charAt(c) == '1';
                if (marked && (a < 0 || b < 0)) r.inconsistentMarks = true;
                int al = (marked && a >= 0 && b >= 0) ? 1 : 0;
                cols.add(new int[]{a, b, al});
            }
        }
        r.chainMismatches = cm[0];
        r.domainMismatches = dm[0];
        r.chainSkipped = cm[1];
        r.domainSkipped = dm[1];
        int first = -1, last = -1;
        for (int i = 0; i < cols.size(); i++)
            if (cols.get(i)[2] == 1) { if (first < 0) first = i; last = i; r.alignedColumns++; }
        if (first < 0) return r;
        r.chainStart = cols.get(first)[0];
        r.chainEnd = cols.get(last)[0];
        r.domainStart = cols.get(first)[1];
        r.domainEnd = cols.get(last)[1];
        if (r.inconsistentMarks) {
            // the marks cannot be trusted column by column: use the printed
            // extent of each side (first and last residue FATCAT printed)
            int fc = -1, lc = -1, fd = -1, ld = -1;
            for (int i = 0; i < cols.size(); i++) {
                if (cols.get(i)[0] >= 0) { if (fc < 0) fc = i; lc = i; }
                if (cols.get(i)[1] >= 0) { if (fd < 0) fd = i; ld = i; }
            }
            r.chainStart = cols.get(fc)[0]; r.chainEnd = cols.get(lc)[0];
            r.domainStart = cols.get(fd)[1]; r.domainEnd = cols.get(ld)[1];
            first = Math.min(fc, fd); last = Math.max(lc, ld);
        }

        // column stream over the whole extent: {chainIdx, domainIdx}, aligned columns as {c, d, 1}
        ArrayList<int[]> stream = new ArrayList<int[]>();
        int chainLen = chainRaf.length() / 7;
        addFlank(stream, chainRaf, 0, r.chainStart - 1, domainRaf, domainFirst, r.domainStart - 1, true);
        int prevC = r.chainStart, prevD = r.domainStart;
        for (int i = first; i <= last; i++) {
            int c = cols.get(i)[0], d = cols.get(i)[1];
            if (c >= 0) { for (int p = prevC + 1; p < c; p++) stream.add(new int[]{p, -1}); prevC = c; }
            if (d >= 0) { for (int p = prevD + 1; p < d; p++) stream.add(new int[]{-1, p}); prevD = d; }
            if (cols.get(i)[2] == 1) stream.add(new int[]{c, d, 1});
            else stream.add(new int[]{c, d});
        }
        addFlank(stream, chainRaf, r.chainEnd + 1, chainLen - 1, domainRaf, r.domainEnd + 1, domainLast, false);

        // compress runs of equal type into gaps; a gap id spans the runs
        // between two aligned columns
        Gap cur = null;
        int gapId = 1;
        boolean sinceAligned = false;
        for (int[] col : stream) {
            if (col.length == 3) {           // aligned column
                if (cur != null) { r.gaps.add(cur); cur = null; }
                if (sinceAligned) gapId++;
                sinceAligned = false;
                continue;
            }
            int type = gapType(chainRaf, col[0], domainRaf, col[1]);
            if (cur != null && cur.gapType != type) { r.gaps.add(cur); cur = null; }
            if (cur == null) { cur = new Gap(); cur.gapType = type; cur.gapId = gapId; }
            if (col[0] >= 0) { if (cur.chainStart == null) cur.chainStart = col[0]; cur.chainEnd = col[0]; }
            if (col[1] >= 0) { if (cur.domainStart == null) cur.domainStart = col[1]; cur.domainEnd = col[1]; }
            sinceAligned = true;
        }
        if (cur != null) r.gaps.add(cur);
        return r;
    }
}
