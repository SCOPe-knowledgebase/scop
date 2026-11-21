/*
 * Software to build and maintain SCOPe, https://scop.berkeley.edu/
 *
 * Copyright (C) 2026 The Regents of the University of California
 *
 * For feedback, mailto:scope@compbio.berkeley.edu
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * Version 2.1 of the License, or (at your option) any later version.
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 * USA
 */
package gov.lbl.scop.app;

import java.io.*;
import java.sql.*;
import gov.lbl.scop.local.LocalSQL;

import java.util.*;

import java.util.logging.FileHandler;
import java.util.logging.SimpleFormatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.nio.file.DirectoryStream;
import org.apache.commons.io.FileUtils;
import org.strbio.mol.*;
import org.strbio.mol.lib.AlignmentParameters;
import org.strbio.mol.lib.GapAffine;
import org.strbio.mol.lib.ScoreID;
import org.strbio.mol.lib.ScoreList;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class FatcatPair {

    /*
    runs one structural comparison (via FATCAT and PostFATCAT) between a specified domain and alphafold model.
    stores the results in SQL tables scop_node_fatcat_af, scop_node_fatcat_af_gap, and scop_node_fatcat_af_matrix.
    also stores copies of the resulting .aln files in /alignments/, deletes superposition files.
     */

    // query = af/pdb/uniprot = seq1 = chain 1
    // hit = scop domain = seq2 = chain 2

    final static String COMMON_DIR = "/lab/proj/astral/fatcat/2.08/";
    final static String TEMP_DIR = "/tmp/";

    final static String SCOP_DIR = "/lab/proj/astral/pdbstyle/2.08/";
    final static String AF_DIR = "/lab/db/model_structures/alphafolddb/v4/pdb/";
    
    final static String ALIGNMENTS_DIR = COMMON_DIR + "alignments/";
    final static String MATRICES_DIR = COMMON_DIR + "matrices/";
    final static String SQL_DIR = COMMON_DIR + "sql/";


    final static String FATCAT = "/lab/app/FATCAT/cur/FATCAT-dist/FATCATMain/FATCAT";
    final static String PostFATCAT = "/lab/app/FATCAT/cur/FATCAT-dist/FATCATMain/PostFATCAT";

    PreparedStatement stmt1;
    PreparedStatement stmt2;
    PreparedStatement stmt3;

    String query; // alphafold
    String hit; // domain
    int query_id;
    int node_id;

    int af_len;
    int scop_len;
    int ini_len;
    double ini_rmsd;
    int opt_equ;
    double opt_rmsd;
    double chain_rmsd;
    double score;
    int align_len;
    int gaps;
    double p_val;
    int afp_num;
    double identity;
    double similarity;
    boolean shortMatch;

    int af_start;
    int af_length;
    String scop_start;
    int scop_start_index;
    int scop_length;

    double[][] rotLst;
    double[] transLst;

    /** positions and gaps of a retained hit: model side as 0-indexed residue
     *  positions (residue number - 1), domain side as 0-indexed RAF positions */
    public ChainDomainRafPositions.Result rafPositions;

    String fatcatSeq; // sequence fatcat reads from the domain
    String rafSeq; // sequence of the domain from raf

    FatcatPairs fcs;
    String afSeq;
    String alnSeq;
    String scopSeq;

    String alignmentFile;
    String superpositionFile;
    String matrixFile;
    String scopFile; // will be different depending on the scop domain you're accessing
    String afFile;
    String subdir;

    String afHash;
    String scopHash;

    boolean comparison_completed;
    boolean matrix_parsed;
    boolean alignment_file_read;
    boolean files_archived;
    boolean temp_dir_cleared;

    String sql;
    Logger logger;
    BufferedWriter SQLWriter;

    static String align1 = "Align.*\\.pdb (\\d+) with .*\\.ent (\\d+)";
    static Pattern pAlign1 = Pattern.compile(align1); // af_len, scop_len
    static String align2 = ".* ini-len (\\d+) ini-rmsd (\\d+\\.\\d+) opt-equ (\\d+) opt-rmsd (\\d+\\.\\d+) chain-rmsd (\\d+\\.\\d+) Score (\\d+\\.\\d+) align-len (\\d+) gaps (\\d+)"; // opt-rmsd, score, align-len, gaps
    static Pattern pAlign2 = Pattern.compile(align2); // ini_len, ini_rmsd, opt_equ, opt_rmsd, chain_rmsd, score, align_len, gaps
    static String align3 = "P-value ([0-9]+.[0-9]+e[+-]*[0-9]+) Afp-num (\\d+) Identity (\\d+\\.\\d+)% Similarity (\\d+\\.\\d+)%";
    static Pattern pAlign3 = Pattern.compile(align3); // p-val, afp_num, identity, similarity
    static String alignIndex = "Chain \\d+:\\s+(-*\\d+[A-Z]*) [-\\w]+";
    static Pattern pAlignIndex = Pattern.compile(alignIndex);
    static Pattern pAlignUniprot = Pattern.compile("AF-([\\w\\d]+)-F[\\d]+");

    public FatcatPair(int query_id, PreparedStatement stmt1, PreparedStatement stmt2, PreparedStatement stmt3, FatcatPairs fcs) throws Exception {

        String af = FatcatPairs.alphaFolds.get(query_id);

        this.query = getAfFromPath(af); // pdb or uniprot or chain 1
        this.hit = fcs.domain; // scop domain or chain 2
        this.query_id = query_id;
        this.node_id = fcs.node_id;

        this.stmt1 = stmt1;
        this.stmt2 = stmt2;
        this.stmt3 = stmt3;

        // fatcat stats
        this.af_len = 0;
        this.scop_len = 0;
        this.ini_len = 0;
        this.ini_rmsd = 0;
        this.opt_equ = 0;
        this.opt_rmsd = 0;
        this.chain_rmsd = 0;
        this.score = 0;
        this.align_len = 0;
        this.gaps = 0;
        this.p_val = 0;
        this.afp_num = 0;
        this.identity = 0;
        this.similarity = 0;
        this.shortMatch = false;

        this.af_start = 0;
        this.af_length = 0;
        this.scop_start = "0";
        this.scop_start_index = -1;
        this.scop_length = 0;

        // matrix coordinates
        this.rotLst = new double[3][3];
        this.transLst = new double[3];

        // gap information for both the query pdb and hit scop domain
        this.rafPositions = null;

        // progress checks
        this.comparison_completed = false;
        this.matrix_parsed = false;
        this.alignment_file_read = false;
        this.files_archived = false;
        this.temp_dir_cleared = false;

        this.fcs = fcs;
        this.afSeq = "";
        this.alnSeq = "";
        this.scopSeq = "";

        Matcher matcher = pAlignUniprot.matcher(this.query);
        String uniprot = null;
        if (matcher.find()) {
            uniprot = matcher.group(1);
        } else if (!matcher.find()) {
            throw new Exception("no regex match uniprot, line " + 173);
        }

        this.sql = "";
        this.SQLWriter = fcs.SQLWriter;
        this.logger = fcs.logger;

        this.afHash = uniprot.substring(uniprot.length() - 3, uniprot.length() - 1);
        this.scopHash = fcs.scopHash;


        // shortcuts for the aln file name, superposed file name, and file name of the scop domain
        this.alignmentFile = this.query + "_" + this.hit + ".aln"; // in /temp/domain-0/ during, then moved to alignments_dir
        this.matrixFile = this.query + "_" + this.hit + ".matrix"; // in /temp/domain-0/ during, then moved to matrices_dir
        this.scopFile = String.valueOf(fcs.tempDomain); // copy from /pdbstyle/. deleted at the end
        this.afFile = AF_DIR + afHash + "/" + this.query + "-model_v4.pdb";
        this.subdir = fcs.subdir; // subdir between *_DIR and [af]_[scop].*

    }

    public FatcatPair(int query_id, FatcatPairs fcs) throws Exception {

        String af = FatcatPairs.alphaFolds.get(query_id);

        this.query = getAfFromPath(af); // pdb or uniprot or chain 1
        this.hit = fcs.domain; // scop domain or chain 2
        this.query_id = query_id;
        this.node_id = fcs.node_id;

        // fatcat stats
        this.af_len = 0;
        this.scop_len = 0;
        this.ini_len = 0;
        this.ini_rmsd = 0;
        this.opt_equ = 0;
        this.opt_rmsd = 0;
        this.chain_rmsd = 0;
        this.score = 0;
        this.align_len = 0;
        this.gaps = 0;
        this.p_val = 0;
        this.afp_num = 0;
        this.identity = 0;
        this.similarity = 0;
        this.shortMatch = false;

        this.af_start = 0;
        this.af_length = 0;
        this.scop_start = "0";
        this.scop_start_index = -1;
        this.scop_length = 0;

        // matrix coordinates
        this.rotLst = new double[3][3];
        this.transLst = new double[3];

        // gap information for both the query pdb and hit scop domain
        this.rafPositions = null;

        // progress checks
        this.comparison_completed = false;
        this.matrix_parsed = false;
        this.alignment_file_read = false;
        this.files_archived = false;
        this.temp_dir_cleared = false;

        this.fcs = fcs;
        this.afSeq = "";
        this.alnSeq = "";
        this.scopSeq = "";

        Matcher matcher = pAlignUniprot.matcher(this.query);
        String uniprot = null;
        if (matcher.find()) {
            uniprot = matcher.group(1);
        } else if (!matcher.find()) {
            throw new Exception("no regex match uniprot, line " + 173);
        }

        this.afHash = uniprot.substring(uniprot.length() - 3, uniprot.length() - 1);
        this.scopHash = fcs.scopHash;

        this.sql = "";
        this.SQLWriter = fcs.SQLWriter;
        this.logger = fcs.logger;

        // shortcuts for the aln file name, superposed file name, and file name of the scop domain
        this.alignmentFile = this.query + "_" + this.hit + ".aln"; // in /temp/domain-0/ during, then moved to alignments_dir
        this.matrixFile = this.query + "_" + this.hit + ".matrix"; // in /temp/domain-0/ during, then moved to matrices_dir
        this.scopFile = String.valueOf(fcs.tempDomain); // copy from /pdbstyle/. deleted at the end
        this.afFile = AF_DIR + afHash + "/" + this.query + "-model_v4.pdb";
        this.subdir = fcs.subdir; // subdir between *_DIR and [af]_[scop].*

    }

    /*
    Run FATCAT and PostFATCAT, read .aln and .matrix files.
    */
    public void compare() throws Exception {

        File alnFile = new File(fcs.TempDir + this.alignmentFile);
        File matrixFile = new File(fcs.TempDir + this.matrixFile);

        if (!new File(this.afFile).exists()) {
            throw new Exception("AF_DIR + " + this.afFile + " does not exist");
        }
        if (!new File(this.scopFile).exists()) {
            throw new Exception("SCOP_DIR + " + this.scopFile + " does not exist");
        }
        if (alnFile.exists()) {
            comparison_completed = true;
            this.readAln(alnFile);

        } else {

            String[] command = {
                    FATCAT,
                    "-p1", this.afFile,
                    "-p2", this.scopFile,
                    "-o", fcs.TempDir + this.alignmentFile.replace(".aln", ""), // (2) run FATCAT to tmp
                    "-r", "-m", "-t"
            };

            // System.out.println(String.join(" ", command));

            logger.info("running command " + String.join(" ", command));
            ProcessBuilder builder = new ProcessBuilder(command);
            builder.directory(new File("/"));
            Process p = builder.start();
            p.waitFor();
            p.destroy();

            if (!alnFile.exists()) {
                throw new Exception(alnFile + " not generated");
            }
            this.readAln(alnFile);
            this.comparison_completed = true;
            if (shortMatch) {
                return;
            }

            String[] command2 = {
                    PostFATCAT,
                    "-a", fcs.TempDir + this.alignmentFile,
                    "-p1", this.afFile,
                    "-p2", this.scopFile,
                    "-n1", this.query + "-model_v4.pdb",
                    "-n2", this.hit + ".ent",
                    "-c", fcs.TempDir + this.superpositionFile,
                    "-m", fcs.TempDir + this.matrixFile
            };

            // System.out.println(String.join(" ", command2));

            logger.info("running command " + String.join(" ", command2));
            builder = new ProcessBuilder(command2);
            builder.directory(new File("/"));
            Process p2 = builder.start();
            p2.waitFor();
            p2.destroy();

            // Read the matrix file
            if (!matrixFile.exists()) {
                throw new Exception(matrixFile + " not generated\n" + String.join(" ", command2));
            }
            this.readMatrix(matrixFile);

            int i = 0;
            while (transLst[0] == 0.0 && transLst[1] == 0.0 && transLst[2] == 0.0 && rotLst[0][0] == 0.0 && rotLst[1][1] == 0.0) { // make sure PostFATCAT runs
                if (i > 5) {
                    throw new Exception("PostFATCAT failed for " + this.query + " x " + this.scopFile);
                }
                System.out.println("repeating PostFATCAT");
                p2 = builder.start();
                p2.waitFor();
                p2.destroy();
                i++;
            }

            this.matrix_parsed = true;
        }
    }


    /*
    Read FATCAT output .aln file, store values in instance variables.
     */
    public void readAln(File result) throws Exception {

        logger.info("readAln " + result.toString());

        // get fatcatSeq
        Scanner s2 = new Scanner(result);
        while (s2.hasNextLine()) {
            String line2 = s2.nextLine();
            if (line2.startsWith("Short match")) {
                shortMatch = true;
                logger.info("readAln found short match");
            }
            else if (line2.startsWith("Chain 1")) {
                s2.nextLine();
                fatcatSeq += s2.nextLine().substring(13);
                logger.info("readAln found seq info");
            }
        }
        s2.close();

        if (!shortMatch) {
            Scanner s = new Scanner(result);
            logger.info("parsing aln");

            int i = 0;
            while (s.hasNextLine()) {
                String line = s.nextLine();

                if (i == 0) {
                    Matcher match = pAlign1.matcher(line);
                    if (match.find()) {
                        this.af_len = Integer.parseInt(match.group(1));
                        this.scop_len = Integer.parseInt(match.group(2));
                    } else if (!match.find()) {
                        throw new Exception("no regex match lens, line " + 300);
                    }
                } else if (i == 1) {
                    Matcher match = pAlign2.matcher(line);
                    if (match.find()) {
                        this.ini_len = Integer.parseInt(match.group(1));
                        this.ini_rmsd = Double.parseDouble(match.group(2));
                        this.opt_equ = Integer.parseInt(match.group(3));
                        this.opt_rmsd = Double.parseDouble(match.group(4));
                        this.chain_rmsd = Double.parseDouble(match.group(5));
                        this.score = Double.parseDouble(match.group(6));
                        this.align_len = Integer.parseInt(match.group(7));
                        this.gaps = Integer.parseInt(match.group(8));
                    } else if (!match.find()) {
                        throw new Exception("no regex match stats, line " + 319);
                    }
                } else if (i == 2) {
                    Matcher match = pAlign3.matcher(line);
                    if (match.find()) {
                        this.p_val = Double.parseDouble(match.group(1));
                        this.afp_num = Integer.parseInt(match.group(2));
                        this.identity = Double.parseDouble(match.group(3));
                        this.similarity = Double.parseDouble(match.group(4));
                    } else if (!match.find()) {
                        throw new Exception("no regex match stats2, line " + 329);
                    }
                }
                if (line.startsWith("Chain 1")) {
                    String line_gap = s.nextLine();
                    String line_scop = s.nextLine();
                    s.nextLine();
                    parseLines(line, line_gap, line_scop);
                }
                i++;
            }
            s.close();
        }
        s2.close();
    }

    /*
    Read PostFATCAT output [AF]_[domain].matrix file, storing the transformation matrix coordinates into rotLst and transLst.
     */
    public void readMatrix(File result) throws Exception {

        logger.info("readMatrix " + result.toString());

        Scanner s = new Scanner(result);
        double[] Tfatcat = new double[3];
        String[] lst;

        while (s.hasNextLine()) {
            String line = s.nextLine();
            if (line.startsWith("rotate") || line.startsWith("translate")) {

                lst = line.split("\\s");
                int l = 1;

                if (line.startsWith("rotate")) {
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            rotLst[i][j] = Double.parseDouble(lst[l]);
                            l++;
                        }
                    }
                } else if (line.startsWith("translate")) { // + extra transformations to get the correct translation vector
                    for (int i = 0; i < 3; i++) {
                        Tfatcat[i] = Double.parseDouble(lst[l]);
                        l++;
                    }
                    for (int i = 0; i < 3; i++) { // Tcorrect = R x -Tfatcat
                        transLst[i] = 0;
                        for (int j = 0; j < 3; j++) {
                            transLst[i] += rotLst[i][j] * -Tfatcat[j];
                        }
                    }
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            transLst[i] = Math.round(transLst[i] * 1000.0) / 1000.0;
        }
        s.close();
    }

    /*
    Write insert statements into a .sql file (one per job)
     */
    private boolean shouldStoreDetail() {
        // These are the shared public viewer/API thresholds.  A FATCAT
        // P-value reported as zero is an underflow/serialization case and is
        // therefore retained as the most significant value, not discarded.
        return p_val >= 0.0 && p_val <= 0.05 && ini_rmsd < 4.0 && score > 0.0;
    }

    public void storeSQL() throws Exception {

        boolean detail = shouldStoreDetail() && rafPositions != null && rafPositions.hasAlignment();

        // stmt1 (scop_node_fatcat_af); the domain positions are stored for retained hits only
        String stmt1 = String.format(
                "INSERT INTO scop_node_fatcat_af (af_id, node_id, af_len, scop_len, ini_len, ini_rmsd, opt_equ, opt_rmsd, " +
                "chain_rmsd, score, align_len, gaps, p_val, afp_num, identity, similarity, af_start, af_length, " +
                "domain_start, domain_end) VALUES " +
                "(%d, %d, %d, %d, %d, %.6f, %d, %.6f, %.6f, %.6f, %d, %d, %s, %d, %.6f, %.6f, %d, %d, %s, %s);",
                query_id,
                node_id,
                af_len,
                scop_len,
                ini_len,
                ini_rmsd,
                opt_equ,
                opt_rmsd,
                chain_rmsd,
                score,
                align_len,
                gaps,
                Double.toString(p_val),
                afp_num,
                identity,
                similarity,
                af_start,
                af_length,
                detail ? String.valueOf(rafPositions.domainStart) : "NULL",
                detail ? String.valueOf(rafPositions.domainEnd) : "NULL"
        );

        String detailSQL = "";
        if (detail) {
            // store auto-incremented id from scop_node_fatcat_af insertion
            String auto_id = "SET @hit_id = LAST_INSERT_ID();";

            // stmt2 (scop_node_fatcat_af_gap): model residue numbers, domain RAF positions
            StringBuilder sb = new StringBuilder();
            if (!rafPositions.gaps.isEmpty()) {
                sb.append("INSERT INTO scop_node_fatcat_af_gap ");
                sb.append("(hit_id, gap_id, gap_type, gap_start_af, gap_end_af, gap_start_domain, gap_end_domain) VALUES\n");
                for (int n = 0; n < rafPositions.gaps.size(); n++) {
                    ChainDomainRafPositions.Gap g = rafPositions.gaps.get(n);
                    sb.append(String.format(
                            "(@hit_id, %d, %d, %s, %s, %s, %s)",
                            g.gapId,
                            afGapType(g.gapType),
                            g.chainStart == null ? "NULL" : String.valueOf(g.chainStart + 1),
                            g.chainEnd == null ? "NULL" : String.valueOf(g.chainEnd + 1),
                            g.domainStart == null ? "NULL" : String.valueOf(g.domainStart),
                            g.domainEnd == null ? "NULL" : String.valueOf(g.domainEnd)
                    ));
                    sb.append(n < rafPositions.gaps.size() - 1 ? ",\n" : ";");
                }
            }
            String stmt2 = sb.toString();

            // stmt3 (scop_node_fatcat_af_matrix)
            sb = new StringBuilder();
            sb.append("INSERT INTO scop_node_fatcat_af_matrix ");
            sb.append("(hit_id, r11, r12, r13, r21, r22, r23, r31, r32, r33, t1, t2, t3) VALUES\n");
            sb.append("(@hit_id");
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sb.append(String.format(", %.10f", rotLst[i][j]));
                }
            }
            for (int m = 0; m < 3; m++) {
                sb.append(String.format(", %.10f", transLst[m]));
            }
            sb.append(");");
            String stmt3 = sb.toString();
            detailSQL = auto_id + "\n";
            if (!stmt2.isEmpty()) {
                detailSQL += stmt2 + "\n";
            }
            detailSQL += stmt3 + "\n";
        }

        // formatting
        this.sql = stmt1 + "\n" + detailSQL;
        SQLWriter.write(sql);
        SQLWriter.flush();
    }

    /*
    Store the .aln file and superposed pdb file into their archive directories, and clear the temp directory
     */
    public void archiveAndClear() throws Exception {

        Path oldAlnPath = Paths.get(fcs.TempDir + this.alignmentFile);
        Path newAlnPath = Paths.get(fcs.AlnDir + this.alignmentFile);

        logger.info("archiveAndClear--moving aln " + oldAlnPath + " to " + newAlnPath);
        Files.move(oldAlnPath, newAlnPath, StandardCopyOption.REPLACE_EXISTING);

        if (Files.exists(newAlnPath)) {
            if (Files.exists(oldAlnPath)) {
                Files.delete(oldAlnPath);
                if (!Files.exists(oldAlnPath)) {
                    logger.info("File was successfully moved to " + newAlnPath);
                }
                else {
                    logger.severe("failed to delete " + oldAlnPath);
                }
            }
        } else {
            logger.severe("File move failed: " + newAlnPath + " does not exist.");
        }

        if (!shortMatch) {

            Path oldMatrixPath = Paths.get(fcs.TempDir + this.matrixFile);
            Path newMatrixPath = Paths.get(fcs.MatrixDir + this.matrixFile);
            logger.info("archiveAndClear--moving matrix " + oldMatrixPath + " to " + newMatrixPath);
            Files.move(oldMatrixPath, newMatrixPath, StandardCopyOption.REPLACE_EXISTING);

            if (Files.exists(newMatrixPath)) {
                if (Files.exists(oldMatrixPath)) {
                    Files.delete(oldMatrixPath);
                    if (!Files.exists(oldMatrixPath)) {
                        logger.info("File was successfully moved to " + newMatrixPath);
                    } else {
                        logger.severe("failed to delete " + oldMatrixPath);
                    }
                }
            } else {
                logger.severe("File move failed: " + newMatrixPath + " does not exist.");
            }

            // clear the temp dir
            logger.info("clearing temp dir " + fcs.TempDir);

            String logFileName = hit + "-" + fcs.segment_id + ".log";

            try (DirectoryStream<Path> stream = Files.newDirectoryStream(Paths.get(fcs.TempDir))) {
                for (Path path : stream) {
                    if (Files.isRegularFile(path)) {
                        String fileName = path.getFileName().toString();

                        // Skip log file, .ent file, and NFS ghost files
                        if (fileName.equals(logFileName) || fileName.endsWith(".ent") || fileName.startsWith(".nfs")) {
                            continue;
                        }

                        Files.delete(path);
                    }
                }
            }
        }
    }

    /*
    collect scopSeq, alnSeq, afSeq, prepare for gap processing.
     */
    public void parseLines(String af_line, String line2, String scop_line) throws Exception {

        logger.info("parseLines\n" + af_line + "\n" + line2 + "\n" + scop_line);

        // get starts of alignments in terms of AF and SCOP residue indices
        Matcher matchIndex_af = pAlignIndex.matcher(af_line);
        Matcher matchIndex_scop = pAlignIndex.matcher(scop_line);

        int af_index;
        int scop_index;
        String scop_start_str;

        if (matchIndex_af.find()) {
            af_index = Integer.parseInt(matchIndex_af.group(1));
        } else {
            throw new Exception("no regex match af_index, line " + 502);
        }
        if (matchIndex_scop.find()) {
            scop_start_str = matchIndex_scop.group(1);
            scop_index = new ArrayList<>(fcs.rafMap.keySet()).indexOf(scop_start_str);
            if (scop_start_index < 0) {
                scop_start_index = scop_index;
            }
        } else {
            throw new Exception("no regex match scop_index, line " + 507);
        }

        // initialize indices
        if (this.af_start == 0 || this.scop_start == "0") {
            this.af_start = af_index;
            this.scop_start = scop_start_str;
        }

        char[] af_line_arr = af_line.toCharArray();
        char[] scop_line_arr = scop_line.toCharArray();

        for (int i = 14; i < af_line.length(); i++) {
            char c_af = af_line_arr[i];
            char c_scop = scop_line_arr[i];
            if (c_af != '-') {
                af_index++;
            }
            if (c_scop != '-') {
                scop_index++;
            }
        }
        this.afSeq += new String(af_line_arr).substring(14).trim();
        this.alnSeq += line2.substring(14);
        this.scopSeq += new String(scop_line_arr).substring(14).trim();

        this.af_length = af_index - af_start;
        this.scop_length = scop_index - scop_start_index;
    }


    public void clearTemp() throws Exception {
        FileUtils.cleanDirectory(new File(TEMP_DIR));
    }


    public static void connectDB() throws Exception {
        LocalSQL.connectRW();
    }

    /*
    returns the alphafold identifier given the path
     */
    public String getAfFromPath(String path) throws Exception {

        Pattern pattern = Pattern.compile("AF-[A-Za-z0-9]+-F\\d+");
        Matcher matcher = pattern.matcher(path);
        if (matcher.find()) {
            return matcher.group(0);
        }
        return null;
    }

    /** scop_node_fatcat_af_gap_type id for a gap type of
     *  ChainDomainRafPositions (model in the place of the chain).  A model
     *  has no unobserved residues or residues of unknown type. */
    static int afGapType(int type) throws Exception {
        switch (type) {
            case 1: return 1;    // AlphaFold only
            case 2: return 2;    // SCOP domain only
            case 3: return 3;    // both
            case 6: return 4;    // domain residue of unknown type, domain only
            case 7: return 5;    // domain residue of unknown type, both
            case 10: return 6;   // unobserved domain residue, domain only
            case 11: return 7;   // unobserved domain residue, both
        }
        throw new Exception("gap type " + type + " cannot occur for a model");
    }

    /** RAF-style body for a model: one record per residue of the file, all
     *  observed, so that position i is residue number i + 1. */
    static String modelRafBody(ArrayList<String[]> residues) throws Exception {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < residues.size(); i++) {
            String[] r = residues.get(i);
            if (!r[0].equals(String.valueOf(i + 1)))
                throw new Exception("model residues are not numbered 1..N (residue " + r[0] + " at position " + (i + 1) + ")");
            char c = Character.toLowerCase(r[1].charAt(0));
            sb.append(String.format("%4s %c%c", r[0], c, c));
        }
        return sb.toString();
    }

    /** Derive the positions and gaps of a retained hit from the .aln rows
     *  and the RAF of the domain, as ChainDomainFatcatPair does for chains. */
    public void deriveRafPositions() throws Exception {
        rafPositions = null;
        if (shortMatch || !shouldStoreDetail()) return;
        java.util.List<String> lines = Files.readAllLines(Paths.get(fcs.TempDir + this.alignmentFile),
                                                          java.nio.charset.StandardCharsets.ISO_8859_1);
        ArrayList<ChainDomainRafPositions.Row> rows = ChainDomainRafPositions.parseAln(lines, new boolean[1]);
        ArrayList<String[]> modelRes = ChainDomainRafPositions.fileResidues(Paths.get(this.afFile));
        if (fcs.domainResidues == null)
            fcs.domainResidues = ChainDomainRafPositions.fileResidues(Paths.get(this.scopFile));
        int[] range = ChainDomainRafPositions.domainRange(fcs.rafBody, fcs.domainDescription);
        rafPositions = ChainDomainRafPositions.derive(modelRafBody(modelRes), modelRes, fcs.rafBody,
                                                      fcs.domainResidues, range[0], range[1], rows);
        int mismatches = rafPositions.chainMismatches + rafPositions.domainMismatches;
        if (mismatches > 0)
            throw new Exception(mismatches + " residue letters of the FATCAT alignment of " +
                                this.query + " x " + this.hit + " disagree with the structure files");
    }

    public void runComparison() throws Exception {
        
        this.compare();
        this.deriveRafPositions();
        this.archiveAndClear();
        this.storeSQL();
        
    }

    public static void main(String[] args) throws Exception {

//        List<Integer> nodeIds = Arrays.asList(
//                2711103, 2711113, 2711118, 2711128, 2711136,
//                2711140, 2711145, 2711147, 2711152, 2711162,
//
//                2711176, 2711179, 2711183, 2711189, 2711196,
//                2711207, 2711216, 2711255, 2711332, 2711352
//        );
//
////        int idx = Integer.parseInt(args[0]);
//        int idx = 0;
//
//        LocalSQL.connectRW();
//
//        ResultSet rs = LocalSQL.createStatement().executeQuery("SHOW GRANTS");
//        while (rs.next()) {
//            System.out.println(rs.getString(1));
//        }
//
//        for (int i = 0; i < 10; i++) {
//            FatcatPairs.runComparisons(nodeIds.get(idx), i);
//        }

        
        LocalSQL.connectRW();
        int nodeID = LocalSQL.lookupNodeBySid("d1ryha_", 19);
        FatcatPairs.runTestComparisons(nodeID, 2, 1);
    }
}
