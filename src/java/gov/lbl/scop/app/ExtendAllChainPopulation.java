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

import gov.lbl.scop.local.LocalSQL;
import java.sql.*;
import java.util.*;

/**
   Adds chains imported after an all-chain population was built to
   that population and to its 90% identity set, without changing
   anything already there:

   - a new chain whose sequence is already in the population joins
     that sequence's cluster;
   - otherwise it joins the cluster of the representative it is most
     similar to (at least 90% identical), if there is one;
   - otherwise it becomes a new representative, placed after all
     existing ones, and its FATCAT comparisons are queued.

   No candidate, representative or cluster assignment is ever removed
   or reassigned.  Run "prepare", then the BLAST workers
   (AllChainBlast worker), then "assign".
*/
public class ExtendAllChainPopulation {
    final public static int PCT_ID = 90;
    final public static String ASSIGNMENT_VERSION = "append-closest-representative-v1";

    final private static String ELIGIBLE_CHAIN =
        "ac.source_id=2 and c.is_polypeptide=1 and q.is_reject=0 and length(q.seq)>=30";

    private static class Assignment {
        int repChainID, repSeqID;
        long repBindingID;
        double pctIdentity;
        String provenance;
    }

    private static class Chain {
        int chainID, seqID;
        long bindingID;
        double aerospaci;
        String sid;
    }

    private static int baseRelease(Statement stmt, int populationID) throws Exception {
        ResultSet rs = stmt.executeQuery("select base_scop_release_id, status from astral_chain_all_population where id="+populationID);
        if (!rs.next())
            throw new Exception("population not found: "+populationID);
        int releaseID = rs.getInt(1);
        String status = rs.getString(2);
        rs.close();
        if (!"validated".equals(status))
            throw new Exception("population is not validated: "+status);
        return releaseID;
    }

    private static int count(Statement stmt, String query) throws Exception {
        ResultSet rs = stmt.executeQuery(query);
        rs.next();
        int rv = rs.getInt(1);
        rs.close();
        return rv;
    }

    /**
       Add every ready, eligible chain that is not yet a candidate,
       and queue BLAST searches for sequences new to the population.
       Chains of entries that are already obsolete are not added, as
       when the population was first built.
    */
    final public static void prepare(int populationID) throws Exception {
        Statement stmt = LocalSQL.createStatement();
        int releaseID = baseRelease(stmt, populationID);

        String newBindings =
            "from astral_chain_stable_binding b " +
            "join astral_chain ac on ac.id=b.astral_chain_id " +
            "join astral_seq q on q.id=ac.seq_id " +
            "join pdb_chain c on c.id=b.pdb_chain_id " +
            "join pdb_entry e on e.id=b.pdb_entry_id " +
            "where b.base_scop_release_id="+releaseID+" and b.status='ready' and "+ELIGIBLE_CHAIN+
            " and e.obsolete_date is null and not exists (select 1 from astral_chain_all_candidate x " +
            "where x.population_id="+populationID+" and x.stable_binding_id=b.id)";

        int noScore = count(stmt, "select count(*) "+newBindings+
                            " and not exists (select 1 from aerospaci a where a.pdb_entry_id=b.pdb_entry_id and a.release_id=b.base_scop_release_id)");
        if (noScore != 0)
            throw new Exception(noScore+" ready chains have no AEROSPACI row for release "+releaseID);

        int added = stmt.executeUpdate(
            "insert into astral_chain_all_candidate " +
            "(population_id,stable_binding_id,astral_chain_id,sid,seq_id,pdb_entry_id,pdb_release_id,origin,aerospaci) " +
            "select "+populationID+",b.id,b.astral_chain_id,b.sid,ac.seq_id,b.pdb_entry_id,b.pdb_release_id,b.origin," +
            "(select a.aerospaci from aerospaci a where a.pdb_entry_id=b.pdb_entry_id and a.release_id=b.base_scop_release_id) " +
            newBindings);
        stmt.executeUpdate("update astral_chain_all_population set candidate_count=(select count(*) from astral_chain_all_candidate where population_id="+populationID+"), sequence_count=(select count(distinct seq_id) from astral_chain_all_candidate where population_id="+populationID+") where id="+populationID);
        stmt.close();

        int queued = AllChainBlast.queueAdded(populationID);
        System.out.println("population "+populationID+": added "+added+" candidate chains; queued "+queued+" BLAST searches");
    }

    /**
       Give every candidate that has no cluster yet its cluster, and
       queue FATCAT jobs for the new representatives.  With dryRun,
       report what would be done and change nothing.
    */
    final public static void assign(int populationID, boolean dryRun) throws Exception {
        Statement stmt = LocalSQL.createStatement();
        baseRelease(stmt, populationID);

        int unfinished = count(stmt, "select count(*) from astral_chain_all_blast_coverage where population_id="+populationID+" and status<>'complete'");
        if (unfinished != 0)
            throw new Exception(unfinished+" BLAST searches of population "+populationID+" are not complete");
        int unsearched = count(stmt, "select count(distinct c.seq_id) from astral_chain_all_candidate c where c.population_id="+populationID+
                               " and not exists (select 1 from astral_chain_all_blast_coverage v where v.population_id=c.population_id and v.query_seq_id=c.seq_id)");
        if (unsearched != 0)
            throw new Exception(unsearched+" sequences of population "+populationID+" have not been searched; run prepare");

        // chains to place, best first: the order in which representatives are chosen
        ArrayList<Chain> added = new ArrayList<Chain>();
        ResultSet rs = stmt.executeQuery(
            "select c.astral_chain_id, c.seq_id, c.stable_binding_id, c.aerospaci, c.sid from astral_chain_all_candidate c " +
            "left join astral_chain_all_cluster_member m on m.population_id=c.population_id and m.pct_identical="+PCT_ID+
            " and m.member_stable_binding_id=c.stable_binding_id where c.population_id="+populationID+
            " and m.population_id is null order by round(c.aerospaci,2) desc, c.sid collate latin1_bin asc");
        while (rs.next()) {
            Chain c = new Chain();
            c.chainID = rs.getInt(1);
            c.seqID = rs.getInt(2);
            c.bindingID = rs.getLong(3);
            c.aerospaci = rs.getDouble(4);
            c.sid = rs.getString(5);
            added.add(c);
        }
        rs.close();
        if (added.size() == 0) {
            System.out.println("population "+populationID+": no chains to assign");
            stmt.close();
            return;
        }

        // cluster of every sequence already placed
        HashMap<Integer,Assignment> bySeq = new HashMap<Integer,Assignment>();
        rs = stmt.executeQuery(
            "select distinct member_seq_id, representative_astral_chain_id, representative_seq_id, representative_stable_binding_id, " +
            "representative_pct_identity, identity_provenance from astral_chain_all_cluster_member where population_id="+populationID+
            " and pct_identical="+PCT_ID);
        while (rs.next()) {
            Assignment a = new Assignment();
            a.repChainID = rs.getInt(2);
            a.repSeqID = rs.getInt(3);
            a.repBindingID = rs.getLong(4);
            a.pctIdentity = rs.getDouble(5);
            a.provenance = rs.getString(6);
            if (bySeq.put(rs.getInt(1), a) != null)
                throw new Exception("sequence "+rs.getInt(1)+" is in more than one cluster");
        }
        rs.close();

        // representatives, by sequence
        HashMap<Integer,Chain> repBySeq = new HashMap<Integer,Chain>();
        rs = stmt.executeQuery(
            "select c.astral_chain_id, c.seq_id, c.stable_binding_id, c.aerospaci, c.sid from astral_chain_all_subset_id s " +
            "join astral_chain_all_candidate c on c.population_id=s.population_id and c.astral_chain_id=s.astral_chain_id " +
            "where s.population_id="+populationID+" and s.pct_identical="+PCT_ID);
        while (rs.next()) {
            Chain c = new Chain();
            c.chainID = rs.getInt(1);
            c.seqID = rs.getInt(2);
            c.bindingID = rs.getLong(3);
            c.aerospaci = rs.getDouble(4);
            c.sid = rs.getString(5);
            repBySeq.put(c.seqID, c);
        }
        rs.close();
        int rank = count(stmt, "select coalesce(max(representative_rank),0) from astral_chain_all_subset_id where population_id="+populationID+" and pct_identical="+PCT_ID);

        // identity of a sequence to every other, from searches in either direction
        PreparedStatement edges = LocalSQL.prepareStatement(
            "select seq2_id, max(pct_identical) from astral_chain_all_blast_edge where population_id=? and seq1_id=? and pct_identical>=? group by seq2_id " +
            "union all select seq1_id, max(pct_identical) from astral_chain_all_blast_edge where population_id=? and seq2_id=? and pct_identical>=? group by seq1_id");

        ArrayList<Chain> newReps = new ArrayList<Chain>();
        HashMap<Chain,Assignment> placed = new LinkedHashMap<Chain,Assignment>();
        int sameSeq = 0, joined = 0;
        for (Chain c : added) {
            Assignment a = bySeq.get(c.seqID);
            if (a != null)
                sameSeq++;
            else {
                Chain best = null;
                double bestPct = 0.0;
                edges.setInt(1, populationID); edges.setInt(2, c.seqID); edges.setInt(3, PCT_ID);
                edges.setInt(4, populationID); edges.setInt(5, c.seqID); edges.setInt(6, PCT_ID);
                rs = edges.executeQuery();
                while (rs.next()) {
                    Chain rep = repBySeq.get(rs.getInt(1));
                    if (rep == null)
                        continue;
                    double pct = rs.getDouble(2);
                    // closest representative; ties as in MakeAllChainClusters
                    int cmp = (best == null) ? 1 : Double.compare(pct, bestPct);
                    if (cmp == 0)
                        cmp = Double.compare(Math.round(rep.aerospaci*100.0), Math.round(best.aerospaci*100.0));
                    if (cmp == 0)
                        cmp = best.sid.compareTo(rep.sid);
                    if (cmp == 0)
                        cmp = Integer.compare(best.chainID, rep.chainID);
                    if (cmp > 0) {
                        best = rep;
                        bestPct = pct;
                    }
                }
                rs.close();

                a = new Assignment();
                if (best == null) {
                    // matches no representative: starts a new cluster
                    repBySeq.put(c.seqID, c);
                    newReps.add(c);
                    a.repChainID = c.chainID;
                    a.repSeqID = c.seqID;
                    a.repBindingID = c.bindingID;
                    a.pctIdentity = 100.0;
                    a.provenance = "same_seq";
                }
                else {
                    joined++;
                    a.repChainID = best.chainID;
                    a.repSeqID = best.seqID;
                    a.repBindingID = best.bindingID;
                    a.pctIdentity = Math.min(bestPct, 100.0);
                    a.provenance = "blastpgp-2.2.18";
                }
                bySeq.put(c.seqID, a);
            }
            placed.put(c, a);
        }
        edges.close();

        System.out.println("population "+populationID+": "+added.size()+" chains to assign; "+sameSeq+
                           " with a sequence already placed, "+joined+" joining a cluster by BLAST, "+
                           newReps.size()+" new representatives");
        if (dryRun) {
            for (Chain c : newReps)
                System.out.println("  new representative "+c.sid);
            stmt.close();
            return;
        }

        LocalSQL.setAutoCommit(false);
        try {
            PreparedStatement addRep = LocalSQL.prepareStatement(
                "insert into astral_chain_all_subset_id (population_id,pct_identical,astral_chain_id,representative_rank) values (?,?,?,?)");
            for (Chain c : newReps) {
                addRep.setInt(1, populationID); addRep.setInt(2, PCT_ID);
                addRep.setInt(3, c.chainID); addRep.setInt(4, ++rank);
                addRep.addBatch();
            }
            addRep.executeBatch();
            addRep.close();

            PreparedStatement addMember = LocalSQL.prepareStatement(
                "insert into astral_chain_all_cluster_member (population_id,pct_identical,member_stable_binding_id,member_astral_chain_id," +
                "representative_stable_binding_id,representative_astral_chain_id,member_seq_id,representative_seq_id," +
                "representative_pct_identity,identity_provenance,assignment_version,created_at) values (?,?,?,?,?,?,?,?,?,?,?,now())");
            for (Map.Entry<Chain,Assignment> e : placed.entrySet()) {
                Chain c = e.getKey();
                Assignment a = e.getValue();
                addMember.setInt(1, populationID); addMember.setInt(2, PCT_ID);
                addMember.setLong(3, c.bindingID); addMember.setInt(4, c.chainID);
                addMember.setLong(5, a.repBindingID); addMember.setInt(6, a.repChainID);
                addMember.setInt(7, c.seqID); addMember.setInt(8, a.repSeqID);
                addMember.setDouble(9, a.pctIdentity); addMember.setString(10, a.provenance);
                addMember.setString(11, ASSIGNMENT_VERSION);
                addMember.addBatch();
            }
            addMember.executeBatch();
            addMember.close();

            stmt.executeUpdate("update astral_chain_all_population set representative_count=(select count(*) from astral_chain_all_subset_id where population_id="+populationID+" and pct_identical="+PCT_ID+") where id="+populationID);
            LocalSQL.commit();
        }
        catch (Exception e) {
            LocalSQL.rollback();
            throw e;
        }
        finally {
            LocalSQL.setAutoCommit(true);
        }

        // structure comparisons for the new clusters; queued after the
        // commit, because LocalSQL.newJob locks tables
        if (newReps.size() > 0) {
            int runID = ChainDomainFatcatPairs.createRun(populationID);
            int queued = 0;
            for (Chain c : newReps)
                for (int segment=0; segment<2; segment++)
                    if (LocalSQL.newJob(25, c.chainID, segment+" "+populationID+" "+runID, stmt) != 0)
                        queued++;
            System.out.println("queued "+queued+" FATCAT jobs in run "+runID);
        }
        stmt.close();
    }

    final public static void main(String argv[]) {
        try {
            boolean dryRun = (argv.length == 3) && argv[2].equals("--dry-run");
            if (!((argv.length == 2) || dryRun) || !(argv[0].equals("prepare") || argv[0].equals("assign")) ||
                (dryRun && !argv[0].equals("assign")))
                throw new Exception("Usage: ExtendAllChainPopulation prepare population_id | assign population_id [--dry-run]");
            LocalSQL.connectRW();
            int populationID = Integer.parseInt(argv[1]);
            if (argv[0].equals("prepare"))
                prepare(populationID);
            else
                assign(populationID, dryRun);
        }
        catch (Exception e) {
            System.err.println("Exception: "+e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
}
