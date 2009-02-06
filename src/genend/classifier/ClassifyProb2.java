package genend.classifier;

import genend.util.KmerDistribProcessor;
import genend.util.SeqUtils;
import genend.util.container.ClassifyResultObj;
import genend.util.config.BioSQLConnectionConfig;
import genend.util.config.BioSQLConfigException;

import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.Date;
import java.text.DateFormat;
import java.util.concurrent.*;


//import org.biojava.bio.seq.db.biosql.BioSQLSequenceDB;
//import org.biojava.bio.seq.Sequence;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojavax.bio.db.biosql.BioSQLRichSequenceDB;
import org.biojavax.bio.seq.RichSequence;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.*;

public class ClassifyProb2
{
    private int kmer_min, kmer_max, num_threads;
    private String output_dir, input_dir, conf_name;
    private boolean gen_distrib;
    private int powers_of_four[] = new int[20];
    private int[] piece_sizes;
    private KmerDistribProcessor kmer_distrib_proc = null;
    private final int iterations = 50;

    public ClassifyProb2(int[] piece_sizes, int kmer_min, int kmer_max, String input_dir,
                         String output_dir, boolean gen_distrib, int num_threads,
                         String conf_name)
    {
        this.piece_sizes = piece_sizes;
        this.kmer_min = kmer_min;
        this.kmer_max = kmer_max;
        this.input_dir = input_dir;
        this.output_dir = output_dir;
        this.gen_distrib = gen_distrib;
        this.num_threads = num_threads;
        this.conf_name = conf_name;
        kmer_distrib_proc = new KmerDistribProcessor(input_dir, output_dir, kmer_min, kmer_max, num_threads);

        for(int i = 0; i < powers_of_four.length; i++)
            powers_of_four[i] = (int)Math.pow(4.0, (double)i);
    }

    public void execute()
    {
        ThreadPoolExecutor tpe = null;
        Vector<ClassifyResultObj> ret_vector = new Vector<ClassifyResultObj>();

        if (gen_distrib)
            kmer_distrib_proc.processFiles();

        File input_dir_h = new File(input_dir);
        File[] input_files = input_dir_h.listFiles();

        for (int piece_size: piece_sizes)
        {
            String str_date = DateFormat.getDateInstance(DateFormat.SHORT).format(new Date()).replace("/", ".");
            String filename = "phylo-" + Integer.toString(piece_size) + "-" + str_date + ".dat";
            FileOutputStream fout = null;
            try { fout = new FileOutputStream(filename); }
            catch (Exception e) { e.printStackTrace(); }

            for (int kmer_size = kmer_min; kmer_size <= kmer_max; kmer_size++)
            {
                tpe = new ThreadPoolExecutor(num_threads, num_threads,
                                             25000L, TimeUnit.MILLISECONDS,
                                             new LinkedBlockingQueue<Runnable>());

                RelDistReader dist_rdr = new RelDistReader(output_dir, kmer_size);
                HashMap<String, HashMap<String, Double>> rel_distribs = dist_rdr.get();

                for (int i = 0; i < input_files.length; i++)
                {
                    System.out.println("conf_name: " + conf_name);
                    Runnable r = new StatRun(iterations, input_files[i].toString(),
                                             piece_size, kmer_size, rel_distribs,
                                             ret_vector, tpe, conf_name);

                    tpe.execute(r);
                }

                tpe.shutdown();
                try { tpe.awaitTermination(10000000L, TimeUnit.SECONDS); }
                catch (Exception e) { e.printStackTrace(); }

                System.out.println("\tKMER: " + kmer_size + "... done.");
                if (ret_vector.size() > 0)
                {
                    double[] vals = processKmer(ret_vector, (double)(iterations*ret_vector.size()));
                    String write_ln = Integer.toString(kmer_size) + "\t" +
                        Double.toString(vals[0]) + "\t" + Double.toString(vals[1]) +
                        "\t" + Double.toString(vals[2]) + "\t" + Double.toString(vals[3]) +
                        "\t" + Double.toString(vals[4]) + "\t" + Double.toString(vals[5]) +
                        "\t" + Double.toString(vals[6]);
                    try { new PrintStream(fout).println(write_ln); }
                    catch (Exception e) { e.printStackTrace(); }
                }
                else System.out.println("No results in result vector.\n");
            }

            System.out.println("PIECE: " + piece_size + "... done.");
        }
    }

    double[] processKmer(Vector<ClassifyResultObj> results, double divisor)
    {
        int size = results.size();
        int piece_size = results.get(0).getPieceSize();
        int kmer_size = results.get(0).getKmerSize();
        double[] match_val = new double[7];

        for (int i = 0; i < size; i++)
        {
            ClassifyResultObj item = results.get(i);
            match_val[0] += (double)item.getNameMatchVal();
            match_val[1] += (double)item.getKingdomMatchVal();
            match_val[2] += (double)item.getPhylumMatchVal();
            match_val[3] += (double)item.getClassMatchVal();
            match_val[4] += (double)item.getOrderMatchVal();
            match_val[5] += (double)item.getFamilyMatchVal();
            match_val[6] += (double)item.getGenusMatchVal();
        }

        for (int i = 0; i < match_val.length; i++)
            match_val[i] /= divisor;
        results.clear();
        return match_val;
    }

    String[] parse_genome_file(String path)
    {
        String ret[] = new String[2];
        BufferedReader br = null;
        StringBuilder str_builder = new StringBuilder();
        try
        {
            br = new BufferedReader(new FileReader(path));
            //ret[0] = parse_header(br.readLine());
            br.readLine();
            ret[0] = (new File(path)).getName();
            String line = "";
            while ((line = br.readLine()) != null)
                 str_builder.append(line.toUpperCase());
            ret[1] = str_builder.toString();
            br.close();
        }
        catch (Exception e) { }
        finally { try { if (br != null) br.close(); } catch (Exception e) {} }

        return ret;
    }

    String parse_header(String header)
    {
        String[] tokens = header.split(" ");
        String ret_str = tokens[1];
        int count = 2;
        while (tokens[count].endsWith(",") == false)
        {
            ret_str = ret_str + "_" + tokens[count];
            count++;
        }
        ret_str = ret_str + "_" + tokens[count];
        if (ret_str.contains(","))
            ret_str = ret_str.substring(0, ret_str.length() - 1);
        return ret_str;
    }

    private class StatRun extends Thread
    {
        private int total_frags, iterations, piece_size, kmer_size;

        private int[] match_vals = new int[7];
        private Vector<ClassifyResultObj> ret_vector;
        private String input_file, spec_name, conf_name;
        private HashMap<String, HashMap<String, Double>> rel_distribs;
        private Random rg;
        private ThreadPoolExecutor tpe;
        private BioSQLConnectionConfig biosql = null;

        StatRun(int iterations, String input_file, int piece_size,
                int kmer_size, HashMap<String, HashMap<String, Double>> rel_distribs,
                Vector<ClassifyResultObj> ret_vector, ThreadPoolExecutor tpe, String conf_name)
        {
            this.iterations = iterations;
            this.kmer_size = kmer_size;
            this.piece_size = piece_size;
            this.input_file = input_file;
            this.rel_distribs = rel_distribs;
            this.ret_vector = ret_vector;
            this.tpe = tpe;
            rg = new Random();
            try {
                biosql = new BioSQLConnectionConfig(conf_name, "default");
            }
            catch (BioSQLConfigException e) { e.printStackTrace(); }
        }

        public void run()
        {
            String[] cur_org_info = parse_genome_file(input_file);
            spec_name = cur_org_info[0];
            String org_seq = cur_org_info[1];
            int org_len = org_seq.length();
            double B, Bl, jp_factor, jp_summant;

            System.out.println("\tSPECIES: "+spec_name+
                               " | KMER: "+kmer_size+" | PIECE: "+piece_size+
                               " | THREADS: "+tpe.getActiveCount()+
                               " | TASKS LEFT: "+
                               (tpe.getTaskCount()-tpe.getCompletedTaskCount()));

            if (org_len < 500000)
                return;

            B = powers_of_four[kmer_size - 1];
            Bl = B / 2.0;
            jp_factor = (double)org_len / ((double)org_len + Bl);
            jp_summant = 1.0 / (2.0 * ((double)org_len + Bl));

            File f = new File(output_dir);
            File[] files = f.listFiles();

            for (int i = 0; i < iterations; i++)
            {
                HashMap<String, ArrayList<Double>> prob_set =
                    getEmptyProbSet(rel_distribs.keySet());

                String str_kmer = getFragment(org_seq, org_len, piece_size);
                int count = 0, str_kmer_len;

                while (!SeqUtils.matchBases(str_kmer) && count < 1000)
                {
                    str_kmer = getFragment(org_seq, org_len, piece_size);
                    count++;
                }
                str_kmer_len = str_kmer.length();

                for (int ind = 1; ind < str_kmer_len - kmer_size + 1; ind++)
                {
                    String tmp_frag = str_kmer.substring(ind, ind + kmer_size);

                    for (int j = 0; j < files.length; j++)
                    {
                        if ((files[j].getName()).endsWith("-"+Integer.toString(kmer_size)))
                        {
                            String org_str = files[j].getName();
                            org_str = org_str.substring(0, (int)(org_str.length() - 2));

                            double log_tmp, tmp_val;
                            log_tmp =
                                Math.log(SeqUtils.jp(rel_distribs.get(org_str).get(tmp_frag),
                                                     jp_factor, jp_summant));
                            tmp_val =
                                ((Double)prob_set.get(org_str).get(ind - 1)).doubleValue() + log_tmp;
                            prob_set.get(org_str).add(tmp_val);
                        }
                    }
                }

                String high_spec_name = getHighestProb(prob_set);

                if (high_spec_name.equals(spec_name))
                {
                    for (int e = 0; e < 8; e++)
                        match_vals[e]++;
                }
                else match_vals = compareTaxonomies(spec_name, high_spec_name);
            }

            ret_vector.add(new ClassifyResultObj(kmer_size, piece_size, match_vals[0],
                                                 match_vals[1], match_vals[2], match_vals[3],
                                                 match_vals[4], match_vals[5], match_vals[6]));
        }

        String[] fetchTaxonomy(String id)
        {
            String[] tax = new String[7];
            String url = "jdbc:mysql://" + biosql.getServer() + ":3306/biosql";
            int index = 0;
            // load the mysql driver
            try { Class.forName("com.mysql.jdbc.Driver"); }
            catch (ClassNotFoundException ex) {
                System.out.println("Cannot find DB driver, is it on your classpath?");
            }

            try {
                Connection conn = DriverManager.getConnection(url,
                        biosql.getUsername(), biosql.getPassword());

                String cur_taxon_id = null;
                String query = "SELECT taxon_id FROM bioentry WHERE name='" + id + "'";
                Statement stmt = conn.createStatement();
                ResultSet rs = stmt.executeQuery(query);

                if (rs.next()) cur_taxon_id = rs.getString(1);
                else return null;

                boolean done = false;

                while (done == false)
                {
                    String rank = "", next_taxon_id;
                    query = "SELECT node_rank, parent_taxon_id " +
                            "FROM taxon WHERE taxon_id='" + cur_taxon_id +"'";

                    rs = stmt.executeQuery(query);
                    if (rs.next()) {
                        
                        rank = rs.getString(1);

                        if (rank.equals("genus") &&
                            rank.equals("species") &&
                            rank.equals("order") &&
                            rank.equals("family") &&
                            rank.equals("phylum") &&
                            rank.equals("class")) continue;

                        if (rank.equals("kingdom")) break;

                        next_taxon_id = rs.getString(2);

                        query = "SELECT name FROM taxon_name WHERE " +
                                "name_class='scientific name' AND " +
                                "taxon_id='" + cur_taxon_id;

                        stmt.execute(query);

                        if (rs.next())
                        {
                            tax[index] = rs.getString(1);
                        }
                        else return null;

                        cur_taxon_id = next_taxon_id;
                    }
                    else return null;
                }

                rs.close();
                stmt.close();
                conn.close();

            } catch (SQLException ex) {
                Logger.getLogger(ClassifyProb2.class.getName()).log(Level.SEVERE, null, ex);
            }

            return tax;
        }

        int[] compareTaxonomies(String cur_spec, String high_spec)
        {
            int[] match_vals = new int[7];
            String[] cur_spec_tax, high_spec_tax;

            cur_spec_tax = fetchTaxonomy(cur_spec);
            high_spec_tax = fetchTaxonomy(high_spec);

            for (int i = 6; i >= 0; i--)
            {
                if (cur_spec_tax[i].equals(high_spec_tax[i]))
                    match_vals[(i-6)*-1]++;
            }

            return match_vals;
        }

        String getHighestProb(HashMap<String, ArrayList<Double>> prob_set)
        {
            Map.Entry tmp_pair = null, high_pair = null;
            Iterator itr = prob_set.entrySet().iterator();
            high_pair = (Map.Entry) itr.next();
            while (itr.hasNext())
            {
                tmp_pair = (Map.Entry) itr.next();
                ArrayList<Double> high_arr = (ArrayList<Double>)high_pair.getValue();
                ArrayList<Double> tmp_arr = (ArrayList<Double>)tmp_pair.getValue();
                Double tmp_dbl = tmp_arr.get(tmp_arr.size() - 1);
                Double high_dbl = high_arr.get(high_arr.size() - 1);
                //                System.out.println("tmp_dbl: "+tmp_dbl + " | high_dbl: "+high_dbl);
                if (tmp_dbl.compareTo(high_dbl) >= 0)
                    high_pair = tmp_pair;
            }
            return (String)high_pair.getKey();
        }

        HashMap<String, ArrayList<Double>> getEmptyProbSet(Set<String> spec_set)
        {
            HashMap<String, ArrayList<Double>> ret_set =
                new HashMap<String, ArrayList<Double>>();

            Iterator itr = spec_set.iterator();

            while(itr.hasNext())
            {
                String next = (String)itr.next();

                ArrayList<Double> list = new ArrayList<Double>();
                list.add(new Double(0));
                ret_set.put(next, list);
            }

            return ret_set;
        }

        String getFragment(String seq_str, int seq_len, int frag_size)
        {
            int frags_per_piece = piece_size - frag_size + 1;
            int local_pos = total_frags % frags_per_piece;
            int piece_pos = 0;

            if (local_pos == 0)
                piece_pos = rg.nextInt(seq_len - piece_size - frag_size + 1) + 1;

            total_frags++;
            return seq_str.substring(piece_pos + local_pos,
                                     piece_pos + local_pos + frag_size);
        }

    }

    private class RelDistReader
    {
        String input_dir;
        int kmer_size;

        public RelDistReader(String input_dir, int kmer_size)
        {
            this.input_dir = input_dir;
            this.kmer_size = kmer_size;
        }

        private Vector<String> getStringVector()
        {
            Vector<String> ret_vect = new Vector<String>();
            for (int i = 0; i < powers_of_four[kmer_size]; i++)
                ret_vect.add(SeqUtils.numToStr(i, kmer_size));
            return ret_vect;
        }

        public HashMap<String, HashMap<String, Double>> get()
        {
            File input_dir_h = new File(input_dir);
            File files[] = input_dir_h.listFiles();
            HashMap<String, HashMap<String, Double>> ret_table =
                new HashMap<String, HashMap<String, Double>>();
            Double zero = new Double(0.0);
            int count = 0;

            Vector<String> kmer_strings = getStringVector();

            for (int i = 0; i < files.length; i++)
            {
                if ((files[i].getName()).endsWith("-"+Integer.toString(kmer_size)))
                {
                    String spec_name = files[i].getName();
                    spec_name = spec_name.substring(0, (int)(spec_name.length() - 2));
                    ret_table.put(spec_name, new HashMap<String, Double>());
                    HashMap<String, Double> cur_spec_tbl = ret_table.get(spec_name);
                    try
                    {
                        BufferedReader br = new BufferedReader(new FileReader(files[i]));
                        String line = null;
                        int ind = 0;
                        while((line = br.readLine()) != null)
                        {
                            String[] tokens = line.split("\t");
                            String cur_string = kmer_strings.get(ind++);
                            if (Double.valueOf(tokens[1]) == 0)
                                cur_spec_tbl.put(cur_string, zero);
                            else
                                cur_spec_tbl.put(cur_string, Double.valueOf(tokens[1]));
                        }
                    }
                    catch (Exception e) { e.printStackTrace(); }
                }
            }
            System.out.println("Done loading distributions for K-mer size: " + kmer_size);
            return ret_table;
        }
    }
}
