package genend.classifier;
import genend.util.SeqUtils;
import genend.util.container.ResultObj;

import java.io.*;
import java.util.*;
import java.text.DateFormat;
import java.util.concurrent.*;



public class OverlapBayesianClassifier
{
    private int kmer_min, kmer_max, num_threads;
    private String output_dir, input_dir;
    private boolean gen_distrib;
    private int powers_of_four[] = new int[20];
    private int[] piece_sizes;
    private final int offset = 5000;


    public OverlapBayesianClassifier(int[] piece_sizes, int kmer_min, int kmer_max, String input_dir,
                     String output_dir, boolean gen_distrib, int num_threads)
    {
        this.piece_sizes = piece_sizes;
        this.kmer_min = kmer_min;
        this.kmer_max = kmer_max;
        this.input_dir = input_dir;
        this.output_dir = output_dir;
        this.gen_distrib = gen_distrib;
        this.num_threads = num_threads;

        for(int i = 0; i < powers_of_four.length; i++)
            powers_of_four[i] = (int)Math.pow(4.0, (double)i);
    }

    public void execute()
    {
        ThreadPoolExecutor tpe = null;
        Vector<ResultObj> ret_vector = new Vector<ResultObj>();

        File input_dir_h = new File(input_dir);
        File[] input_files = input_dir_h.listFiles();

        for (int piece_size: piece_sizes)
        {
            String str_date = DateFormat.getDateInstance(DateFormat.SHORT).format(new Date()).replace("/", ".");
            String filename = "piece-" + Integer.toString(piece_size) + "-" + str_date + ".dat";
            String match_filename = "match-" + Integer.toString(piece_size) + "-" + str_date + ".dat";
            FileOutputStream fout = null, match_fout = null;
            try {
                fout = new FileOutputStream(filename);
                match_fout = new FileOutputStream(match_filename);
            }
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
                    Runnable r = new StatRun(offset, input_files[i].toString(),
                                             piece_size, kmer_size, rel_distribs,
                                             ret_vector, tpe);

                    tpe.execute(r);
                }

                tpe.shutdown();
                try { tpe.awaitTermination(10000000L, TimeUnit.SECONDS); }
                catch (Exception e) { e.printStackTrace(); }

                System.out.println("\tKMER: " + kmer_size + "... done.");

                // output matches to file
                PrintStream printer;
                try {
                    printer = new PrintStream(match_fout);
                    for (int f = 0; f < ret_vector.size(); f++) {
                        ResultObj tmp_obj = ret_vector.get(f);

                        printer.println(kmer_size + "\t" + tmp_obj.getCurSpec()
                                + "\t" + tmp_obj.getHighSpec() + "\t" + tmp_obj.getScore());
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }


                double val = processKmer(ret_vector, (double)(ret_vector.size()));
                String write_ln = Integer.toString(kmer_size) + "\t" + Double.toString(val);
                try { new PrintStream(fout).println(write_ln); }
                catch (Exception e) { e.printStackTrace(); }
            }

            System.out.println("PIECE: " + piece_size + "... done.");
        }
    }

    double processKmer(Vector<ResultObj> results, double divisor)
    {
        int size = results.size();
        int piece_size = results.get(0).getPieceSize();
        int kmer_size = results.get(0).getKmerSize();
        double match_val = 0.0;

        for (int i = 0; i < size; i++)
        {
            ResultObj item = results.get(i);
            match_val += (double)item.getMatchVal();
        }

        results.clear();
        return match_val / divisor;
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
            while ((line = br.readLine()) != null) {
                if (line.contains(">")) break;
                str_builder.append(line.toUpperCase());
            }
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
        private int total_frags, offset, piece_size, kmer_size, org_matches;
        private Vector<ResultObj> ret_vector;
        private String input_file, spec_name;
        private HashMap<String, HashMap<String, Double>> rel_distribs;
        private Random rg;
        private ThreadPoolExecutor tpe;

        StatRun(int offset, String input_file, int piece_size,
                int kmer_size, HashMap<String, HashMap<String, Double>> rel_distribs,
                Vector<ResultObj> ret_vector, ThreadPoolExecutor tpe)
        {
            this.offset = offset;
            this.kmer_size = kmer_size;
            this.piece_size = piece_size;
            this.input_file = input_file;
            this.rel_distribs = rel_distribs;
            this.ret_vector = ret_vector;
            this.tpe = tpe;
            org_matches = 0;
            rg = new Random();
        }

        @Override
        public void run()
        {
            String[] cur_org_info = parse_genome_file(input_file);
            spec_name = cur_org_info[0];
            ProbObj high_spec_obj = null;
            String high_spec_name = "";
            String org_seq = cur_org_info[1];
            int org_len = org_seq.length();
            double B, Bl, jp_factor, jp_summant;

            System.out.println("\tSPECIES: "+spec_name.replace('_',' ')+
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

            for (int i = 0; i < org_len; org_len += offset)
            {
                HashMap<String, ArrayList<Double>> prob_set =
                    getEmptyProbSet(rel_distribs.keySet());

                String str_kmer = "";
                int str_kmer_len;

                if (i + piece_size < org_len)
                    str_kmer = org_seq.substring(i, i+piece_size);
                else
                    str_kmer = org_seq.substring(i);

                str_kmer_len = str_kmer.length();

                for (int ind = 1; ind < str_kmer_len - kmer_size + 1; ind++)
                {
                    String tmp_frag = str_kmer.substring(ind, ind + kmer_size);
                    if (!SeqUtils.matchBases(tmp_frag)) continue;

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

                high_spec_obj = getHighestProb(prob_set);
                high_spec_name = high_spec_obj.getName();

                if (high_spec_name.equals(spec_name))
                    org_matches++;
                
                ret_vector.add(new ResultObj(kmer_size, piece_size, org_matches,
                    spec_name, high_spec_name, high_spec_obj.getScore()));
            }

            // moved ret_vector above from here to where it is now.
        }

        ProbObj getHighestProb(HashMap<String, ArrayList<Double>> prob_set)
        {
            Map.Entry tmp_pair = null, high_pair = null;
            Iterator itr = prob_set.entrySet().iterator();
            Double high_dbl = null, tmp_dbl = null;
            high_pair = (Map.Entry) itr.next();
            while (itr.hasNext())
            {
                tmp_pair = (Map.Entry) itr.next();
                ArrayList<Double> high_arr = (ArrayList<Double>)high_pair.getValue();
                ArrayList<Double> tmp_arr = (ArrayList<Double>)tmp_pair.getValue();
                tmp_dbl = tmp_arr.get(tmp_arr.size() - 1);
                high_dbl = high_arr.get(high_arr.size() - 1);
                //                System.out.println("tmp_dbl: "+tmp_dbl + " | high_dbl: "+high_dbl);
                if (tmp_dbl.compareTo(high_dbl) >= 0)
                    high_pair = tmp_pair;
            }
            return new ProbObj((String)high_pair.getKey(),
                    high_dbl.doubleValue());
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
            System.out.println("input_dir: "+input_dir);
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
                    /*
                    for (int k = 0; k < powers_of_four[kmer_size]; k++)
                    {
                        String str = SeqUtils.numToStr(k, kmer_size);
                        if (cur_spec_tbl.containsKey(str) == false)
                            cur_spec_tbl.put(str, zero);
                    }
                    */
                }
            }
            System.out.println("Done loading distributions for K-mer size: " + kmer_size);
            return ret_table;
        }
    }

    private class ProbObj {

        private String name;
        private double score;

        public ProbObj(String name, double score) {
            this.name = name;
            this.score = score;
        }

        public String getName() {
            return name;
        }

        public double getScore() {
            return score;
        }
    }

}
