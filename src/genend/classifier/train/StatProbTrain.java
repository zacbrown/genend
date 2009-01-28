package genend.classifier.train;
import genend.util.SeqUtils;
import genend.util.container.ResultKeeper;
import genend.util.container.ResultObj;

import java.io.*;
import java.util.*;
import java.text.DateFormat;
import java.util.concurrent.*;

public class StatProbTrain
{
    private int kmer_low, kmer_high, num_threads;
    private double perc_train, divisor;
    private String output_dir, input_dir, genome_dir;
    private int powers_of_four[] = new int[20];
    private final int iterations = 1000;
    private boolean debug = false;
    private ResultKeeper my_keeper = null;

    public StatProbTrain(double perc_train, int kmer_low, int kmer_high, double divisor,
    		String genome_dir, String input_dir, String output_dir, int num_threads)
    {
        this.kmer_low = kmer_low;
        this.kmer_high = kmer_high;
        this.input_dir = input_dir;
        this.genome_dir = genome_dir;
        this.output_dir = output_dir;
        this.num_threads = num_threads;
        this.perc_train = perc_train;
        this.divisor = divisor;
        this.my_keeper = new ResultKeeper(kmer_low, kmer_high);
        
        for(int i = 0; i < powers_of_four.length; i++)
            powers_of_four[i] = (int)Math.pow(4.0, (double)i);
    }

    public void setDebug() { debug = !debug; }
    public boolean getDebug() { return debug; }
    
    public void execute()
    {
        ThreadPoolExecutor tpe = null;
        HashMap<String, Integer> org_len_map = getOrgLenMap(); 

        File input_dir_h = new File(input_dir);
        FilenameFilter filter = new FilenameFilter()
        {
        	public boolean accept(File dir, String name)
        	{
        		return name.startsWith("FRAG");
        	}
        };
        File[] input_filenames = input_dir_h.listFiles(filter);
        
        String str_date = DateFormat.getDateInstance(DateFormat.SHORT).format(new Date()).replace("/", ".");
        String filename = "train-" + (int)(perc_train*100) + "-" + 
        	(new File(genome_dir)).listFiles().length + str_date + ".dat";
        BufferedWriter bw = null;
        try { bw = new BufferedWriter(new FileWriter(filename)); }
        catch (Exception e) { e.printStackTrace(); }

        for (int kmer_size = kmer_low; kmer_size <= kmer_high; kmer_size++)
        {
        	tpe = new ThreadPoolExecutor(num_threads, num_threads,
        			25000L, TimeUnit.MILLISECONDS,
        			new LinkedBlockingQueue<Runnable>());

        	RelDistReader dist_rdr = new RelDistReader(output_dir, kmer_size);
        	HashMap<String, HashMap<String, Double>> rel_distribs = dist_rdr.get();

        	for (int i = 0; i < input_filenames.length; i++)
        	{
        		File cur_file = input_filenames[i];

        		Runnable r = new StatRun(iterations, cur_file,
        				kmer_size, rel_distribs, my_keeper, tpe, org_len_map);
        		tpe.execute(r);
        	}

        	tpe.shutdown();
        	try { tpe.awaitTermination(10000000L, TimeUnit.SECONDS); }
        	catch (Exception e) { e.printStackTrace(); }

        	System.out.println("\tKMER: " + kmer_size + "... done.");
        	
        	String write_ln = Integer.toString(kmer_size) + "\t" + 
        		Double.toString(my_keeper.get(kmer_size)) + "\n";
        	try { bw.write(write_ln); }
        	catch (Exception e) { e.printStackTrace(); }
        }
        
        try { bw.write(Double.toString(my_keeper.get_div(kmer_low))+"\n"); bw.close(); }
        catch (Exception e) { e.printStackTrace(); }
    }
    
    HashMap<String, Integer> getOrgLenMap()
    {
    	HashMap<String, Integer> ret_map = new HashMap<String, Integer>();
    	File input_dir_h = new File(genome_dir);
    	File[] input_files = input_dir_h.listFiles();
    	
    	for (int i = 0; i < input_files.length; i++)
    	{
    		String[] cur_org = parse_genome_file(input_files[i].toString());
    		ret_map.put(cur_org[0], Integer.valueOf(cur_org[1].length()));
    	}
    	
    	return ret_map;
    }
    
    String[] parse_genome_file(String path)
    {
        String ret[] = new String[2];
        BufferedReader br = null;
        StringBuilder str_builder = new StringBuilder();
        try
        {
            br = new BufferedReader(new FileReader(path));
            ret[0] = parse_header(br.readLine());
            String line = "";
            while ((line = br.readLine()) != null)
                 str_builder.append(line.toUpperCase());
            ret[1] = str_builder.toString();
            br.close();
        }
        catch (Exception e) { e.printStackTrace(); }
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
        private int kmer_size;
        private ResultKeeper my_keeper = null;
        private String input_file, spec_name, short_filename;
        private HashMap<String, HashMap<String, Double>> rel_distribs;
        private HashMap<String, Integer> org_len_map;
        private ThreadPoolExecutor tpe;

        StatRun(int iterations, File input_file, int kmer_size, 
        		HashMap<String, HashMap<String, Double>> rel_distribs,
                ResultKeeper my_keeper, ThreadPoolExecutor tpe,
                HashMap<String, Integer> org_len_map)
        {
            this.kmer_size = kmer_size;
            this.input_file = input_file.toString();
            this.rel_distribs = rel_distribs;
            this.short_filename = input_file.getName();
            this.org_len_map = org_len_map;
            this.my_keeper = my_keeper;
            this.tpe = tpe;
        }

        public void run()
        {
            double B, Bl, jp_factor, jp_summant;
            String[] filename_toks = short_filename.split(",");
            spec_name = filename_toks[2].replace(",", " ");
            spec_name = spec_name.replace("#", "/");
            int org_len;
            if (org_len_map.containsKey(spec_name.replace(" ", "_")))
            	org_len = org_len_map.get(spec_name.replace(" ", "_")).intValue();
            else return;
            if (debug)
            	System.out.println("PROCESSING: "+spec_name+
                               	" | THREADS: "+tpe.getActiveCount()+
                               	" | TASKS LEFT: "+
                               	(tpe.getTaskCount()-tpe.getCompletedTaskCount()));

            B = powers_of_four[kmer_size - 1];
            Bl = B / 2.0;
            jp_factor = (double)org_len / ((double)org_len + Bl);
            jp_summant = 1.0 / (2.0 * ((double)org_len + Bl));

            String line = "";
            StringBuilder str_bldr = new StringBuilder();
            BufferedReader br = null;
            try
            {
            	br = new BufferedReader(new FileReader(input_file));
            	while((line = br.readLine()) != null) str_bldr.append(line);
            	br.close();
            }
            catch(Exception e) { e.printStackTrace(); }
            
            String str_kmer = str_bldr.toString().trim();
            int str_kmer_len = str_kmer.length();
            
            HashMap<String, ArrayList<Double>> prob_set =
            	getEmptyProbSet(rel_distribs.keySet());
            for (int ind = 1; ind < str_kmer_len - kmer_size + 1; ind++)
            {
            	String tmp_frag = str_kmer.substring(ind - 1, ind - 1 + kmer_size);
            	Iterator itr = rel_distribs.entrySet().iterator();

            	while (itr.hasNext())
            	{
            		Map.Entry cur_train_org = (Map.Entry)itr.next();
            		String org_str = (String)cur_train_org.getKey();

            		HashMap<String, Double> cur_org_distribs = 
            			(HashMap<String, Double>)cur_train_org.getValue();

            		double log_tmp, tmp_val;
            		log_tmp =
            			Math.log(SeqUtils.jp(cur_org_distribs.get(tmp_frag),
            					jp_factor, jp_summant));
            		tmp_val =
            			((Double)prob_set.get(org_str).get(ind - 1)).doubleValue() + log_tmp;
            		prob_set.get(org_str).add(tmp_val);
            	}
            }

            String high_spec_name = getHighestProb(prob_set);
            if (high_spec_name.equals(spec_name))
            	my_keeper.increment(kmer_size);
            	
            my_keeper.increment_div(kmer_size);
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

            Iterator<String> itr = spec_set.iterator();

            while(itr.hasNext())
            {
                String next = itr.next();

                ArrayList<Double> list = new ArrayList<Double>();
                list.add(new Double(0));
                ret_set.put(next, list);
            }

            return ret_set;
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
            FilenameFilter filter = new FilenameFilter()
            {
            	public boolean accept(File dir, String name)
            	{
            		return name.startsWith("TRAIN");
            	}
            };
            File files[] = input_dir_h.listFiles(filter);
            HashMap<String, HashMap<String, Double>> ret_table =
                new HashMap<String, HashMap<String, Double>>();
            Double zero = new Double(0.0);

            Vector<String> kmer_strings = getStringVector();

            for (int i = 0; i < files.length; i++)
            {
                if ((files[i].getName()).endsWith(Integer.toString(kmer_size)))
                {
                    String[] file_name_toks = files[i].getName().split(",");
                    System.out.println("name: "+files[i].getName());
                    String spec_name = file_name_toks[2].replace(",", " ");
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
                            assert cur_string.equals(tokens[0]);
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
