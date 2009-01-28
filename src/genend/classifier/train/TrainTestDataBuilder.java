package genend.classifier.train;

import genend.util.MersenneTwisterFast;
import genend.util.SeqUtils;
import genend.util.container.SeqObj;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;

public class TrainTestDataBuilder
{
    private int kmer_low, kmer_high, num_threads;
    private String output_dir, input_dir;
    private int powers_of_four[] = new int[20];
    private int piece_size;
    private boolean rep;
    private double perc_train, divisor;


    public TrainTestDataBuilder(int piece_size, int kmer_low, int kmer_high, double perc_train, 
    		String input_dir, String output_dir, boolean rep, int num_threads)
    {
        this.piece_size = piece_size;
        this.kmer_low = kmer_low;
        this.kmer_high = kmer_high;
        this.input_dir = input_dir;
        this.output_dir = output_dir;
        this.perc_train = perc_train;
        this.rep = rep;
        this.num_threads = num_threads;

        for(int i = 0; i < powers_of_four.length; i++)
            powers_of_four[i] = (int)Math.pow(4.0, (double)i);
    }

    public double getDivisor() { return divisor; }
    
    public void execute()
    {
        ThreadPoolExecutor tpe = null;
        ConcurrentLinkedQueue<SeqObj> frag_queue = new ConcurrentLinkedQueue<SeqObj>();
        ConcurrentLinkedQueue<SeqObj> org_queue = new ConcurrentLinkedQueue<SeqObj>();
          
        File input_dir_h = new File(input_dir);
        File[] input_files = input_dir_h.listFiles();

        System.out.println("STARTING SAMPLING");
        /* Sample the genomes */
        tpe = new ThreadPoolExecutor(num_threads, num_threads,
        		25000L, TimeUnit.MILLISECONDS,
        		new LinkedBlockingQueue<Runnable>());

        for (int i = 0; i < input_files.length; i++)
        {
        	Runnable r = new FragGetter(input_files[i].toString(), perc_train, piece_size,
        			rep, kmer_low, kmer_high);
        	tpe.execute(r);
        }
        tpe.shutdown();
        try { tpe.awaitTermination(10000000L, TimeUnit.SECONDS); }
        catch (Exception e) { e.printStackTrace(); }
        System.out.println("DONE SAMPLING\n");
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

    private class FragGetter extends Thread
    {
    	private String input_file = null;
    	private double perc_test;
    	private int piece_size, kmer_low, kmer_high;
    	private boolean rep;
    	private ConcurrentLinkedQueue<SeqObj> frag_queue = null, org_queue = null;
    	private MersenneTwisterFast rg = null;
    	
    	FragGetter(String input_file, double perc_train, int piece_size, boolean rep, int kmer_low, int kmer_high)
    	{
    		this.input_file = input_file;
    		this.kmer_low = kmer_low;
    		this.kmer_high = kmer_high;
    		this.rep = rep;
    		perc_test = 1.0 - perc_train;
            this.piece_size = piece_size;
    		rg = new MersenneTwisterFast();
    	}
    	
    	public void run()
    	{
    		String[] cur_org_info = parse_genome_file(input_file);
    		String org_name = cur_org_info[0];
    		StringBuffer org_seq = new StringBuffer(cur_org_info[1]);
    		int org_len = org_seq.length();
    		int num_frags = org_len / piece_size;
    		int num_sample = (int)(perc_test * num_frags);
    		Vector<SeqObj> tmp_frags = new Vector<SeqObj>();

    		if (org_len < 500000) return;
    		
    		System.out.println("FRAGS: "+org_name);
    		
    		for (int i = 0; i < num_sample; i++)
    		{
    			String frag = getFragment(org_seq, org_len, piece_size, rep);
    		    tmp_frags.add(new SeqObj(org_name, frag, org_len, SeqObj.Type.FRAG, i));
    		}
    	
            SeqObj tmp_full_obj = new SeqObj(org_name, org_seq.toString(), org_len, SeqObj.Type.FULL, -1);
            if (tmp_full_obj.getSeq() != null)
            {
        		new DistribProc(tmp_full_obj, kmer_low, kmer_high, "TRAIN").run();
                int size = tmp_frags.size();

                for (int i = 0; i < size; i++)
                	new FragWriter(tmp_frags.get(i)).run();
            }
    	}
    	
    	String getFragment(StringBuffer seq_str, int seq_len, int frag_size, boolean rep)
        {
    		int piece_pos;
            while (true)
            {
            	piece_pos = rg.nextInt(seq_len - frag_size + 1) + 1; // swap piece_size with frag_size ?

            	if (SeqUtils.matchBases(seq_str.substring(piece_pos,
            			piece_pos + frag_size)) == true)
            		break;
            }
            
            int start_pos = piece_pos;
            int end_pos = piece_pos + frag_size;
            String ret_str = new String(seq_str.substring(start_pos,end_pos));
            
            if (rep == false)
            {
            	StringBuilder str_builder = new StringBuilder();
            	for (int i = start_pos; i <= end_pos; i++)
            		str_builder.append("X");
            	seq_str.replace(start_pos, end_pos+1, str_builder.toString());
            }
            return ret_str;
        }
    	
    }
    
    private class FragWriter extends Thread
    {
    	private SeqObj seq_obj = null;
    	
    	FragWriter(SeqObj seq_obj)
    	{
    		this.seq_obj = seq_obj;
    	}
    	
    	public void run()
    	{
    		String spec_name = seq_obj.getOrganism();
    		int spec_index = seq_obj.getIndex();
    		
    		BufferedWriter bw = null;
            try
        	{
        		String filename = "FRAG," + spec_index + "," + 
        			spec_name.replace(" ", "_");
        		filename = filename.replace("/", "#");
        		bw = new BufferedWriter(new FileWriter(output_dir+"/"+filename));
        		bw.write(seq_obj.getSeq());
        		bw.close();
        	}
        	catch (Exception e) { e.printStackTrace(); }
    	}
    }
    
    private class DistribProc extends Thread
    {
    	Hashtable<String, Double> distrib_tbl = null;
    	String spec_name = null, seq_str = null, prefix = null;
    	SeqObj.Type seq_type;
    	int kmer_low, kmer_high, seq_len, spec_index;
    	
    	DistribProc(SeqObj seq_obj, int kmer_low, int kmer_high, String prefix)
    	{
    		this.distrib_tbl = new Hashtable<String, Double>();
    		seq_str = seq_obj.getSeq();
    		spec_name = seq_obj.getOrganism();
    		spec_index = seq_obj.getIndex();
    		seq_type = seq_obj.getType();
    		if (seq_type == SeqObj.Type.FULL) 
    			seq_len = seq_obj.getOrgLen();
    		else
    			seq_len = seq_str.length();
    		this.prefix = prefix;
    		this.kmer_low = kmer_low;
    		this.kmer_high = kmer_high;
    	}
    	
    	public void run()
    	{
    		if (seq_str == null) return;
    		else seq_str = seq_str.toUpperCase();
    		System.out.println("DISTRIB PROC: "+prefix+"-"+spec_name);
    		for (int kmer_size = kmer_low; kmer_size <= kmer_high; kmer_size++)
    		{
    			/* building distributions */
    			
                for (int ind = 0; ind < seq_len - kmer_size + 1; ind++)
                {
                    String frag = seq_str.substring(ind, ind + kmer_size);
                    if (SeqUtils.matchBases(frag))
                    {
                        if (distrib_tbl.containsKey(frag))
                        {
                            double new_val = distrib_tbl.get(frag).doubleValue() + 1.0;
                            distrib_tbl.put(frag, new Double(new_val));
                        }
                        else distrib_tbl.put(frag, new Double(1.0));
                    }
                }

                Iterator itr = distrib_tbl.entrySet().iterator();

                while(itr.hasNext())
                {
                    Map.Entry frag_entry = (Map.Entry)itr.next();
                    String frag = (String)frag_entry.getKey();
                    double old_val = ((Double)frag_entry.getValue()).doubleValue();
                    double new_val = old_val / ((double)seq_len);
                    distrib_tbl.put(frag, new_val);
                }
                
                /* build relative distributions */
                int count = 0;
                double tmp_val = 0.0, total_val = 0.0;
                ArrayList<String> tmp_kmer_list = new ArrayList<String>();
                ArrayList<Double> tmp_val_list = new ArrayList<Double>();
              
                for (int j = 0; j < powers_of_four[kmer_size]; j++)
                {
                    String kmer_str = SeqUtils.numToStr(j, kmer_size);
                    tmp_kmer_list.add(kmer_str);

                    if (distrib_tbl.containsKey(kmer_str))
                        tmp_val = distrib_tbl.get(kmer_str).doubleValue();
                    else
                    {
                        tmp_val = 0.0;
                        distrib_tbl.put(kmer_str, new Double(0.0));
                    }

                    total_val += tmp_val;
                    tmp_val_list.add(new Double(tmp_val));

                    if (count == 3)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            if (total_val == 0.0)
                                distrib_tbl.put(tmp_kmer_list.get(k), new Double(0.0));
                            else
                                distrib_tbl.put(tmp_kmer_list.get(k),
                                                 new Double(tmp_val_list.get(k) / total_val));
                        }

                        count = 0;
                        total_val = 0.0;
                        tmp_kmer_list.clear();
                        tmp_val_list.clear();
                    }
                    else
                        count++;
                }
                
                /* write distribs to file */
                BufferedWriter bw = null;
                try
            	{
            		String filename = prefix + "," + spec_index + "," + 
            			spec_name.replace(" ", "_") + "," + kmer_size;
            		filename = filename.replace("/", "_");
            		bw = new BufferedWriter(new FileWriter(output_dir+"/"+filename));
            		for (int j = 0; j < powers_of_four[kmer_size]; j++)
                	{
            			String frag = SeqUtils.numToStr(j, kmer_size);
            			Double val = distrib_tbl.get(frag);
            			bw.write(frag + "\t" + val.toString()+"\n");
                	}
            		bw.close();
            	}
            	catch (Exception e) { e.printStackTrace(); }
    		}
    	}
    }
}
