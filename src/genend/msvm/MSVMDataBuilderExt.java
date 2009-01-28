package genend.msvm;
import genend.util.SeqUtils;

import java.util.*;
import java.io.*;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.concurrent.*;



public class MSVMDataBuilderExt
{
	private String input_dir = null, output_dir = null;
	private String output_filename = null;
    private int powers_of_four[] = new int[20];
    private int kmer_size, num_threads, num_samples;
    private int[] piece_sizes;
	
	public MSVMDataBuilderExt(String input_dir, String output_dir, String output_filename, 
			int kmer_size, int[] piece_sizes, int num_samples, int num_threads)
	{
		this.input_dir = input_dir;
		this.output_dir = output_dir;
		this.output_filename = output_filename;
		this.kmer_size = kmer_size;
		this.piece_sizes = piece_sizes;
		this.num_threads = num_threads;
		this.num_samples = num_samples;
		for(int i = 0; i < powers_of_four.length; i++)
            powers_of_four[i] = (int)Math.pow(4.0, (double)i);
	}
	
	public void execute()
	{
		ThreadPoolExecutor tpe = null;
		ConcurrentLinkedQueue<String> ret_queue = new ConcurrentLinkedQueue<String>();
        BufferedWriter bw = null;
        String[] num_to_str = new String[powers_of_four[kmer_size]];
        for (int i = 0; i < powers_of_four[kmer_size]; i++)
        	num_to_str[i] = SeqUtils.numToStr(i, kmer_size);
        try { bw = new BufferedWriter(new FileWriter(output_dir+"/"+output_filename+"-spec")); }
        catch (Exception e) { e.printStackTrace(); }
        
        File input_dir_h = new File(input_dir);
        File[] input_files = input_dir_h.listFiles();
        
        Runnable svm_writer = new SVMWriter(output_dir+"/"+output_filename, ret_queue);
        tpe = new ThreadPoolExecutor(num_threads, num_threads,
        		25000L, TimeUnit.MILLISECONDS,
        		new LinkedBlockingQueue<Runnable>());
        tpe.execute(svm_writer);
        
        for (int i = 0; i < input_files.length; i++)
        {
        	try { bw.write(input_files[i].getName() + "\t" + i + "\n"); }
        	catch (Exception e) { e.printStackTrace(); }
        	Runnable r = new SamplerThread(input_files[i].toString(), piece_sizes,
        			kmer_size, num_samples, i, ret_queue, tpe, num_to_str);

        	tpe.execute(r); 	
        }

        tpe.shutdown();
        while (tpe.getCompletedTaskCount() - input_files.length != 0){}
        ((SVMWriter)svm_writer).setDone();
        try { tpe.awaitTermination(10000000L, TimeUnit.SECONDS); }
        catch (Exception e) { e.printStackTrace(); }
        
        try { bw.close(); }
    	catch (Exception e) { e.printStackTrace(); }
	}
	
	private String[] parse_genome_file(String path)
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
        catch (Exception e) { }
        finally { try { if (br != null) br.close(); } catch (Exception e) {} }

        return ret;
    }

    private String parse_header(String header)
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
	
    private class SVMWriter extends Thread
    {
    	private String output_filename = null;
    	private ConcurrentLinkedQueue<String> queue = null;
    	private boolean done = false;
    	
    	SVMWriter(String output_filename, ConcurrentLinkedQueue<String> queue)
    	{
    		this.output_filename = output_filename;
    		this.queue = queue;
    	}
    	
    	public void setDone()
    	{
    		done = true;
    	}
    	
    	public void run()
    	{
    		BufferedWriter bw = null;
    		try 
    		{
    			bw = new BufferedWriter(new FileWriter(output_filename));
    			
    			while (!done)
    			{
    				if (queue.size() > 250)
    				{
    					while (queue.size() > 0) bw.write((String)queue.poll()+"\n");
    				}
    			}
    			while (queue.size() > 0) bw.write((String)queue.poll()+"\n");
    			bw.close();
    		}
    		catch (Exception e) { e.printStackTrace(); }
    	}
    	
    }
    
	private class SamplerThread extends Thread
    {
        private int kmer_size, spec_id, total_frags = 0, sample_size;
        private int[] piece_sizes;
        private ConcurrentLinkedQueue<String> queue;
        private String input_file, spec_name;
        private HashMap<String, Double> distribs;
        private Random rg;
        private ThreadPoolExecutor tpe;
        private String[] num_to_str = null;

        SamplerThread(String input_file, int[] piece_sizes, int kmer_size, int sample_size, 
        		int spec_id, ConcurrentLinkedQueue<String> ret_queue, ThreadPoolExecutor tpe, String[] num_to_str)
        {
            this.kmer_size = kmer_size;
            this.piece_sizes = piece_sizes;
            this.input_file = input_file;
            this.spec_id = spec_id;
            this.sample_size = sample_size;
            this.distribs = new HashMap<String, Double>();
            this.queue = ret_queue;
            this.num_to_str = num_to_str;
            this.tpe = tpe;
            rg = new Random();
        }

        public void run()
        {
            String[] cur_org_info = parse_genome_file(input_file);
            NumberFormat formatter = new DecimalFormat("0.###############");
            spec_name = cur_org_info[0];
            String org_seq = cur_org_info[1];
            int org_len = org_seq.length();
            
            System.out.println("PROCESSING: "+spec_name);
            for (int piece_size : piece_sizes)
            {
            	for (int i = 0; i < sample_size; i++)
            	{
            		String str_kmer = getFragment(org_seq, org_len, piece_size, piece_size);
                    int count = 0, str_kmer_len;

            		while (!SeqUtils.matchBases(str_kmer) && count < 1000)
                    {
                        str_kmer = getFragment(org_seq, org_len, kmer_size, piece_size);
                        count++;
                    }
                    str_kmer_len = str_kmer.length();

                    for (int ind = 1; ind < str_kmer_len - kmer_size + 1; ind++)
                    {
                        String tmp_frag = str_kmer.substring(ind, ind + kmer_size);
                        double value = 0;
                        if (distribs.containsKey(tmp_frag))
                        	value = (distribs.get(tmp_frag)).doubleValue();
                        value += 1;
                        distribs.put(tmp_frag, Double.valueOf(value));
                    }
                    
                    StringBuilder vect_str_bldr = new StringBuilder();
                    vect_str_bldr.append(spec_id + " ");
                    String line = "";
                    for (int k = 0; k < powers_of_four[kmer_size]; k++)
                    {
                    	String test_kmer = num_to_str[k];
                    	if (distribs.containsKey(test_kmer) == false)
                    		line = (k+1) + ":" + 0 + " ";
                    	else
                    	{
                    		double value = (distribs.get(test_kmer)).doubleValue();
                    		value /= org_len;
                    		line = (k+1) + ":" + formatter.format(value) + " ";
                    	}
                    	vect_str_bldr.append(line);
                    }

                    queue.add(vect_str_bldr.toString());
            	}
            	
            	total_frags = 0;
            }
        }

        String getFragment(String seq_str, int seq_len, int frag_size, int piece_size)
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
	
	public static void main(String[] args) 
	{
		int num_samples_train = 100;
        int num_samples_test = 25;
		int[] piece_sizes_train = {200};		
        int[] piece_sizes_test = {1000};
		
		if (args.length < 5) { System.out.println("ERROR: need 3 arguments"); return; }
		MSVMDataBuilderExt test = new MSVMDataBuilderExt(args[0], args[1], args[2],
				Integer.valueOf(args[3]), piece_sizes_train, num_samples_train,
				Integer.valueOf(args[4]));
		test.execute();
        test = new MSVMDataBuilderExt(args[0], args[1], "test"+args[2], Integer.valueOf(args[3]),
            piece_sizes_test, num_samples_test, Integer.valueOf(args[4]));
        test.execute();
	}

}
