package genend.msvm;
import genend.util.SeqUtils;

import java.util.*;
import java.io.*;
import java.text.*;
import java.util.concurrent.*;



public class MSVMDataBuilderFullGen
{
	private String input_dir = null, output_dir = null;
	private String output_filename = null;
    private int powers_of_four[] = new int[20];
    private int kmer_size, num_threads;
	
	public MSVMDataBuilderFullGen(String input_dir, String output_dir, String output_filename, 
			int kmer_size, int num_threads)
	{
		this.input_dir = input_dir;
		this.output_dir = output_dir;
		this.output_filename = output_filename;
		this.kmer_size = kmer_size;
		this.num_threads = num_threads;
		for(int i = 0; i < powers_of_four.length; i++)
            powers_of_four[i] = (int)Math.pow(4.0, (double)i);
	}
	
	public void execute()
	{
		ThreadPoolExecutor tpe = null;
        Vector<String> ret_vector = new Vector<String>();
        String[] num_to_str = new String[powers_of_four[kmer_size]];
        for (int i = 0; i < powers_of_four[kmer_size]; i++)
        	num_to_str[i] = SeqUtils.numToStr(i, kmer_size); 
        BufferedWriter bw = null;
        try { bw = new BufferedWriter(new FileWriter(output_dir+"/"+output_filename+"-spec")); }
        catch (Exception e) { e.printStackTrace(); }
        
        File input_dir_h = new File(input_dir);
        File[] input_files = input_dir_h.listFiles();
        System.out.println(input_files);
        tpe = new ThreadPoolExecutor(num_threads, num_threads,
        		25000L, TimeUnit.MILLISECONDS,
        		new LinkedBlockingQueue<Runnable>());

        for (int i = 0; i < input_files.length; i++)
        {
        	try { bw.write(input_files[i].getName() + "\t" + i + "\n"); }
        	catch (Exception e) { e.printStackTrace(); }
        	Runnable r = new SamplerThread(input_files[i].toString(),
        			kmer_size, i, ret_vector, tpe, num_to_str);

        	tpe.execute(r); 	
        }

        tpe.shutdown();
        try { tpe.awaitTermination(10000000L, TimeUnit.SECONDS); }
        catch (Exception e) { e.printStackTrace(); }
        
        try { bw.close(); }
    	catch (Exception e) { e.printStackTrace(); }
        writeFile(ret_vector);
	}

	private void writeFile(Vector<String> process_vect)
	{
		BufferedWriter bw = null;
		System.out.println("Writing output to: "+output_dir+"/"+output_filename);
		try
		{
			bw = new BufferedWriter(new FileWriter(output_dir+"/"+output_filename));
			for (int i = 0; i < process_vect.size(); i++)
				bw.write(process_vect.get(i)+"\n");
			bw.close();
		}
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
	
	private class SamplerThread extends Thread
    {
        private int kmer_size, spec_id;
        private Vector<String> ret_vector;
        private String input_file, spec_name;
        private HashMap<String, Double> distribs;
        private String[] num_to_str = null;
        private ThreadPoolExecutor tpe;

        SamplerThread(String input_file, int kmer_size, int spec_id, 
        		Vector<String> ret_vector, ThreadPoolExecutor tpe, String[] num_to_str)
        {
            this.kmer_size = kmer_size;
            this.input_file = input_file;
            this.spec_id = spec_id;
            this.distribs = new HashMap<String, Double>();
            this.ret_vector = ret_vector;
            this.num_to_str = num_to_str;
            this.tpe = tpe;
        }

        public void run()
        {
            String[] cur_org_info = parse_genome_file(input_file);
            NumberFormat formatter = new DecimalFormat("0.###############");
            spec_name = cur_org_info[0];
            String org_seq = cur_org_info[1];
            int org_len = org_seq.length();
            
            System.out.println("PROCESSING: "+spec_name);

            for (int ind = 1; ind < org_len - kmer_size + 1; ind++)
            {
            	String tmp_frag = org_seq.substring(ind, ind + kmer_size);
            	double value = 0;
            	if (distribs.containsKey(tmp_frag))
            	{
            		value = (distribs.get(tmp_frag)).doubleValue();
            	}
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

            ret_vector.add(vect_str_bldr.toString());      	
        }
    }
	
	public static void main(String[] args) 
	{
		int num_threads = 12;		
		
		if (args.length < 4) { System.out.println("ERROR: need 3 arguments"); return; }
		MSVMDataBuilderFullGen test = new MSVMDataBuilderFullGen(args[0], args[1], args[2],
				Integer.valueOf(args[3]), num_threads);
		test.execute();
	}

}
