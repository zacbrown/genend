package genend.util;

import genend.util.container.KmerObj;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;



public class KmerDistribProcessor
{
    private String genome_dir, output_dir;
    private int kmer_low, kmer_high, num_threads;
    private int powers_of_four[] = new int[20];
    private FileOutputStream index_file_h;
    private KmerObjHandler kmer_obj_h;
    private KmerRelDistribHandler kmer_rel_distrib_h;
    private String index_file_name;

    public KmerDistribProcessor(String genome_dir, String output_dir,
                                int kmer_low, int kmer_high, int num_threads)
    {
        this.genome_dir = genome_dir;
        this.output_dir = output_dir;
        this.kmer_low = kmer_low;
        this.kmer_high = kmer_high;
        this.num_threads = num_threads;
        index_file_name = this.output_dir + ".index";
        try
        {
            File test_index_file = new File(index_file_name);
            if (!test_index_file.exists())
                index_file_h = new FileOutputStream(index_file_name);
            kmer_obj_h = new KmerObjHandler(output_dir);
            kmer_rel_distrib_h =
                new KmerRelDistribHandler(kmer_low, kmer_high, output_dir);
        }
        catch (Exception e) { e.printStackTrace(); }

        for(int i = 0; i < powers_of_four.length; i++)
            powers_of_four[i] = (int)Math.pow(4.0, (double)i);
    }

    protected void finalize() throws Throwable
    {
        index_file_h.close();
        super.finalize();
    }

    public void processFiles()
    {
        HashSet<String> spec_set = new HashSet<String>();

        File gen_dir_h = new File(genome_dir);
        File files[] = gen_dir_h.listFiles();
        String spec_name = "";
        String seq_str = null, organism = null;
        StringBuilder seq_builder = null;
        Hashtable<String, Double> rel_distribs = null;
        Hashtable<String, Double> frag_distribs = null;
        BufferedReader br = null;

        for (int i = 0; i < files.length; i++)
        {
            File filename = files[i];
            System.out.println("Processing file \'" + filename + "\'...");

            System.gc();

            seq_builder = new StringBuilder();
            int seq_len = 0;
            try
            {
                br = new BufferedReader(new FileReader(filename));
                organism = SeqUtils.parseOrgName(br.readLine());

                while ((seq_str = br.readLine()) != null)
                    seq_builder.append(seq_str.toUpperCase());
                br.close();

                seq_str = seq_builder.toString();
                seq_len = seq_str.length();
            }
            catch (Exception e) { e.printStackTrace(); }

            for (int kmer_size = kmer_low; kmer_size <= kmer_high; kmer_size++)
            {
                System.out.println("Generating "+kmer_size+"-mer distributions...");

                frag_distribs = new Hashtable<String, Double>();

                for (int ind = 0; ind < seq_len - kmer_size; ind++)
                {
                    String frag = seq_str.substring(ind, ind + kmer_size);

                    if (SeqUtils.matchBases(frag))
                    {
                        if (frag_distribs.containsKey(frag))
                        {
                            double new_val = frag_distribs.get(frag).doubleValue() + 1.0;
                            frag_distribs.put(frag, new Double(new_val));
                        }
                        else frag_distribs.put(frag, new Double(1.0));
                    }
                    else continue;
                }

                Iterator itr = frag_distribs.entrySet().iterator();

                while(itr.hasNext())
                {
                    Map.Entry frag_entry = (Map.Entry)itr.next();
                    String frag = (String)frag_entry.getKey();
                    double old_val = ((Double)frag_entry.getValue()).doubleValue();
                    double new_val = old_val / ((double)seq_len);
                    frag_distribs.put(frag, new_val);
                }

                spec_name = kmer_obj_h.add(organism, seq_str, kmer_size, frag_distribs);
            }

            // Relative distribution calculations.
            KmerObj my_obj = kmer_obj_h.get(spec_name);

            for (int kmer_size = kmer_low; kmer_size <= kmer_high; kmer_size++)
            {
                System.out.println("Generating "+kmer_size+"-mer relative distributions...");

                rel_distribs = new Hashtable<String, Double>();

                int count = 0;
                double tmp_val = 0.0, total_val = 0.0;
                ArrayList<String> tmp_kmer_list = new ArrayList<String>();
                ArrayList<Double> tmp_val_list = new ArrayList<Double>();
                Hashtable<String, Double> kmer_dict =
                    my_obj.getDistribs().get(Integer.valueOf(kmer_size));

                for (int j = 0; j < powers_of_four[kmer_size]; j++)
                {
                    String kmer_str = SeqUtils.numToStr(j, kmer_size);
                    tmp_kmer_list.add(kmer_str);

                    if (kmer_dict.containsKey(kmer_str))
                        tmp_val = kmer_dict.get(kmer_str).doubleValue();
                    else
                    {
                        tmp_val = 0.0;
                        kmer_dict.put(kmer_str, new Double(0.0));
                    }

                    total_val += tmp_val;
                    tmp_val_list.add(new Double(tmp_val));

                    if (count == 3)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            if (total_val == 0.0)
                                rel_distribs.put(tmp_kmer_list.get(k), new Double(0.0));
                            else
                                rel_distribs.put(tmp_kmer_list.get(k),
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

                kmer_rel_distrib_h.add(kmer_size, spec_name, rel_distribs);
            }
            spec_set.add(spec_name);
            System.out.print("Done!\n");

            br = null;
            my_obj = null;
            rel_distribs = null;
            frag_distribs = null;
            seq_builder = null;
        }

        ObjectOutputStream out_h = null;
        try
        {
            out_h = new ObjectOutputStream(new FileOutputStream(index_file_name));
            out_h.writeObject(spec_set);
            out_h.close();
        }
        catch (Exception e) { e.printStackTrace(); }
    }

    public HashSet<String> getIndex()
    {
        FileInputStream fis_h = null;
        ObjectInputStream obj_h = null;
        HashSet<String> ret_obj = null;
        System.out.println(index_file_name);
        try
        {
            obj_h = new ObjectInputStream(new FileInputStream(index_file_name));
            ret_obj = (HashSet<String>)obj_h.readObject();
            obj_h.close();
        }
        catch (Exception e) { e.printStackTrace(); }

        return ret_obj;
    }

    public KmerObjHandler get_kmer_obj_h()
    {
        return kmer_obj_h;
    }
}
