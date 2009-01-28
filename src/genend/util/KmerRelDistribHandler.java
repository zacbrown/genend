package genend.util;

import genend.util.container.RelDistribObj;

import java.io.*;
import java.util.Hashtable;


public class KmerRelDistribHandler
{
    private int kmer_low, kmer_high;
    private String output_dir, filename_prefix;

    public KmerRelDistribHandler(int kmer_low, int kmer_high, String output_dir)
    {
        this.kmer_low = kmer_low;
        this.kmer_high = kmer_high;
        this.output_dir = output_dir;
        filename_prefix = output_dir + "/relative-distrib-";
    }

    public int add(int kmer_num, String spec_name,
                   Hashtable<String, Double> distribs)
    {
        ObjectOutputStream out_h = null;
        ObjectInputStream in_h = null;
        String filename = filename_prefix + String.valueOf(kmer_num);

        Hashtable<String, Hashtable<String, Double>> cur_distribs = null;
        File file_h = new File(filename);

        if (file_h.exists())
        {
            try
            {
                in_h = new ObjectInputStream(new FileInputStream(filename));
                cur_distribs = ((RelDistribObj)in_h.readObject()).getDistribs();
                in_h.close();
                cur_distribs.put(spec_name, distribs);
            }
            catch (Exception e) { e.printStackTrace(); }
        }
        else
        {
            cur_distribs = new Hashtable<String, Hashtable<String, Double>>();
            cur_distribs.put(spec_name, distribs);
        }

        try
        {
            out_h = new ObjectOutputStream(new FileOutputStream(filename));
            out_h.writeObject(new RelDistribObj(kmer_num, cur_distribs));
            out_h.close();
        }
        catch (Exception e) { e.printStackTrace(); }

        cur_distribs = null;
        in_h = null;
        out_h = null;
        filename = null;

        return kmer_num;
    }

    public void delete(int kmer_num)
    {
        String filename = filename_prefix + String.valueOf(kmer_num);
        File file_h = new File(filename);
        assert file_h.delete();
    }

    public RelDistribObj get(int kmer_num)
    {
        String filename = filename_prefix + String.valueOf(kmer_num);
        ObjectInputStream file_h = null;
        RelDistribObj ret_obj = null;

        try
        {
            file_h = new ObjectInputStream(new FileInputStream(filename));
            ret_obj = (RelDistribObj)file_h.readObject();
            file_h.close();
        }
        catch (Exception e) { e.printStackTrace(); }

        file_h = null;
        filename = null;

        return ret_obj;
    }
}
