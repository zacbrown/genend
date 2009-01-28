package genend.util;

import genend.util.container.KmerObj;

import java.io.*;
import java.util.Hashtable;


public class KmerObjHandler
{
    private String output_dir;

    public KmerObjHandler(String output_dir)
    {
        this.output_dir = output_dir;
    }

    public String add(String org_name, String seq, int kmer_num,
                      Hashtable<String, Double> kmer_distribs)
    {
        String spec_name = org_name;
        String filename = output_dir + "/" + spec_name;

        ObjectOutputStream out_h = null;
        ObjectInputStream in_h = null;
        Hashtable<Integer, Hashtable<String, Double>> cur_distribs = null;
        File file_h = new File(filename);

        if (file_h.exists())
        {
            try
            {
                in_h = new ObjectInputStream(new FileInputStream(filename));
                cur_distribs = ((KmerObj)in_h.readObject()).getDistribs();
                in_h.close();
                cur_distribs.put(Integer.valueOf(kmer_num), kmer_distribs);
            }
            catch (Exception ex) { ex.printStackTrace(); }
        }
        else
        {
            cur_distribs = new Hashtable<Integer, Hashtable<String, Double>>();
            cur_distribs.put(kmer_num, kmer_distribs);
        }

        try
        {
            out_h = new ObjectOutputStream(new FileOutputStream(filename));
            out_h.writeObject(new KmerObj(spec_name, seq, cur_distribs));
            out_h.close();
        }
        catch (Exception ex) { ex.printStackTrace(); }

        out_h = null;
        in_h = null;
        cur_distribs = null;
        filename = null;

        return spec_name;
    }

    public void delete(String species_name)
    {
        String filename = output_dir + species_name.replace(" ", "_");
        File file_h = new File(filename);
        assert file_h.delete();
    }

    public KmerObj get(String species_name)
    {
        String filename = output_dir + species_name.replace(" ", "_");
        ObjectInputStream file_h = null;
        KmerObj ret_obj = null;

        try
        {
            file_h = new ObjectInputStream(new FileInputStream(filename));
            ret_obj = (KmerObj)file_h.readObject();
            file_h.close();
        }
        catch (Exception e) { e.printStackTrace(); }

        file_h = null;
        filename = null;

        return ret_obj;
    }

    private String parse_fasta(String in_str)
    {
        String[] tokens = in_str.split("\\s");
        return tokens[0] + "_" + tokens[1];
    }
}
