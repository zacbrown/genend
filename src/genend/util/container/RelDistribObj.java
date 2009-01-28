package genend.util.container;
import java.io.Serializable;
import java.util.Hashtable;

public class RelDistribObj implements Serializable
{
    private int kmer_num;
    private Hashtable<String, Hashtable<String, Double>> distribs;

    public RelDistribObj(int kmer_num, Hashtable<String, Hashtable<String, Double>> distribs)
    {
        this.kmer_num = kmer_num;
        this.distribs = distribs;
    }

    public int getKmer()
    {
        return kmer_num;
    }

    public Hashtable<String, Hashtable<String, Double>> getDistribs()
    {
        return distribs;
    }
}
