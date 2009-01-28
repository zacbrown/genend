package genend.util.container;
import java.util.Hashtable;

public class RelDistResultObj
{
    private int kmer_size;
    private Hashtable<String, Double> rel_distribs;
    private String spec_name;

    public RelDistResultObj(String spec_name, int kmer_size,
                            Hashtable<String, Double> rel_distribs)
    {
        this.spec_name = spec_name;
        this.kmer_size = kmer_size;
        this.rel_distribs = rel_distribs;
    }

    public String getSpecName() { return spec_name; }
    public int getKmerSize() { return kmer_size; }
    public Hashtable<String, Double> getRelDistribs() { return rel_distribs; }
}
