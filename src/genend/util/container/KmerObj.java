package genend.util.container;
import java.io.Serializable;
import java.util.Hashtable;

public class KmerObj implements Serializable
{
    private String org_name, seq;
    private Hashtable<Integer, Hashtable<String, Double>> distribs;

    public KmerObj(String org_name, String seq, Hashtable<Integer, Hashtable<String, Double>> distribs)
    {
        this.org_name = org_name;
        this.seq = seq;
        this.distribs = distribs;
    }

    public String getSequence()
    {
        return seq;
    }

    public String getOrgName()
    {
        return org_name;
    }

    public Hashtable<Integer, Hashtable<String, Double>> getDistribs()
    {
        return distribs;
    }
}
