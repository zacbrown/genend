package genend.util.container;
public class ResultKeeper
{
    private double[] nums, divs;
    private int kmer_low;
    public ResultKeeper(int kmer_low, int kmer_high)
    {
        this.kmer_low = kmer_low;
        nums = new double[kmer_high - kmer_low + 1];
        divs = new double[kmer_high - kmer_low + 1];
    }

    public synchronized void increment(int kmer)
    {
        nums[kmer - kmer_low]++;
    }

    public synchronized void increment_div(int kmer)
    {
    	divs[kmer - kmer_low]++;
    }
    
    public double get(int kmer)
    {
        return nums[kmer - kmer_low];
    }
    
    public double get_div(int kmer)
    {
    	return divs[kmer - kmer_low];
    }

    public String toString()
    {
        String ret = "{";
        for (int i = 0; i < nums.length - 1; i++)
        {
            ret = ret + (i+kmer_low) + ":" + nums[i] + ", ";
        }
        ret = ret + (nums.length - 1 + kmer_low) + ":" + nums[nums.length - 1] + "}";

        return ret;
    }
}
