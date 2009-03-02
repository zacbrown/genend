package genend.util.container;
public class ResultObj
{
    private int match_val, kmer_size, piece_size;
    private double score;
    private String cur_spec, high_spec;

    public ResultObj(int kmer_size, int piece_size, int match_val,
            String cur_spec, String high_spec, double score)
    {
        this.kmer_size = kmer_size;
        this.match_val = match_val;
        this.piece_size = piece_size;
        this.cur_spec = cur_spec;
        this.high_spec = high_spec;
        this.score = score;

    }

    public int getKmerSize() { return kmer_size; }
    public int getMatchVal() { return match_val; }
    public int getPieceSize() { return piece_size; }
    public String getCurSpec() { return cur_spec; }
    public String getHighSpec() { return high_spec; }
    public double getScore() { return score; }
    public String toString()
    {
        String ret = "{kmer_size: " + kmer_size;
        ret += ", piece_size: " + piece_size;
        ret += ", match_val: " + match_val;
        ret += ", score: " + score;
        ret += ", cur_spec: " + cur_spec;
        ret += ", high_spec: " + high_spec + "}";
        return ret;
    }
}
