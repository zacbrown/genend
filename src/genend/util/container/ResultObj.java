package genend.util.container;
public class ResultObj
{
    private int match_val, kmer_size, piece_size;

    public ResultObj(int kmer_size, int piece_size, int match_val)
    {
        this.kmer_size = kmer_size;
        this.match_val = match_val;
        this.piece_size = piece_size;
    }

    public int getKmerSize() { return kmer_size; }
    public int getMatchVal() { return match_val; }
    public int getPieceSize() { return piece_size; }
    public String toString()
    {
        String ret = "{kmer_size: " + kmer_size;
        ret += ", piece_size: " + piece_size;
        ret += ", match_val: " + match_val + "}";
        return ret;
    }
}
