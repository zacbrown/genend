package genend.util.container;
public class ClassifyResultObj
{
    private int kmer_size, piece_size, name_val,
        kingdom_val, phylum_val, class_val, order_val, family_val,
        genus_val;

    public ClassifyResultObj(int kmer_size, int piece_size, int name_val,
                             int kingdom_val, int phylum_val, int class_val,
                             int order_val, int family_val, int genus_val)
    {
        this.kmer_size = kmer_size;
        this.piece_size = piece_size;
        this.name_val = name_val;
        this.kingdom_val = kingdom_val;
        this.phylum_val = phylum_val;
        this.class_val = class_val;
        this.order_val = order_val;
        this.family_val = family_val;
        this.genus_val = genus_val;
    }

    public int getKmerSize() { return kmer_size; }
    public int getNameMatchVal() { return name_val; }
    public int getKingdomMatchVal() { return kingdom_val; }
    public int getPhylumMatchVal() { return phylum_val; }
    public int getClassMatchVal() { return class_val; }
    public int getOrderMatchVal() { return order_val; }
    public int getFamilyMatchVal() { return family_val; }
    public int getGenusMatchVal() { return genus_val; }
    public int getPieceSize() { return piece_size; }
    public String toString()
    {
        String ret = "{kmer_size: " + kmer_size;
        ret += ", piece_size: " + piece_size;
        //        ret += ", match_val: " + match_val + "}";
        return ret;
    }
}
