import genend.classifier.ClassifyProb2;
import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;
import org.biojavax.bio.seq.RichSequence;


public class Main
{

    public static void main(String[] args)
    {
        int[] piece_sizes = {36};
        int kmer_min = 3, kmer_max = 5;
        String output_dir = "/home/zbrown/genomics/dists-archaea/";
        String input_dir = "/home/zbrown/genomics/genomes-archaea/";
        String conf_yaml = "/home/zbrown/genomics/genend/config.yml";
        int num_threads = 6;

        /*        StatProb2 test = new StatProb2(piece_sizes, kmer_min, kmer_max,
                  input_dir, output_dir, false, num_threads);*/
        
        ClassifyProb2 test = new ClassifyProb2(piece_sizes, kmer_min, kmer_max,
                                               input_dir, output_dir, false, num_threads,
                                               conf_yaml);

        test.execute();
    }

}
