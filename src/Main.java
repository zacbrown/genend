import genend.classifier.ClassifyProb2;
import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;
import org.biojavax.bio.seq.RichSequence;
import genend.util.GenbankSeqFetcher;


public class Main
{

    public static void main(String[] args)
    {
        int[] piece_sizes = {36};
        int kmer_min = 3, kmer_max = 3;
        String output_dir = "../dists-all/";
        String input_dir = "../genomes-all/";
        String db_name = "taxonomy.db";
        int num_threads = 6;

        /*        StatProb2 test = new StatProb2(piece_sizes, kmer_min, kmer_max,
                  input_dir, output_dir, false, num_threads);*/
        /*
        ClassifyProb2 test = new ClassifyProb2(piece_sizes, kmer_min, kmer_max,
                                               input_dir, output_dir, false, num_threads,
                                               db_name);

        test.execute();*/

        String genomePath = "C:\\Users\\Zac\\Code\\Genend\\data\\test-genomes";
        String configPath = "C:\\Users\\Zac\\Code\\Genend\\config.yml";
        GenbankSeqFetcher myFetcher = 
                new GenbankSeqFetcher(genomePath, 1, configPath, "default");

        myFetcher.execute();
    }

}
