import genend.classifier.ClassifyProb2;
import genend.classifier.StatProb2;

public class Main
{

    public static void main(String[] args)
    {
        int[] piece_sizes = {36};
        int kmer_min = 3, kmer_max = 5;
        String output_dir = "/home/zbrown/genomes/dists-archaea/";
        String input_dir = "/home/zbrown/genomes/genomes-archaea/";
        String conf_yaml = "/home/zbrown/genend/config.yml";
        int num_threads = 6;

        StatProb2 test = new StatProb2(piece_sizes, kmer_min, kmer_max,
                input_dir, output_dir, false, num_threads);
        
        /*ClassifyProb2 test = new ClassifyProb2(piece_sizes, kmer_min, kmer_max,
                                               input_dir, output_dir, false, num_threads,
                                               conf_yaml);*/

        test.execute();
    }

}
