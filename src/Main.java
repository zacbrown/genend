import genend.classifier.ClassifyProb2;
import genend.classifier.BayesianClassifier;

public class Main
{

    public static void main(String[] args)
    {
        int[] piece_sizes = {36,100,200,400,800};
        int kmer_min = 3, kmer_max = 8;
        String output_dir = "/home/zbrown/genomics/dists-all/";
        String input_dir = "/home/zbrown/genomics/dists-all/";
        String conf_yaml = "/home/zbrown/genomics/genend/config.yml";
        int num_threads = 6;


        BayesianClassifier test = new BayesianClassifier(piece_sizes, kmer_min, kmer_max,
                input_dir, output_dir, false, num_threads);
        
        /*ClassifyProb2 test = new ClassifyProb2(piece_sizes, kmer_min, kmer_max,
                                               input_dir, output_dir, false, num_threads,
                                               conf_yaml);*/

        test.execute();
    }

}
