import genend.classifier.OverlapBayesianClassifier;

public class OverlapMain
{

    public static void main(String[] args)
    {
        int[] piece_sizes = {10000};
        int kmer_min = 3, kmer_max = 8;
        String output_dir = "/home/zbrown/genomics/dists-all/";
        String input_dir = "/home/zbrown/genomics/overlap-genomes/";
        String conf_yaml = "/home/zbrown/genomics/genend/config.yml";
        int num_threads = 6;


        OverlapBayesianClassifier test = new OverlapBayesianClassifier(piece_sizes, kmer_min, kmer_max,
                input_dir, output_dir, false, num_threads);

        test.execute();
    }

}
