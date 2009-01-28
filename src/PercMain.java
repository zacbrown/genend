import genend.classifier.train.TrainTestDataBuilder;

public class PercMain
{
    public static void main(String[] args)
    {
        int piece_size = 500;
        int kmer_min = 3, kmer_max = 3;
        double perc_train = 0.90;
        String output_dir = ".\\";
        String input_dir = ".\\genomes-all\\";
        int num_threads = 2;

        /*        StatProb2 test = new StatProb2(piece_sizes, kmer_min, kmer_max,
                  input_dir, output_dir, false, num_threads);*/

        TrainTestDataBuilder test = new TrainTestDataBuilder(piece_size, kmer_min, kmer_max, perc_train,
                                               input_dir, output_dir, false, num_threads);

        test.execute();
    }

}
