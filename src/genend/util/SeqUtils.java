package genend.util;
import java.util.regex.*;

public class SeqUtils
{
    private static final char bases[] = {'A', 'C', 'G', 'T'};
    private static Pattern base_pat = Pattern.compile("[ACGT]*");

    public static String numToStr(int num, int frag_size)
    {
        String ret_str = "";
        int cur_num = num;

        for (int i = 0; i < frag_size; i++)
        {
            ret_str = bases[cur_num % 4] + ret_str;
            cur_num /= 4;
        }

        return ret_str;
    }

    public static double jp(double num, double factor, double summant)
    {
        return num * factor + summant;
    }

    public static boolean matchBases(String str)
    {
        return base_pat.matcher(str).matches();
    }

    public static String parseOrgName(String header)
    {
        String[] tokens = header.split("\\s");
        return tokens[1] + "_" + tokens[2];
    }
}
