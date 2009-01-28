package genend.msvm;
import java.util.*;
import java.io.*;

public class MSVMDataBuilder
{
	private String input_dir = null, output_dir = null;
	private String output_filename = null;
	public MSVMDataBuilder(String input_dir, String output_dir, String output_filename)
	{
		this.input_dir = input_dir;
		this.output_dir = output_dir;
		this.output_filename = output_filename;
	}
	
	public void buildMSVMFile()
	{
		File dir_h = new File(input_dir);
		File[] file_list = dir_h.listFiles();
		String line, write_ln;
		int count = 1, file_count = 1;
		BufferedWriter bw = null, spec_index_bw = null;
		try
		{
			bw = new BufferedWriter(new FileWriter(output_dir+"/"+output_filename));
			spec_index_bw = 
				new BufferedWriter(new FileWriter(output_dir+"/"+output_filename+"-spec_index"));
		}
		catch (Exception e) {}
		
		for (int i = 0; i < file_list.length; i++)
		{
			String filename = file_list[i].toString();
			if (filename == "." || filename == "..")
				continue;
			write_ln = Integer.toString(file_count)+" ";
			try
			{
				BufferedReader br = new BufferedReader(new FileReader(filename));
				while ((line = br.readLine()) != null)
				{
					String[] tokens = line.split("\t");
					write_ln = write_ln + Integer.toString(count++)+":"+tokens[1]+" ";
				}
				br.close();
				bw.write(write_ln+"\n");
				spec_index_bw.write(file_list[i].getName() + "\t" + Integer.toString(file_count++)+"\n");
			}
			catch (Exception e) { e.printStackTrace(); }
			System.out.println("Done processing \'"+filename+"\'");
			count = 0;
		}
		
		try
		{
			bw.close();
			spec_index_bw.close();
		}
		catch (Exception e) { e.printStackTrace(); }
	}
	
	public static void main(String[] args) 
	{
		if (args.length < 3) { System.out.println("ERROR: need 3 arguments"); return; }
		MSVMDataBuilder test = new MSVMDataBuilder(args[0], args[1], args[2]);
		test.buildMSVMFile();
	}

}
