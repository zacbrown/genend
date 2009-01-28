package genend.util.container;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;

public class SeqObj 
{
	public enum Type { FRAG, FULL };
	
	private String organism = null;
	private SymbolList seq = null;
	private int org_len, count;
	private final Type type;
	
	public SeqObj(String organism, String seq_str, int org_len, Type type, int count)
	{
		this.organism = organism;
		this.org_len = org_len;
		this.count = count;
		this.type = type;
		Alphabet dna = DNATools.getDNA();
		SymbolTokenization dnaToke = null;
		try 
		{ 
			dnaToke = dna.getTokenization("token"); 
			seq = new SimpleSymbolList(dnaToke, seq_str);
		}
		catch(Exception e) { e.printStackTrace(); }
	}
	
	public int getIndex() { return count; }
	public int getOrgLen() { return org_len; }
	public String getOrganism() { return organism; }
	public Type getType() { return type; }
	public String getSeq() 
	{
		if (seq == null) return null;
		Alphabet dna = DNATools.getDNA();
		String seq_str = null;
		try
		{
			SymbolTokenization dnaToke = dna.getTokenization("token");
			seq_str = dnaToke.tokenizeSymbolList(seq);
		}
		catch (Exception e) { e.printStackTrace(); }
		return seq_str.toUpperCase(); 
	} 
}
