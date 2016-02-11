/*
*
*
*/

import java.io.*;
import java.util.Hashtable;

public class AminoAcid
{
    private static AminoAcid[] aminoAcids;
	private static Hashtable<Character, Double> aaMonoMasses;
	private static Hashtable<Character, Double> aaOccurences;
	private static Hashtable<Character, Double> aaWeights;
	private static Hashtable<String, Character> fromCodon;
	private static String singleLetterAas;
    private String formula;
    private String name;
	private String abbreviation;
	private char letter;
	private double occurence;
    private String[] codons;
    
    public AminoAcid()
    {}
    
	public static AminoAcid[] getAAs()
    {
		if (aminoAcids == null)
        {
			aminoAcids = initAminoAcids();
	    }
		
        return aminoAcids;
	}
    
    public static AminoAcid[] initAminoAcids()
    {
        AminoAcid alanine = new AminoAcid("C3H5NO", "Alanine", "Ala", 'A', 7.4, new String[]{"GCT", "GCC", "GCA", "GCG"});
		AminoAcid arginine = new AminoAcid("C6H12N4O", "Arginine", "Arg", 'R', 4.2, new String[]{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"});
		AminoAcid asparagine = new AminoAcid("C4H6N2O2", "Asparagine", "Asn", 'N', 4.4, new String[]{"AAT", "AAC"});
		AminoAcid asparticAcid = new AminoAcid("C4H5NO3", "Aspartic Acid", "Asp", 'D', 5.9, new String[]{"GAT", "GAC"});
		AminoAcid cysteine = new AminoAcid("C3H5NOS", "Cysteine", "Cys", 'C', 3.3, new String[]{"TGT", "TGC"});
		AminoAcid glutamicAcid = new AminoAcid("C5H7NO3", "Glutamic Acid", "Glu", 'E', 5.8, new String[]{"GAA", "GAG"});
		AminoAcid glutamine = new AminoAcid("C5H8N2O2", "Glutamine", "Gln", 'Q', 3.7, new String[]{"CAA", "CAG"});
		AminoAcid glycine = new AminoAcid("C2H3NO", "Glycine", "Gly", 'G', 7.4, new String[]{"GGT", "GGC", "GGA", "GGG"});
		AminoAcid histidine = new AminoAcid("C6H7N3O", "Histidine", "His", 'H', 2.9, new String[]{"CAT", "CAC"});
		AminoAcid isoleucine = new AminoAcid("C6H11NO", "Isoleucine", "Ile", 'I', 3.8, new String[]{"ATT", "ATC", "ATA"});
		AminoAcid leucine = new AminoAcid("C6H11NO", "Leucine", "Leu", 'L', 7.6, new String[]{"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"});
		AminoAcid lysine = new AminoAcid("C6H12N2O", "Lysine", "Lys", 'K', 7.2, new String[]{"AAA", "AAG"});
		AminoAcid methionine = new AminoAcid("C5H9NOS", "Methionine", "Met", 'M', 1.8, new String[]{"ATG"});
		AminoAcid phenylalanine = new AminoAcid("C9H9NO", "Phenylalanine", "Phe", 'F', 4.0, new String[]{"TTT", "TTC"});
		AminoAcid proline = new AminoAcid("C5H7NO", "Proline", "Pro", 'P', 5.0, new String[]{"CCT", "CCC", "CCA", "CCG"});
		AminoAcid serine = new AminoAcid("C3H5NO2", "Serine", "Ser", 'S', 8.1, new String[]{"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"});
		AminoAcid threonine = new AminoAcid("C4H7NO2", "Threonine", "Thr", 'T', 6.2, new String[]{"ACT", "ACC", "ACA", "ACG"});
		AminoAcid tryptophan = new AminoAcid("C11H10N2O", "Tryptophan", "Trp", 'W', 1.3, new String[]{"TGG"});
		AminoAcid tyrosine = new AminoAcid("C9H9NO2", "Tyrosine", "Tyr", 'Y', 3.3, new String[]{"TAT", "TAC"});
		AminoAcid valine = new AminoAcid("C5H9NO", "Valine", "Val", 'V', 6.8, new String[]{"GTT", "GTC", "GTA", "GTG"});
        
		AminoAcid[] aas = new AminoAcid[]{alanine, arginine, asparagine, asparticAcid,
                                          cysteine, glutamine, glutamicAcid, glycine,
			               	              histidine, isoleucine, leucine, lysine, methionine,
                                          phenylalanine, proline, serine, threonine,
                                          tryptophan, tyrosine, valine};
        
		return aas;
    }
    
    
    private AminoAcid(String formula, String name, String abbreviation, char letter, double occurence, String[] codons)
    {
            this.formula = formula;
            this.name = name;
			this.abbreviation = abbreviation;
			this.letter = letter;
			this.occurence = occurence/100.0;
			this.codons = codons;
    }
    
    public String getAbbreviation()
    {
		return abbreviation;
	}

	public char getLetter()
    {
		return letter;
	}
    
    public String getFormula()
    {
        return formula;
    }
    
    public String getName()
    {
        return name;
    }    

    private static void main(String[] args)
    {
        AminoAcid myAA = new AminoAcid();
        AminoAcid[] aas = myAA.getAAs();
    }
    
}
