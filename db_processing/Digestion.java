/*
*
*/

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.RandomAccessFile;
import java.io.PrintWriter;
import java.util.List;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import molecules.Enzyme;

public class Digestion 
{
    public Digestion()
    {}
    
    
    public static void digest(String fastaFile, Enzyme enzyme, int missed_cleavages, boolean semispecific, int min_length, int max_length, String outFile)
    {
        // Iterate through fasta file
        String ID = "", DE = "", seq = "";
        long offset = 0;
        
        try
		{   
            //RandomAccessFile reader = new RandomAccessFile(fastaFile, "r");
            BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
            RandomAccessFile writer1 = new RandomAccessFile(outFile+"_protein", "rw");
            RandomAccessFile writer2 = new RandomAccessFile(outFile+"_protein-index", "rw");
            String st = reader.readLine();
            
            //PrintWriter writer1 = new PrintWriter(outFile+"_protein", "UTF-8");
            //PrintWriter writer2 = new PrintWriter(outFile+"_protein-index", "UTF-8");
            //PrintWriter writer3 = new PrintWriter(outFile+"_peptide", "UTF-8");
            //PrintWriter writer4 = new PrintWriter(outFile+"_mass", "UTF-8");
            
            while (st != null)
	        {   
                st = st.trim();
                
                if (st.startsWith(">"))
                {
                    List<String> digestion = getDigestion(seq, enzyme, missed_cleavages, semispecific);
                    for (int i = 0; i < digestion.size(); i++)
                    {
                        String contaminant = "", reverse = "";
    
                        if (digestion.get(i).length() >= min_length && digestion.get(i).length() <= max_length)
                        {
                            String pepseq = digestion.get(i);
                            //weight = calculateMass(pepseq, modifications1, modifications2, modifications3);
                            if ( ID.startsWith("CON_") )
                                contaminant = "+";
                            else if ( ID.startsWith("REV_") )
                                reverse = "+";
                            
                            //writer1.println(ID + "\t" + DE + "\t" + pepseq + "\t" + weight + "\t" + reverse + "\t" + contaminant);
                            //writer2.println(ID + "\t" 
                            writer1.writeChars(pepseq + "\t" + ID + "\t" + DE + "\n");
                            String line = "" + offset + "\t" + pepseq.length() + "\t" + ID + "\n";
                            writer2.writeChars(line);
                            offset = writer1.getFilePointer();
                        }
                    }
                    
                    ID = st.split(" ")[0].substring(1);
                    DE = st.substring(1);
                    seq = "";
                }
                else
                {
                    seq = seq.concat(st);
                }
                
                st = reader.readLine();
            }
            
            // Process last sequence in fasta file
            List<String> digestion = getDigestion(seq, enzyme, missed_cleavages, semispecific);
            for (int i = 0; i < digestion.size(); i++)
            {
                if (digestion.get(i).length() >= min_length && digestion.get(i).length() <= max_length)
                {
                    //writer1.println(digestion.get(i) + "\t" + ID);
                    String pepseq = digestion.get(i);
                    writer1.writeChars(pepseq + "\t" + ID + "\t" + DE + "\n");
                    String line = "" + offset + "\t" + pepseq.length() + "\t" + ID + "\n";
                    writer2.writeChars(line);
                    offset = writer1.getFilePointer();
                }
            }
            
            reader.close();
            writer1.close();
            writer2.close();
        }
        catch (Exception e)
		{ System.err.println(e.toString()); }
    }
    
    public static void readFasta(String fastaFile, Enzyme enzyme, int missed_cleavages, int min_length, String outFile)
    {
        // Iterate through fasta file
        String ID = "", seq = "";
        
        
        try
		{
            //RandomAccessFile reader = new RandomAccessFile(fastaFile, "r");
            BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
            //RandomAccessFile writer = new RandomAccessFile(outFile, "rw");
            PrintWriter writer = new PrintWriter(outFile, "UTF-8");
            String st = reader.readLine();
            
            while (st != null)
	        {   
                st = st.trim();
                
                if (st.startsWith(">"))
                {
                    writer.println(">" + ID);
                    writer.println(seq);
                    ID = st.split(" ")[0].substring(1);
                    seq = "";
                }
                else
                {
                    seq = seq.concat(st);
                }
                
                st = reader.readLine();
            }
            
            writer.println(">" + ID);
            writer.println(seq);
            
            reader.close();
            writer.close();
        }
        catch (Exception e)
		{ System.err.println(e.toString()); }
    }
    
    public static List<String> getDigestion(String sequence,  Enzyme enzyme, int missed_cleavages, boolean semispecific)
    {
        String cutSide = enzyme.getCutSide();
        Pattern pattern = Pattern.compile(enzyme.getCleavageSitePattern());
        Matcher matcher = pattern.matcher(sequence);
        List<String> digestion = new ArrayList<String>(); // try initializing with size: sequence.length() * missed_cleavages / 15
        int lastend = 0, n0 = 0;

        if (cutSide.equals("C-term"))
        {
            while (matcher.find()) 
            {
                String current = sequence.substring(lastend, matcher.end());
                digestion.add(current);
                lastend = matcher.end();
                n0++;
            }

            if (lastend != sequence.length())
            {
                String last = sequence.substring(lastend, sequence.length());
                digestion.add(last);
                n0++;
            }
        }
        
        else if (cutSide.equals("N-term"))
        {
            while (matcher.find()) 
            {
                String current = sequence.substring(lastend, matcher.start());
                digestion.add(current);
                lastend = matcher.start();
                n0++;
            }

            if (lastend != sequence.length())
            {
                String last = sequence.substring(lastend, sequence.length());
                digestion.add(last);
                n0++;
            }
        }
        
        
        // Add miscleavaged
        for (int i = 0; i < n0; i++)
        {
            int j = 1;
            while (j <= missed_cleavages && i+j < n0)
            {
                String s = concatenateItems(digestion.subList(i, i+j+1));
                digestion.add(s);
                j++;
            }
        }
        
        // Add semispecific
        if (semispecific)
        {
            int n1 = digestion.size();
            
            for (int i = 0; i < n1; i++)
            {
                String s1 = digestion.get(i);
                
                for (int j = 1; j < s1.length(); j++)
                {
                    digestion.add(s1.substring(0, j));
                    digestion.add(s1.substring(j));
                }
            }
        }
        
        return digestion;
        
    }
    
    
    public static void printMatches(String text, String regex)
    {
        Pattern pattern = Pattern.compile(regex);
        Matcher matcher = pattern.matcher(text);
        // Check all occurrences
        while (matcher.find()) {
            System.out.print("Start index: " + matcher.start());
            System.out.print(" End index: " + matcher.end());
            System.out.println(" Found: " + matcher.group());
        }
    }
    
    public static String concatenateItems(String[] array)
    {
         String result = "";
        
        for (int i = 0; i < array.length; i++)
        {
            result = result.concat(array[i]);
        }
        
        return result;
    }
    
    public static String concatenateItems(List<String> list)
    {
        String result = "";
        
        for (int i = 0; i < list.size(); i++)
        {
            result = result.concat(list.get(i));
        }
        
        return result;
    }
    
    public static void main(String[] args) 
    {
        Enzyme trypsin = new Enzyme(1);
        String protein = "MMRPGFKQRLIKKTTGSSSSSSSKKKDKEKEKEKSSTTSSTSKKPASASSSSHGTTHSSASSTGSKSTTEKGKQSGSVPSQGKHHSSSTSKTKTATTPSSSSSSSRSSSVSRSGSSSTKKTSSRKPRGQEQRP";
        List<String> digestion = getDigestion(protein, trypsin, 2, false);
        for (int i = 0; i < digestion.size(); i++)
        {
            System.out.println(digestion.get(i));
        }
        
        // try with a fasta file
	Enzyme enzyme = new Enzyme(Integer.parseInt(args[1]));
        digest(args[0], enzyme, Integer.parseInt(args[2]), Boolean.valueOf(args[3]), Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[6]);
        //readFasta(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), args[4]);
    }
}
