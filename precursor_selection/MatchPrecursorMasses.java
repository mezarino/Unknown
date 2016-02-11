package precursor_selection;
/*
*
*/

import java.io.*;
import java.util.List;
import java.util.ArrayList;

import molecules.MassInfo;

public class MatchPrecursorMasses
{
    public MatchPrecursorMasses()
    {}
    
    public static double[][] mgfParser(String mgfFilePath)
    {
        try
		{
            RandomAccessFile reader = new RandomAccessFile(mgfFilePath, "r");
            
            int i = 0;
            int scan = -1, charge = -1;
            double precmz = -1, monomass = -1;
            
            String hd = reader.readLine().trim();
            int n = Integer.parseInt(hd.split("#")[1]); // the number of MS2 scans
            double[][] result = new double[n][2];
            
            String st = reader.readLine();
            
            while (st != null)
	        {   
                st = st.trim();
                
                if (st.startsWith("TITLE="))
                {
                    st = st.split("=")[1];
                    scan = Integer.parseInt(st.split(":")[1]);
                }
                
                else if (st.startsWith("CHARGE="))
                {
                    st = st.split("=")[1];
                    charge = Integer.parseInt(st.substring(0, st.length()-1));
                }
                
                else if (st.startsWith("PEPMASS="))
                {
                    precmz = Double.parseDouble(st.split("=")[1]);
                }
                
                else if (st.startsWith("END IONS")) // calculate monomass and push prec
                {
                    monomass = precmz * charge - MassInfo.proton * charge;
                    result[i][0] = monomass;
                    result[i][1] = scan;
                    i++;
                }
                
                else
                {}
                st = reader.readLine();
            }
            
            reader.close();
            return result;
        }
        
        catch (Exception e)
		{ e.toString(); }
        
        return null;
    }
    
    
    public static double[][] indexPeptideFasta(String fastaFilePath)
    {
        try
		{
            RandomAccessFile reader = new RandomAccessFile(fastaFilePath, "r");
            
            String hd = reader.readLine().trim();
            int n = Integer.parseInt(hd.split("#")[1]); // the number of sequences
            double[][] result = new double[2][n];  // contain masses and offsets of sequences
            
            int i = 0;
            long offset = 0;
            String st = reader.readLine();
            
            while (st != null && i < n)
	        {
                st = st.trim();
                
                // All lines in file are full and ordered by mass. First section is pep seq
                String[] content = st.split("\t");
                result[0][i] = Double.parseDouble(content[1]);
                result[1][i] = offset;
                offset = reader.getFilePointer();
                
                // Binary search
                //Array.binarySearch(masses, precmass);
                st = reader.readLine();
                i++;
            }
            
            reader.close();
            return result;
        }
        
        catch (Exception e)
		{ e.toString(); }
        
        return null;
    }
    
    
    public static void main (String[] args)
    {
        // Reading parameters
        String filePath1 = args[0];
        String filePath2 = args[1];
        
        double[][] pairs = mgfParser(filePath1);
        for (int i = 0; i < pairs.length; i++)
        {
            System.out.println("Precursor mass: " + pairs[i][0] + " MS2 scan: " + pairs[i][1]);
        }
        
        double[][] massAndIdxs = indexPeptideFasta(filePath2);
        for (int i = 0; i < massAndIdxs[0].length; i++)
        {
            System.out.println("Database mass: " + massAndIdxs[0][i] + " Binary offset: " + massAndIdxs[1][i]);
        }
    }
}
