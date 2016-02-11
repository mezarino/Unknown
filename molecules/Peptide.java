package molecules;
/*
*
*
*/

import java.io.*;
import java.util.TreeMap;

import molecules.MassInfo;

public class Peptide
{
    private String sequence;
    private double monomass;
    private double weight;
    
    public Peptide()
    {}
    
    public Peptide(String sequence)
    {
        this.sequence = sequence;
        monomass = MassInfo.getMass(sequence, true);
        weight = MassInfo.getWeight(sequence, true);
    }
    
    public static boolean isSequenceValid(String sequence) 
    {
        for (int i = 0; i < sequence.length(); i++)
        {
            if ( MassInfo.aminoAcids.contains(sequence.charAt(i)) )
                continue;
            
            return false;
        }
        
        return true;
    }
    
    public static TreeMap<Character, Integer> getAAComposition(String sequence)
    {
        TreeMap<Character, Integer> result = MassInfo.getCompositionTemplate();
        
        if (sequence.length() > 0)
        {
            result.put('n', 1); // nterminal
            result.put('c', 1); // cterminal
        }
        
        for (int i = 0; i < sequence.length(); i++)
        {
            Character aa = sequence.charAt(i);
            result.put(aa, result.get(aa) + 1);
        }
        
        return result;
    }
    
    public static void main(String[] args)
    {
        String s = "AGTFAH";
        TreeMap<Character, Integer> comp = getAAComposition(s);
        for (int i = 0; i < s.length(); i++)
            System.out.println(comp.get(s.charAt(i)));
    }
}
