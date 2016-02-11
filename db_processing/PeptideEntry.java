package db_processing;
/*
*
*/

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

import molecules.MassInfo;
import molecules.Mass;
import molecules.Peptide;
import utils.ArrayFunctions;


public class PeptideEntry
{
    private short length;
    private double monomass = Double.NaN;
    private List<Integer> proteinIndices;
    private List<Short> proteinOffsets;
    
    public PeptideEntry(int length, double mass, int protIdx, int protOffset)
    {
        this.length = (short)length;
        this.monomass = mass;
        proteinIndices = new ArrayList<Integer>(1);
        proteinOffsets = new ArrayList<Short>(1);
        proteinIndices.add(protIdx);
        proteinOffsets.add((short)protOffset);
    }
    
    public PeptideEntry(int length, double mass, int nprots)
    {
        this.length = (short)length;
        this.monomass = mass;
        proteinIndices = new ArrayList<Integer>(nprots);
        proteinOffsets = new ArrayList<Short>(nprots);
    }
    
    public PeptideEntry(String sequence, int protIdx, int protOffset)
    {
        length = (short) sequence.length();
        monomass = MassInfo.getMass(sequence, true);
        proteinIndices = new ArrayList<Integer>(1);
        proteinOffsets = new ArrayList<Short>(1);
        proteinIndices.add(protIdx);
        proteinOffsets.add((short)protOffset);
    }
    
    public void addOccurrence(int protIdx, int protOffset)
    {
        proteinIndices.add(protIdx);
        proteinOffsets.add((short)protOffset);
    }
    
    public static PeptideEntry read(RandomAccessFile reader)
    {
        try
        {
            short length = reader.readShort();
            double mass = reader.readDouble();
            int nproteins = reader.readInt();
            PeptideEntry entry = new PeptideEntry(length, mass, nproteins);
            
            for (int i = 0; i < nproteins; i++)
            {
                entry.proteinIndices.add(reader.readInt());
                entry.proteinOffsets.add(reader.readShort());
            }
        
            return entry;
        }
        catch (IOException e)
        {
            return null;
        }
    }
    
    public void write(RandomAccessFile writer) throws IOException
    {
        writer.writeShort(length);
        writer.writeDouble(monomass);
        writer.writeInt(proteinIndices.size());
        for (int i = 0; i < proteinIndices.size(); i++)
        {
            writer.writeInt(proteinIndices.get(i));
            writer.writeShort(proteinOffsets.get(i));
        }
    }
    
    public int proteinIndexAt(int i)
    {
        return proteinIndices.get(i);
    }
    
    public int proteinOffsetAt(int i)
    {
        return proteinOffsets.get(i);
    }
    
    public List<Integer> getProteinIndices()
    {
        return proteinIndices; //proteinIndices.toArray(new Integer[0]);  
    }
    
    public List<Short> getProteinOffsets()
    {
        return proteinOffsets; //proteinOffsets.toArray(new Short[0]);
    }
    
    public short getLength()
    {
        return length;
    }
    
    public double getMonoisotopicMass()
    {
        return monomass;
    }
    
    public int getNumberOfOccurrences()
    {
        return proteinOffsets.size();
    }
    
    public HashSet<Integer> getIndicesOfContainingProteins()
    {
        HashSet<Integer> result = new HashSet<Integer>(proteinIndices);
        return result;
    }
    
    public int getNumberOfContainingProteins()
    {
        HashSet<Integer> result = new HashSet<Integer>(proteinIndices);
        return result.size();
    }
    
    public String getSequence(ProteinsIndexation proteins)
    {
        //return proteins.getProtein(proteinIndices.get(0)).getSequence().substring(proteinOffsets.get(0), proteinOffsets.get(0)+length);
        return proteins.getSubsequenceFromProtein(proteinIndices.get(0), proteinOffsets.get(0), length);
    }
    
    public char getResidueBefore(int index, ProteinsIndexation proteins)
    {
        int pos = proteinOffsets.get(index) - 1;
        return proteins.getResidueAtProteinPos(proteinIndices.get(index), pos);
    }

    public char getResidueAfter(int index, ProteinsIndexation proteins)
    {
        int pos = proteinOffsets.get(index) + length;
        return proteins.getResidueAtProteinPos(proteinIndices.get(index), pos);
    }
    
    public boolean isAtProteinNterm(ProteinsIndexation proteins)
    {
        for (int i = 0; i < proteinOffsets.size(); i++)
        {
            if (proteinOffsets.get(i) == 0)
                return true;
            else if (proteinOffsets.get(i) == 1 && proteins.getResidueAtProteinPos(proteinIndices.get(i), 0) == 'M')
                return true;
        }
        
        return false;
    }

    public boolean isAtProteinCterm(ProteinsIndexation proteins)
    {
        try
        {
            for (int i = 0; i < proteinOffsets.size(); i++)
            {
                int pos = proteinOffsets.get(i) + length;
                int protLen = proteins.getLength(proteinIndices.get(i));
                if (pos >= protLen - 1)
                    return true;
            }
            
            return false;
        }
        catch (IOException e)
        {
            System.err.println(e.getMessage());
            return false;
        }
    }
    
    public void applyFixedModifications(String pepseq, Map<Character, Double> deltaMasses)
    {
        double deltaMass = 0.0;
        
        if (deltaMasses.containsKey('n'))
            deltaMass += deltaMasses.get('n');
        if (deltaMasses.containsKey('c'))
            deltaMass += deltaMasses.get('c');
        
        for (int i = 0; i < pepseq.length(); i++)
        {
            char aa = pepseq.charAt(i);
            if (deltaMasses.get(aa) != null)
                deltaMass += deltaMasses.get(aa);
        }
        
        monomass += deltaMass;
    }
    
    public static Mass[] applyVariableModifications(String pepseq, double pepmass, Map<Character, Double> deltaMasses, long pepidx, int maxmods)
    {
        List<Mass> result = new ArrayList<Mass>();
        result.add(new Mass(pepmass, pepidx, 0)); // The non-modified version of the peptide; 0 is the number of modifications.;
        
        TreeMap<Character, Integer> seqComp = Peptide.getAAComposition(pepseq);
        Character[] sites = deltaMasses.values().toArray(new Character[deltaMasses.size()]);
        
        for (int l = 1; l <= maxmods; l++)
        {
            String[] modCombinations = ArrayFunctions.getCombinations(l, sites); // All modifications combinations of size l (note is a proper combination, not a permutation)
            for (int n = 0; n < modCombinations.length;  n++)
            {
                TreeMap<Character, Integer> modCombComp = ArrayFunctions.getCharComposition(modCombinations[n]);
                Character[] sites2 = modCombComp.values().toArray(new Character[modCombComp.size()]);
                double modifiedMass = pepmass;
                boolean flag = true;
                
                for (int j = 0; j < sites2.length; j++)
                {
                    if ( modCombComp.get(sites2[j]) > seqComp.get(sites2[j]) )
                    {
                        flag = false;
                        break;
                    }
                    
                    modifiedMass += deltaMasses.get(sites2[j]) * modCombComp.get(sites2[j]);
                }
                
                if (flag)
                    result.add(new Mass(modifiedMass, pepidx, modCombinations[n].length()));
            }
        }
        
        return result.toArray(new Mass[result.size()]);
    }
    
    public void extend(PeptideEntry[] peps)
    {
        for (int i = 0; i < peps.length; i++)
        {
            if (peps[i] != null)
            {
                this.proteinIndices.addAll(peps[i].getProteinIndices());
                this.proteinOffsets.addAll(peps[i].getProteinOffsets());
            }
        }
    }
    
    public static PeptideEntry merge(PeptideEntry[] peps) // The user must make sure the entries have the same sequence
    {
        PeptideEntry result = null;
        
        for (int i = 0; i < peps.length; i++)
        {
            if (peps[i] == null)
                continue;
            else if (result == null)
                result = new PeptideEntry(peps[i].getLength(), peps[i].getMonoisotopicMass(), peps[i].getNumberOfOccurrences());
            else
               break;
        }
        
        result.extend(peps);
        
        return result;
    }
    
    public static PeptideEntry copyOf(PeptideEntry pep)
    {
        if (pep == null)
            return null;
        
        int n = pep.getNumberOfOccurrences();
        PeptideEntry entry = new PeptideEntry(pep.getLength(), pep.getMonoisotopicMass(), n);
        
        for (int i = 0; i < n; i++)
            entry.addOccurrence(pep.proteinIndexAt(i), pep.proteinOffsetAt(i));
        
        return entry;
        
    }
}



