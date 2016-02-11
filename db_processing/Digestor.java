package db_processing;

/*
*
*/

import java.util.List;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import molecules.Enzyme;
import molecules.MassInfo;
import molecules.Peptide;
import molecules.Modification;
import utils.ArrayFunctions;

public class Digestor 
{
    private List<String> digestion = new ArrayList<String>();
    private List<Short> positions = new ArrayList<Short>();
    private List<Short> lengths = new ArrayList<Short>();
    private List<Double> masses = new ArrayList<Double>();
    
    public Digestor()
    {}
    
    public void digest(String sequence,  Enzyme enzyme, int missed_cleavages, boolean semispecific, int min_length, int max_length)
    {
        ArrayFunctions myFunctions = new ArrayFunctions();
        
        String restrictedAA = enzyme.getRestrictedAA();
        Pattern pattern = Pattern.compile(enzyme.getRegex());
        Matcher matcher = pattern.matcher(sequence);
         
        int laststart = 0, lastend = 0, n0 = 0;
        String last = null;
        
        if (matcher.find())
        {
            last = matcher.group();
            laststart = matcher.start();
            lastend = matcher.end();
        }
        
        while (matcher.find()) 
        {
            String current = matcher.group();
            if ( !(restrictedAA.equals("")) && current.startsWith(restrictedAA) )
            {
                last = last.concat(current);
                continue;
            }
            
            digestion.add(last);
            positions.add((short) laststart);
            last = current;
            laststart = matcher.start();
            lastend = matcher.end();
            n0++;
        }
        
        if (lastend != sequence.length())
        {
            String tail = sequence.substring(lastend, sequence.length());
            
            if ( !restrictedAA.equals("") && tail.startsWith(restrictedAA) )
            {
                digestion.add(last + tail);
                positions.add((short) laststart);
                n0++;
            }
            else if (last != null)
            {
                digestion.add(last);
                positions.add((short) laststart);
                digestion.add(tail);
                positions.add((short) lastend);
                n0 += 2;
            }
            else
            {
                digestion.add(tail);
                positions.add((short) lastend);
                n0++;
            }
        }
        else if (last != null)
        {
            digestion.add(last);
            positions.add((short) laststart);
            n0++;
        }
        
        // Add miscleavaged
        for (int i = 0; i < n0; i++)
        {
            int j = 1;
            while (j <= missed_cleavages && i+j < n0)
            {
                String s = myFunctions.concatenateItems(digestion.subList(i, i+j+1));
                digestion.add(s);
                positions.add(positions.get(i));
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
                    positions.add(positions.get(i));
                    digestion.add(s1.substring(j));
                    positions.add((short)(positions.get(i) + j));
                }
            }
        }

        // Cleaved methionine from those peptides that start at offset 0 and add them to digestion
        int n2 = digestion.size();
        for (int i = 0; i < n2; i++)
        {
            String s2 = digestion.get(i);
            if (positions.get(i) == 0 && s2.startsWith("M"))
            {
                digestion.add(s2.substring(1));
                positions.add((short) 1);
            }
        }            

        // Filter sequences by length
        int i = 0;
        while (i < digestion.size())
        {
            if (digestion.get(i).length() < min_length || digestion.get(i).length() > max_length || !Peptide.isSequenceValid(digestion.get(i)))
            {
                digestion.remove(i);
                positions.remove(i);
                continue;
            }
            
            lengths.add((short) digestion.get(i).length());
            masses.add(MassInfo.getMass(digestion.get(i), true));
            i++;
        }
        
        //return digestion;
    }
    
    public void modifyMasses(Modification[] modifications)
    {
        for (int m = 0; m < masses.size(); m++)
        {
            for (int n = 0; n < modifications.length; n++)
            {
                if ( modifications[n].getPosition().equals("protNterm") )
                {
                    if ( positions.get(m) == 0 ) // || positions.get(m) == 1 )
                        masses.set(m, masses.get(m) +  modifications[n].getDeltaMass());
                }
                else if ( modifications[n].getPosition().equals("pepNterm") )
                {
                    masses.set(m, masses.get(m) +  modifications[n].getDeltaMass());
                }
                else if ( modifications[n].getPosition().equals("pepCterm") )
                {
                    masses.set(m, masses.get(m) +  modifications[n].getDeltaMass());
                }
                else
                {
                    Character[] sites = modifications[n].getSites();
                    TreeMap<Character, Integer> composition = getComposition(digestion.get(m));
                    double deltaMass = 0.0;
                    
                    for (int s = 0; s < sites.length; s++)
                        deltaMass = composition.get(sites[s]) * modifications[n].getDeltaMass();
                    
                    masses.set(m, masses.get(m) + deltaMass);
                }
            }
        }
    }
    
    private double modifyMass(String sequence, double mass, char site, double deltaMass)
    {
        double newmass = mass;
        
        for (int i = 0; i < sequence.length(); i++)
        {
            if ( sequence.charAt(i) == site )
                newmass += deltaMass;
        }
        
        return newmass;
    }
    
    private TreeMap<Character, Integer> getComposition(String sequence)
    {
        TreeMap<Character, Integer> result = MassInfo.getCompositionTemplate();
        
        for (int i = 0; i < sequence.length(); i++)
        {
            Character aa = sequence.charAt(i);
            result.put(aa, result.get(aa) + 1);
        }
        
        return result;
    }
    
    public int getSize()
    {
        return digestion.size();
    }
    
    public List<String> getSequences()
    {
        return digestion;
    }
    
    public List<Short> getPositions()
    {
        return positions;
    }
    
    public List<Short> getLengths()
    {
        return lengths;
    }
    
    public List<Double> getMasses()
    {
        return masses;
    }
    
    public String getSequenceAt(int i)
    {
        return digestion.get(i);
    }
    
    public int getPositionAt(int i)
    {
        return (int)positions.get(i);
    }
    
    public int getLengthAt(int i)
    {
        return (int)lengths.get(i);
    }
    
    public double getMassAt(int i)
    {
        return masses.get(i);
    }
}
