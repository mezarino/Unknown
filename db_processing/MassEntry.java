/*
* LICENSE
*/

package db_processing;

/**
* 
*/

import java.io.*;
import java.util.List;
import java.util.ArrayList;

public class MassEntry
{
    public double monomass;
    public List<Long> peptideIndices;
    
    public MassEntry()
    {}
    
    public MassEntry(double mass, List<Long> indices)
    {
        monomass = mass;
        peptideIndices = indices;
    }
    
    public MassEntry(double mass, int npeptides)
    {
        peptideIndices = new ArrayList<Long>(npeptides);
        monomass = mass;
    }
    
    public static MassEntry read(RandomAccessFile reader)
    {
        try
        {
            double mass = reader.readDouble();
            int npeptides = reader.readInt();
            MassEntry result = new MassEntry(mass, npeptides);
            
            for (int i = 0; i < npeptides; i++)
                result.peptideIndices.add(reader.readLong());
            
            return result;
        }
        catch(IOException e)
        {
            return null;
        }
    }
    
    public void write(RandomAccessFile writer) throws IOException
    { 
        writer.writeDouble(monomass);
        writer.writeInt(peptideIndices.size());
        for (int i = 0; i < peptideIndices.size(); i++)
            writer.writeLong(peptideIndices.get(i));
    }
    
    public void addPeptideIndex(long index)
    {
        peptideIndices.add(index);
    }
    
    public void addPeptideIndices(List<Long> indices)
    {
        peptideIndices.addAll(indices);
    }
    
    public void extend(MassEntry[] modpeps)
    {
        for (int i = 0; i < modpeps.length; i++)
        {
            if (modpeps[i] != null)
                this.peptideIndices.addAll(modpeps[i].peptideIndices);
        }
    }
    
    public static MassEntry merge(MassEntry[] modpeps) // The user must make sure the entries have the same monomass
    {
        MassEntry result = null;
        
        for (int i = 0; i < modpeps.length; i++)
        {
            if (modpeps[i] == null)
                continue;
            else if (result == null)
                result = new MassEntry(modpeps[i].monomass, modpeps[i].peptideIndices);
            else
                result.peptideIndices.addAll(modpeps[i].peptideIndices);
        }
        
        return result;
    }
}
