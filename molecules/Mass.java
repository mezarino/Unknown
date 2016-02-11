/*
* LICENSE
*/

package molecules;

/**
* 
*/

import java.io.*;
import java.util.Arrays;


public class Mass
{
    public double monomass;
    public long peptideIndex;
    public int nMods = 0;
    private Short[] modPositions;
    private Short[] modCodes;
    
    public Mass()
    {}
    
    public Mass(double mass, long idx, int nMods)
    {
        monomass = mass;
        peptideIndex = idx;
        this.nMods = nMods;
    }

    public Mass(double mass, long idx, int pos, int code)
    {
        monomass = mass;
        peptideIndex = idx;
        modPositions = new Short[]{(short)pos};
        modCodes = new Short[]{(short)code};
    }
    
    public Mass modify(double deltaMass, int pos, int code)
    {
        Mass newmass = new Mass();
        newmass.monomass = this.monomass + deltaMass;
        newmass.modPositions = Arrays.copyOf(this.modPositions, this.modPositions.length + 1);
        newmass.modPositions[this.modPositions.length] = (short)pos;
        newmass.modCodes = Arrays.copyOf(this.modCodes, this.modCodes.length + 1);
        newmass.modCodes[this.modCodes.length] = (short)code;
        
        return newmass;
    }
    
    public int getNumberOfMods()
    {
        if (modPositions != null)
            nMods = modPositions.length;
        
        return nMods;
    }
    
    public void writeTemp(RandomAccessFile writer) throws IOException // this is use in MassessIndexation.java to write temporary files
    { 
        writer.writeDouble(monomass);
        writer.writeInt(1);
        writer.writeLong(peptideIndex);
    }
    
    public void writeShort(RandomAccessFile writer) throws IOException // this is use in MassessIndexation.java to write temporary files
    { 
        writer.writeDouble(monomass);
        writer.writeLong(peptideIndex);
        writer.writeInt(nMods);
    }
    
    public void writeFull(RandomAccessFile writer) throws IOException // this is use in MassessIndexation.java to write temporary files
    { 
        writer.writeDouble(monomass);
        writer.writeLong(peptideIndex);
        writer.writeInt(nMods);
        for (int i = 0; i < nMods; i++)
        {
            writer.writeShort(modPositions[i]);
            writer.writeShort(modCodes[i]);
        }
    }
}
