package molecules;
/*
*
*
*/

import java.io.*;
import molecules.MassInfo;

public class Protein
{   
    private String id;
	private String description;
	private String sequence;
	private double weight;
	private boolean reverse;
	private boolean contaminant;

	public Protein(String id, String description, String sequence, boolean reverse, boolean contaminant)
    {
		this.id = id;
		this.description = description;
        this.sequence = sequence;
		weight = MassInfo.getWeight(sequence, true); // the 'true' is to specify if sequence is of aminoacids
		this.reverse = reverse;
		this.contaminant = contaminant;
	}

	public Protein(RandomAccessFile reader) throws IOException
    {
        //try
        //{
	    id = reader.readUTF();
	    description = reader.readUTF();
        sequence = reader.readUTF();
	    weight = reader.readDouble();
	    reverse = reader.readBoolean();
	    contaminant = reader.readBoolean();
        //}
        //catch (IOException e)
		//{ System.out.println(e.toString()); }
	}

	public void write(RandomAccessFile writer)
    {
        try
        {
		    writer.writeUTF(id);
            writer.writeUTF(description);
            writer.writeUTF(sequence);
            writer.writeDouble(weight);
            writer.writeBoolean(reverse);
            writer.writeBoolean(contaminant);
        }
        catch (IOException e)
		{ System.out.println(e.toString()); }
	}
    
    public String getId()
    {
		return id;
	}
    
    public void setId(String value)
    {
		id = value;
	}

	public String getDescription()
    {
		return description;
	}
    
    public void setDescription(String value)
    {
		description = value;
	}

	public String getSequence()
    {
		return sequence;
	}
    
    public void setSequence(String value)
    {
		sequence = value;
	}
    
    public double getWeight()
    {
		return weight;
	}
    
	public boolean isReverse()
    {
		return reverse;
	}

	public boolean isContaminant()
    {
		return contaminant;
	}

	public String getSubsequence(int offset, int length)
    {
		return sequence.substring(offset, offset+length);
	}
    
    public void digest() // args Enzyme enzyme, int missed, int maxLen
    {
        return;
    }

	public void dispose()
    {
		id = null;
		description = null;
        sequence = null;
	}

}
