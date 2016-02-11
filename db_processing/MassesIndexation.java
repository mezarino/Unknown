/*
* LICENSE
*/

package db_processing;

/**
* 
*/

import java.io.*;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Map;
import java.nio.file.Paths;
import java.nio.file.Files;

import utils.ArrayFunctions;
import molecules.Mass;
import molecules.Modification;

public class MassesIndexation
{
    private final int maxIndexLength = 8000000;
    private String dataFile;
    private String indexFile;
    private RandomAccessFile dataReader;
    private Long[] filePointers;
    private Double[] masses;
    protected String directoryName;
    public long gCounter;
    private List<String> tempFileNames = new ArrayList<String>();
    private List<Long> ind = new ArrayList<Long>();
    private List<Double> m = new ArrayList<Double>();
    private int min;
    
    public MassesIndexation()
    {}
    
    public MassesIndexation(String directory) throws FileNotFoundException, IOException
    {
        directoryName = directory;
        dataFile = directoryName + "/modPeptides";
        indexFile = directoryName + "/modPeptides.ind";
        
        RandomAccessFile reader = new RandomAccessFile(indexFile, "r");
        gCounter = reader.readLong();
        int c = reader.readInt();
        min = reader.readInt();
        filePointers = new Long[c];
        masses = new Double[c];
        
        for (int i = 0; i < c; i++)
        {
            filePointers[i] = reader.readLong();
            masses[i] = reader.readDouble();
        }
    }
    
    public MassesIndexation(String directory, PeptidesIndexation peptides, ProteinsIndexation proteins, Modification[] varModifications, int maxMods, double maxPepMass) throws FileNotFoundException, IOException
    {
        directoryName = directory;
        dataFile = directoryName + "/modPeptides";
        indexFile = directoryName + "/modPeptides.ind";
        
        Map<Character, Double> deltaMasses = Modification.arrayToDict(varModifications);
        
        for (long pepIdx = 0; pepIdx < peptides.getSize(); pepIdx++)
        {
            if (pepIdx%10000 == 9999){
				System.out.println("Modifying peptide " + (pepIdx + 1) + "/" + peptides.getSize());
			}
            
            PeptideEntry p = peptides.getPeptide(pepIdx);
            String pepSeq = p.getSequence(proteins);
            double origmass = p.getMonoisotopicMass();
            
            if ( p.isAtProteinNterm(proteins) )
                pepSeq = "-" + pepSeq;
            if ( p.isAtProteinCterm(proteins) )
                pepSeq = pepSeq + "_";
            
            Mass[] entries = p.applyVariableModifications(pepSeq, origmass, deltaMasses, pepIdx, maxMods);
            
            if (ind.size() >= 10000000)
		        dumpTempFile();
            
            for (int i = 0; i < entries.length; i++)
            {
                ind.add(entries[i].peptideIndex);
	            m.add(entries[i].monomass);
            }
        }
    }
    
    public void writeFile() throws FileNotFoundException, IOException
    {
		if (m.size() > 0){
			dumpTempFile();
		}
        
        m.clear();
        ind.clear();
		m = null;
		ind = null;
        
		int n = tempFileNames.size();
		RandomAccessFile[] readers = new RandomAccessFile[n];
		MassEntry[] nextPeptides = new MassEntry[n];
		dataFile = directoryName + "/modPeptides";
		indexFile = directoryName + "/modPeptides.ind";
		String tmpIndexFile = directoryName + "/modPeptides.ind_tmp";
		gCounter = 0;
		RandomAccessFile writer = new RandomAccessFile(dataFile, "rw");
		RandomAccessFile tmpIndexWriter = new RandomAccessFile(tmpIndexFile, "rw");
        
		for (int i = 0; i < n; i++){
			readers[i] = new RandomAccessFile(tempFileNames.get(i), "r");
			nextPeptides[i] = MassEntry.read(readers[i]);
		}
		
		for (;;)
        {
			double smallest = getSmallest(nextPeptides);
			int[] indices = getIndices(smallest, nextPeptides);
			//MassEntry[] toBeWritten = ArrayFunctions.getSubarray(nextPeptides, indices);
            MassEntry mp = MassEntry.merge(ArrayFunctions.getSubarray(nextPeptides, indices));
            
            if (n > 1 && indices.length == nextPeptides.length)
            {
                //System.out.println("Gotcha!");
                nextPeptides[0] = mp;
                
                for (int i = 1; i < n; i++)
                    nextPeptides[i] = MassEntry.read(readers[i]);
                
                continue;
            }
            
			//for (int i = 0; i < toBeWritten.length; i++)
            //{
				//tmpIndexWriter.writeLong(writer.getFilePointer());
				//tmpIndexWriter.writeDouble(toBeWritten[i].monomass);
				//toBeWritten[i].write(writer);
				//gCounter++;
			//}
            
            tmpIndexWriter.writeLong(writer.getFilePointer());
            tmpIndexWriter.writeDouble(mp.monomass);
            mp.write(writer);
            gCounter++;
            
			for(int i = 0; i < indices.length; i++){
				nextPeptides[indices[i]] = MassEntry.read(readers[indices[i]]);
			}
            
			boolean finished = true;
			for (int i = 0; i < n; i++)
            {
				if (nextPeptides[i] != null)
                {
					finished = false;
					break;
				}
			}
            
			if (finished)
				break;
		}
        
		writer.close();
		tmpIndexWriter.close();
        
		for (int i = 0; i < tempFileNames.size(); i++){
			readers[i].close();
		}
        for (int i = 0; i < tempFileNames.size(); i++){
			deleteFile(tempFileNames.get(i));
		}
        
		tempFileNames = null;
		min = (int) Math.min(gCounter, maxIndexLength);
		filePointers = new Long[min];
		masses = new Double[min];
		RandomAccessFile tmpIndexReader = new RandomAccessFile(tmpIndexFile, "r");
        
		if (min == gCounter)
        {
			for (int i = 0; i < min; i++)
            {
				long s = tmpIndexReader.readLong();
				double mass = tmpIndexReader.readDouble();
				filePointers[i] = s;
				masses[i] = mass;
			}
		} 
        else
        {
			long previousEntryIndex = -1;
			for (int i = 0; i < min; i++){
				long s = Long.MAX_VALUE;
				double mass = Double.NaN;
                long entryIndex = (gCounter * i) / min;
                
				for (int j = 0; j < entryIndex - previousEntryIndex; j++)
                {
					s = tmpIndexReader.readLong();
					mass = tmpIndexReader.readDouble();
				}
				filePointers[i] = s;
				masses[i] = mass;
				previousEntryIndex = entryIndex;
			}
		}
        
		tmpIndexReader.close();
		deleteFile(tmpIndexFile);
	}
     
    public void writeMemoryIndexFile() throws FileNotFoundException, IOException
    {
        RandomAccessFile writer = new RandomAccessFile (indexFile, "rw");
        writer.writeLong(gCounter);
        writer.writeInt(filePointers.length);
        writer.writeInt(min);
        
        for (int i = 0; i < filePointers.length; i++)
        {
            writer.writeLong(filePointers[i]);
            writer.writeDouble(masses[i]);
        }
        
        writer.close();   
    }
    
    private void dumpTempFile() throws FileNotFoundException, IOException
    {
		Double[] massList = new Double[m.size()];
		for (int i = 0; i < massList.length; i++){
			massList[i] = m.get(i);
		}
        
		Integer[] o = ArrayFunctions.getSortedIndexes(massList);
		String fileName = directoryName + "/tmpModPeptides" + tempFileNames.size();
		tempFileNames.add(fileName);
		RandomAccessFile writer = new RandomAccessFile(fileName, "rw");
		
        for (int i = 0; i < massList.length; i++){
			Mass p = new Mass(m.get(o[i]), ind.get(o[i]), 0);
			p.writeTemp(writer);
		}
        
		writer.close();
		ind.clear();
		m.clear();
	}
    
    public void addEntry(MassEntry entry)
    {
        try
        {
	        if (ind.size() >= 10000000)
		        dumpTempFile();
        }
        catch (Exception e)
        {
            System.err.println(e.getMessage());
            
        }
	    
	    ind.addAll(entry.peptideIndices);
	    m.add(entry.monomass);
    }
    
    public void addEntry(Mass entry)
    {
        try
        {
	        if (ind.size() >= 10000000)
		        dumpTempFile();
        }
        catch (Exception e)
        {
            System.err.println(e.getMessage());
            
        }
	    
	    ind.add(entry.peptideIndex);
	    m.add(entry.monomass);
    }
    
    private static double getSmallest(MassEntry[] peptides)
    {
		double min = Double.NaN;
        
		for (int i = 0; i < peptides.length; i++)
        {
			if (peptides[i] == null)
				continue;
			
			if (Double.isNaN(min))
            {
				min = peptides[i].monomass;
				continue;
			}
            
			if (min > peptides[i].monomass)
				min = peptides[i].monomass;
		}
        
		return min;
	}
    
    private static int[] getIndices(double smallest, MassEntry[] peptides)
    {
		int[] result = new int[peptides.length];
        int size = 0;
        
		for (int i = 0; i < peptides.length; i++)
        {
			if (peptides[i] != null && peptides[i].monomass == smallest)
            {
				result[size] = i;
                size++;
			}
		}
        
		return Arrays.copyOf(result, size);
	}
    
    private static void deleteFile(String fileName)
    {   
        try
        {
            Files.delete(Paths.get(fileName));
        }
        catch (Exception e)
        {
            System.err.println("Not possible to delete file " + fileName);
            System.err.println(e);
        }
    }
    
    public void dispose() throws IOException
    {
	    if (dataReader != null)
        {
		    dataReader.close();
		    dataReader = null;
	    }
        
	    filePointers = null;
	    masses = null;
        
	    if(tempFileNames != null)
        {
		    tempFileNames.clear();
		    tempFileNames = null;
	    }
        
	    if(ind != null)
        {
		    ind.clear();
		    ind = null;
	    }
        
	    if(m != null){
		    m.clear();
		    m = null;
	    }
    }

    
    
}
