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
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;
import java.nio.file.Paths;
import java.nio.file.Files;

import molecules.Protein;
import molecules.Enzyme;
import molecules.Modification;
import db_processing.Digestor;
import utils.ArrayFunctions;


public class PeptidesIndexation
{
    private final int maxIndexLength = 10000000;
    
    private String dataFile;
    private String indexFile;
    private RandomAccessFile dataReader;
    private Long[] fileOffsets;
    private Long[] entryIndices;
    protected String directoryName;
    public long gCounter;
    
    private List<String> tempFileNames = new ArrayList<String>();
    private Map<String, PeptideEntry> dbPepDict = new TreeMap<String, PeptideEntry>(); // what is better a HashMap or a TreeMap (sorted Map)?
    
    protected PeptidesIndexation()
    {}
    
    public PeptidesIndexation(String directory) throws FileNotFoundException, IOException
    {
        directoryName = directory;
        dataFile = directoryName + "/peptides";
        indexFile = directoryName + "/peptides.ind";
        RandomAccessFile reader = new RandomAccessFile(indexFile, "r");
        
        gCounter = reader.readLong();
        int c = reader.readInt();
        fileOffsets = new Long[c];
        entryIndices = new Long[c];
        
        for (int i = 0; i < c; i++)
        {
            fileOffsets[i] = reader.readLong();
            entryIndices[i] = reader.readLong();
        }
    }
    
    public PeptidesIndexation(String directory, ProteinsIndexation proteins, Enzyme enzyme, int missedCleavages, boolean semispecific, int minPepLen, int maxPepLen, Modification[] fixedModifications) throws Exception // enzymatic digestion
    {
        directoryName = directory;
        dataFile = directoryName + "/peptides";
        indexFile = directoryName + "/peptides.ind";
        
        for (int protIdx = 0; protIdx < proteins.getSize(); protIdx++)
        {
            if (protIdx%1000 == 999)
                System.out.println("Digesting " + (protIdx + 1) + "/" + proteins.getSize());
            
            Protein p = proteins.getProtein(protIdx);
            String protSeq = p.getSequence();
            Digestor myDigestor = new Digestor();
            
            if (enzyme == null)
            {
                throw new Exception("Enzyme name has not been specified.");
            }
            else
            {
                if (semispecific)
                    myDigestor.digest(protSeq, enzyme, missedCleavages, true, minPepLen, maxPepLen);
                else
                    myDigestor.digest(protSeq, enzyme, missedCleavages, false, minPepLen, maxPepLen);
                
                //myDigestor.modifyMasses(fixedModifications); // better not yet
            }
            
            if (dbPepDict.size() >= 4000000 )
                dumpTempFile();
            
            for (int i = myDigestor.getSize() - 1; i >= 0; i--)
            {
                if ( dbPepDict.containsKey(myDigestor.getSequenceAt(i)) )
                    dbPepDict.get(myDigestor.getSequenceAt(i)).addOccurrence(protIdx, myDigestor.getPositionAt(i));
                else
                    dbPepDict.put(myDigestor.getSequenceAt(i), new PeptideEntry(myDigestor.getLengthAt(i), myDigestor.getMassAt(i), protIdx, myDigestor.getPositionAt(i)));
            }
        }
        
        System.out.println("Finishing digestion.");

        writeFile(proteins, fixedModifications);  // Note: If there's any protein N-/C-term fixed modification, this subroutine ignores it. MS experiments don't have prot-terminal fixed modifications so is OK. 
        writeMemoryIndexFile();
    }
    
    private void dumpTempFile() throws FileNotFoundException, IOException
    {
        String fileName = directoryName + "/tmpPeptides" + tempFileNames.size();
        tempFileNames.add(fileName);
        RandomAccessFile writer = new RandomAccessFile(fileName, "rw");
        
        //for (Map.Entry<String, PeptideEntry> entry : dbPepDict.entrySet())
            //entry.getValue().write(writer);   
        String[] keys = dbPepDict.keySet().toArray(new String[0]);
        for (int p = 0; p < keys.length; p++)
            dbPepDict.get(keys[p]).write(writer);
        
        writer.close();
        dbPepDict.clear();
    }
    
    private void writeFile(ProteinsIndexation proteins, Modification[] fixedMods) throws FileNotFoundException, IOException
    {
        if (dbPepDict.size() > 0)
            dumpTempFile();
        
        dbPepDict = null;
        int n = tempFileNames.size();
        RandomAccessFile[] readers = new RandomAccessFile[n];
        PeptideEntry[] nextPeptides = new PeptideEntry[n];
        String[] sequences = new String[n];
        dataFile = directoryName + "/peptides";
        String tmpIndexFile = directoryName + "/peptides.ind_tmp";
        RandomAccessFile writer = new RandomAccessFile(dataFile, "rw");
        RandomAccessFile tmpIndexWriter = new RandomAccessFile(tmpIndexFile, "rw");
        
        for (int i = 0; i < n; i++)
        {
            readers[i] = new RandomAccessFile(tempFileNames.get(i), "r");
            nextPeptides[i] = PeptideEntry.read(readers[i]);
            sequences[i] = nextPeptides[i] != null ? nextPeptides[i].getSequence(proteins) : null;
        }
        
        Map<Character, Double> deltaMasses = Modification.arrayToDict(fixedMods);
        ArrayFunctions myFunctions = new ArrayFunctions();
        
        for (;;)
        {
            String first = myFunctions.getSmallestString(sequences);
            int[] indices = myFunctions.getIndicesOf(first, sequences);
            PeptideEntry p = PeptideEntry.merge(myFunctions.getSubarray(nextPeptides, indices));
            
            if (n > 1 && indices.length == sequences.length)
            {
                //System.out.println("Gotcha!");
                nextPeptides[0] = p;
                
                for (int i = 1; i < n; i++)
                {
                    nextPeptides[i] = PeptideEntry.read(readers[i]);
                    sequences[i] = nextPeptides[i] != null ? nextPeptides[i].getSequence(proteins) : null;
                }
                
                continue;
            }
            
            tmpIndexWriter.writeLong(writer.getFilePointer());
            p.applyFixedModifications(sequences[indices[0]], deltaMasses);
            p.write(writer);
            gCounter++;
            
            for (int i = 0; i < indices.length; i++){
                nextPeptides[indices[i]] = PeptideEntry.read(readers[indices[i]]);
                sequences[indices[i]] = nextPeptides[indices[i]] != null ? nextPeptides[indices[i]].getSequence(proteins) : null;
            }
            
            boolean finished = true;
            for (int i = 0; i < n; i++){
                if (nextPeptides[i] != null){
                    finished = false;
                    break;
                }
            }
            if (finished){
                break;
            }
        }
        
        writer.close();
        tmpIndexWriter.close();
        
        for (int i = 0; i < n; i++){
            readers[i].close();
        }
        for (int i = 0; i < n; i++){
            deleteFile(tempFileNames.get(i));
        }
        
        tempFileNames = null;
        int min = (int) Math.min(gCounter, maxIndexLength);
        fileOffsets = new Long[min];
        entryIndices = new Long[min];
        RandomAccessFile tmpIndexReader = new RandomAccessFile(tmpIndexFile, "r");
        
        if (min == gCounter)
        {
            for (int i = 0; i < min; i++)
            {
                fileOffsets[i] = tmpIndexReader.readLong();;
                entryIndices[i] = (long)i;
            }
        } 
        else
        {
            long previousEntryIndex = -1;
            for (int i = 0; i < min; i++)
            {
                entryIndices[i] = (gCounter * i) / min;
                long s = Long.MAX_VALUE;
                for (int j = 0; j < entryIndices[i] - previousEntryIndex; j++)
                    s = tmpIndexReader.readLong();
                
                fileOffsets[i] = s;
                previousEntryIndex = entryIndices[i];
            }
        }
        
        tmpIndexReader.close();
        deleteFile(tmpIndexFile);
        
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
    
    public void writeMemoryIndexFile() throws FileNotFoundException, IOException
    {
        RandomAccessFile writer = new RandomAccessFile(indexFile, "rw");
        writer.writeLong(gCounter);
        writer.writeLong(fileOffsets.length);
        
        for (int i = 0; i < fileOffsets.length; i++){
            writer.writeLong(fileOffsets[i]);
            writer.writeLong(entryIndices[i]);
        }
        writer.close();
    }
    
    public long getSize()
    {
        return gCounter;
    }
    
    public PeptideEntry getPeptide(long index)
    {
        try
        {
	        if (dataReader == null)
		        dataReader = new RandomAccessFile(dataFile, "r");
	
	        int i = ArrayFunctions.floorIndex(entryIndices, index);
	        dataReader.seek(fileOffsets[i]);
	        int o = (int) (index - entryIndices[i]);
            
	        for (int k = 0; k < o; k++)
		        PeptideEntry.read(dataReader);
	
	        return PeptideEntry.read(dataReader);
        }
        catch (FileNotFoundException e)
        {
            System.err.println(e.getMessage());
            return null;
        }
        catch (IOException e)
        {
            System.err.println("IOException while seeking file pointer.");
            System.err.println(e.getMessage());
            return null;
        }
	}

	public PeptideEntry[] getPeptides(Set<Long> indices)
    {
        PeptideEntry[] result = new PeptideEntry[indices.size()]; // not found indices will remain null in the array
		//indOut = new int[indices.size()];
        
        try
        {
	        if (dataReader == null)
		        dataReader = new RandomAccessFile(dataFile, "r");
	        else
		        dataReader.seek(0);
            
	        int count = 0;
            
	        for (long i = 0; count < indices.size() && i < gCounter; i++)
            {
		        PeptideEntry pep = PeptideEntry.read(dataReader);
		        if (indices.contains(i))
                {
			        result[count] = pep;
			        //indOut[count] = i;
			        count++;
		        }
	        }
            
	        return result;
        }
        catch (FileNotFoundException e)
        {
            System.err.println(e.getMessage());
            return null;
        }
        catch (IOException e)
        {
            System.err.println("IOException while seeking file pointer.");
            System.err.println(e.getMessage());
            return null;
        }
	}
    
    public void dispose() throws IOException
    {
		if (dataReader != null){
			dataReader.close();
			dataReader = null;
		}
		fileOffsets = null;
		entryIndices = null;
		if (tempFileNames != null){
			tempFileNames.clear();
			tempFileNames = null;
		}
		if (dbPepDict != null){
			dbPepDict.clear();
			dbPepDict = null;
		}
	}
    
    
    // JUST FOR TESTING
    public static void main (String[] args)
    {
        String directory = "/home/vsolis/Documents/GC/Javacode/XLID/data";
        String[] mainFastaFiles = {"/home/vsolis/Documents/GC/Javacode/XLID/data/UniProt_reviewed_human_canonical_201508.fasta"}; // human_short.fasta
        String conFastaFile = "/home/vsolis/Documents/GC/Javacode/XLID/data/Contaminants.fasta"; // contaminants_short.fasta
        String specialAAs = "KR";
        boolean includeContaminants = true;
        
        try
        {
            ProteinsIndexation proteins = new ProteinsIndexation(directory, mainFastaFiles, conFastaFile, specialAAs, includeContaminants);
            //ProteinsIndexation proteins = new ProteinsIndexation(directory);
            proteins.writeMemoryIndexFile();
            
            Enzyme trypsin = new Enzyme(1);
            Modification[] fixedModifications = { new Modification(3, true), new Modification(4, true) };
            PeptidesIndexation peptides = new PeptidesIndexation(directory, proteins, trypsin, 2, false, 7, 25, fixedModifications);
            
            System.out.println("TOTAL NUMBER OF PEPTIDES: " + peptides.getSize());
                       
            PeptideEntry expep = peptides.getPeptide(924405); // upper bound case 25
            if (expep != null)
            {
                System.out.println("sequence: " + expep.getSequence(proteins));
                System.out.println("aa before: " + expep.getResidueBefore(0, proteins));
                System.out.println("aa after: " + expep.getResidueAfter(0, proteins));
                System.out.println("length: " + expep.getLength());
                System.out.println("mass: " + expep.getMonoisotopicMass());
                for (int o = 0; o < expep.getNumberOfOccurrences(); o++)
                {
                    System.out.println("Protein index: " + expep.proteinIndexAt(o));
                    System.out.println("Protein offset: " + expep.proteinOffsetAt(o));
                }
            }
            
            Set<Long> myIndices = new TreeSet<Long>(Arrays.asList(new Long[]{2L, 10L, 7L, 3L, 934500L}));
            PeptideEntry[] ps = peptides.getPeptides(myIndices);
            for (int p = 0; p < ps.length; p++)
            {
                if (ps[p] != null)
                {
                    System.out.println("sequence: " + ps[p].getSequence(proteins));
                    System.out.println("aa before: " + ps[p].getResidueBefore(0, proteins));
                    System.out.println("aa after: " + ps[p].getResidueAfter(0, proteins));
                    System.out.println("length: " + ps[p].getLength());
                    System.out.println("mass: " + ps[p].getMonoisotopicMass());
                    for (int o = 0; o < ps[p].getNumberOfOccurrences(); o++)
                    {
                        System.out.println("Protein index: " + ps[p].proteinIndexAt(o));
                        System.out.println("Protein offset: " + ps[p].proteinOffsetAt(o));
                    }
                }
            }
        }
        catch (Exception e)
        {
            System.err.println(e.getMessage());
        }
    }
    
}
