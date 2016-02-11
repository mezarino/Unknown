/*
* LICENSE
*/

package db_processing;

/**
* Creates an indexed protein list from a FASTA file.
* The list is sorted by protein name. Each entry (i.e., row) in the list contains
* the following fields (i.e., columns) for a protein: name, description, sequence
* molecular weight, is_sequence_reversed and is_entry_contaminant.
*
* The list is stored in disc; the index is loaded in memory.
*
*/

import java.io.*;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashSet;
import java.nio.file.Paths;
import java.nio.file.Files;

import molecules.Protein;
import utils.ArrayFunctions;


public class ProteinsIndexation
{
    private int blockSize;
    private Long[] offsets;
    private String[] names;
    private RandomAccessFile dataReader;
    private RandomAccessFile indexReader;

    private String directoryName;
    private List<String> tempFileNames = new ArrayList<String>();
    private List<Integer> lengths2 = new ArrayList<Integer>();
    private List<String> names2 = new ArrayList<String>();
    private List<Long> offsets2 = new ArrayList<Long>();

    private String dataPath;
    private String longIndexPath;
    private String shortIndexPath;
    public int gCounter;
    
    public ProteinsIndexation()
    {}
    
    /**
    * Loads into memory the index file <em>proteins.inds</em> located in the
    * specify directory.
    *
    * The files <em>protein.ind</em> and <em>proteins</em> must be also in the
    * specify directory. Thus the first called to this <tt>Constructor</tt>
    * should be done after having created an instance using the second
    * <tt>Constructor</tt> in this class.
    * @param  directory directory path where the list and index files are located.
    * @throws FileNotFoundException if the list and index files are absent.
    * @throws IOException if an error occurs while reading the index file.
    */
    public ProteinsIndexation(String directory) throws FileNotFoundException, IOException
    {
        shortIndexPath = directory + "/proteins.inds";
        longIndexPath = directory + "/proteins.ind";
        dataPath = directory + "/proteins";
        
        //try
        //{
            RandomAccessFile reader = new RandomAccessFile(shortIndexPath, "r");
            gCounter = reader.readInt();
            blockSize = reader.readInt();
            int c = reader.readInt();
            offsets = new Long[c];
            names = new String[c];
            for (int i = 0; i < c; i++)
            {
                offsets[i] = reader.readLong();
                names[i] = reader.readUTF();
            }
            
            reader.close();
            
            dataReader = new RandomAccessFile(dataPath, "r");
            indexReader = new RandomAccessFile(longIndexPath, "r");
        //}
        //catch(FileNotFoundException e)
        //{
            //System.err.println("FileNotFoundException caught: " + e.getMessage());
        //}
        //catch(IOException e)
        //{
            //System.err.println("IOException caught: " + e.getMessage());
        //}
    }
    
    public ProteinsIndexation(String directory, String[] mainFastaFiles, String conFastaFile, String specialAAs, boolean includeContaminants) throws FileNotFoundException, IOException
    {
        String updateStatus = null;
        
        shortIndexPath = directory + "/proteins.inds";
        longIndexPath = directory + "/proteins.ind";
        dataPath = directory + "/proteins";
        this.directoryName = directory;
        HashSet<Character> keep = new HashSet<Character>(specialAAs.length());
        for (int i = 0; i < specialAAs.length(); i++)
            keep.add(specialAAs.charAt(i));
        
        //try
        //{
            RandomAccessFile writer = new RandomAccessFile(dataPath, "rw");
            for (int i = 0; i < mainFastaFiles.length; i++)
                addFastaFile(mainFastaFiles[i], writer, updateStatus, false, keep);
                
            if (includeContaminants)
                addFastaFile(conFastaFile, writer, updateStatus, true, keep);
            
            
            if (updateStatus != null)
                updateStatus = "Finishing proteins.";
            
            writer.close();
            writeIndexFile();
            writeMemoryIndexFile();
            
            dataReader = new RandomAccessFile(dataPath, "r");
            indexReader = new RandomAccessFile(longIndexPath, "r");
        //}
        //catch(FileNotFoundException e)
        //{
            //System.err.println("FileNotFoundException caught: " + e.getMessage());
        //}
        //catch(IOException e)
        //{
            //System.err.println("IOException caught: " + e.getMessage());
        //}
    }

    public void addFastaFile(String fastaFileName, RandomAccessFile writer, String updateStatus, boolean contaminant, HashSet<Character> keep) throws FileNotFoundException, IOException
    {
        // Iterate through fasta file
        String id = null, description = null, sequence = "", lastline="";
        long offset = 0;
        
        //try
        //{
            RandomAccessFile reader = new RandomAccessFile(fastaFileName, "r");
            String line = reader.readLine();
            int counter = 0;
            
            while (line != null)
            {   
                line = line.trim();
                
                if (line.length() == 0)
                {//do nothing and pass to next line;
                    //System.out.println("Blank line found!");
                }
                else if (line.startsWith(">"))
                {
                    counter++;
                    
                    if (sequence.length() > 0)
                    {
                        addRecord(id, description, sequence, writer, false, contaminant);
                        
                        String revSequence = revertSequence(sequence, keep);
                        String revId = "REV__" + id;
                        addRecord(revId, description, revSequence, writer, true, contaminant);
                    }
                    
                    id = line.split(" ")[0].substring(1);
                    description = line.substring(1);
                    sequence = "";
                    
                    if (contaminant)
                        id = "CON__" + id;
                    
                    if (updateStatus != null && counter%10000 == 0){
                        updateStatus = "Parsing protein " + counter + ".";
                        System.out.println(updateStatus);
                    }
                }
                else
                {
                    sequence = sequence.concat(line);
                }
                
                line = reader.readLine();
                lastline = line;
            }
            
            // Process last sequence in fasta file
            if (sequence.length() > 0)
            {
                addRecord(id, description, sequence, writer, false, contaminant);
                
                String revSequence = revertSequence(sequence, keep);
                String revId = "REV__" + id;
                addRecord(revId, description, revSequence, writer, true, contaminant);
            }
            
            reader.close();
        //}
        //catch(FileNotFoundException e)
        //{
            //System.err.println("FileNotFoundException caught: " + e.getMessage());
            //System.err.println("The file " + fastaFileName + " could not be read.");
        //}
        //catch(IOException e)
        //{
            //System.err.println("IOException caught: " + e.getMessage());
            //System.err.println("Something went wrong while reading " + fastaFileName + " at sequence: " + id);
        //}  
    }

    private void addRecord(String id, String description, String sequence, RandomAccessFile writer, boolean reverse, boolean contaminant) throws FileNotFoundException, IOException
    {
        if (offsets2.size() > 1000000){
            dumpTempFile();
        }
        
        if (sequence.endsWith("*")){
            sequence = sequence.substring(0, sequence.length() - 1);
        }
        
        Protein p = new Protein(id, description, sequence, reverse, contaminant);
        //try {
        offsets2.add(writer.getFilePointer());
        //}
        //catch (IOException e) {System.err.println(e.getMessage()); } 
        p.write(writer);
        names2.add(p.getId());
        lengths2.add(sequence.length());
    }
    
    private void dumpTempFile() throws FileNotFoundException, IOException
    {
        String filename = directoryName + "/tmpProtInd" + tempFileNames.size();
        tempFileNames.add(filename);
        
        ArrayFunctions myFunctions = new ArrayFunctions();
        String[] names2Array = new String[names2.size()];
        names2.toArray(names2Array);
        Integer[] idxs = myFunctions.getSortedIndexes(names2Array);
        names2Array = ArrayFunctions.getSubarray(names2Array, idxs);
        
        //try
        //{
            RandomAccessFile writer = new RandomAccessFile(filename, "rw");
            
            for (int i = 0; i < idxs.length; i++){
                if (i == 0 || !names2Array[i].equals(names2Array[i - 1])){
                    writer.writeUTF(names2Array[i]);
                    writer.writeInt(lengths2.get(idxs[i]));
                    writer.writeLong(offsets2.get(idxs[i]));
                    gCounter++;
                }
            }
            
            writer.close();
        //}
        //catch(FileNotFoundException e)
        //{
            //System.err.println("FileNotFoundException caught: " + e.getMessage());
        //}
        //catch(IOException e)
        //{
            //System.err.println("IOException caught: " + e.getMessage());
        //}
        
        names2.clear();
        lengths2.clear();
        offsets2.clear();
    }
    
    private static String revertSequence(String sequence, HashSet<Character> keep)
    {
        char[] rev = new char[sequence.length()];
        
        for (int i = sequence.length(); i > 0; i--){
            rev[i-1] = sequence.charAt(i - 1);
        }
        
        // Move cleavage sites
        for (int i = 1; i < sequence.length(); i++){
            if (keep.contains(rev[i])){
                char c = rev[i - 1];
                rev[i - 1] = rev[i];
                rev[i] = c;
            }
        }
        
        return new String(rev);
    }
    
    
    private void writeIndexFile() throws FileNotFoundException, IOException
    {
        if (offsets2.size() > 0)
            dumpTempFile();
        
        blockSize = gCounter/1000000 + 1;
        offsets2 = null;
        names2 = null;
        lengths2 = null;
        
        //try
        //{
            int n = tempFileNames.size();
            RandomAccessFile[] readers = new RandomAccessFile[n];
            long[] nextFilePos = new long[n];
            String[] nextNames = new String[n];
            int[] nextLengths = new int[n];
            RandomAccessFile writer = new RandomAccessFile(longIndexPath, "rw");
        
            for (int i = 0; i < n; i++)
            {
                readers[i] = new RandomAccessFile(tempFileNames.get(i), "r");
            
                try
                {
                    nextNames[i] = readers[i].readUTF();
                    nextLengths[i] = readers[i].readInt();
                    nextFilePos[i] = readers[i].readLong();
                } 
                catch (IOException e)
                {
                    nextNames[i] = null;
                    nextLengths[i] = Integer.MAX_VALUE;
                    nextFilePos[i] = Long.MAX_VALUE;
                    System.err.println("IOException caught: " + e.getMessage());
                }
            }
            
            List<Long> blockFilePos = new ArrayList<Long>();
            List<Integer> blockLengths = new ArrayList<Integer>();
            List<String> blockNames = new ArrayList<String>();
            List<Long> filePosShort2 = new ArrayList<Long>();
            List<String> namesShort2 = new ArrayList<String>();
            
            for (;;)
            {
                int smallestInd = getSmallestInd(nextNames);
                for (int i = 0; i < n; i++){
                    if (i == smallestInd){
                        continue;
                    }
                    if (nextNames[i] == null){
                        continue;
                    }
                    if (nextNames[i].equals(nextNames[smallestInd])){
                        try
                        {
                            nextNames[i] = readers[i].readUTF();
                            nextLengths[i] = readers[i].readInt();
                            nextFilePos[i] = readers[i].readLong();
                        } 
                        catch (IOException e)
                        {
                            nextNames[i] = null;
                            nextLengths[i] = Integer.MAX_VALUE;
                            nextFilePos[i] = Long.MAX_VALUE;
                            System.err.println("IOException caught: " + e.getMessage());
                        }
                    }
                }
                
                blockFilePos.add(nextFilePos[smallestInd]);
                blockLengths.add(nextLengths[smallestInd]);
                blockNames.add(nextNames[smallestInd]);
                if (blockFilePos.size() == blockSize){
                    writeBlock(blockFilePos, blockLengths, blockNames, writer, filePosShort2, namesShort2);
                }
            
                try
                {
                    nextNames[smallestInd] = readers[smallestInd].readUTF();
                    nextLengths[smallestInd] = readers[smallestInd].readInt();
                    nextFilePos[smallestInd] = readers[smallestInd].readLong();
                } 
                catch (IOException e) // end of temp file.
                {
                    nextNames[smallestInd] = null;
                    nextLengths[smallestInd] = Integer.MAX_VALUE;
                    nextFilePos[smallestInd] = Long.MAX_VALUE;
                    //System.err.println("IOException caught: " + e.getMessage());
                }
                
                boolean finished = true;
                for (int i = 0; i < n; i++)
                {
                    if (nextNames[i] != null)
                    {
                        finished = false;
                        break;
                    }
                }
                
                if (finished) // end of all temp files
                    break;
            }
            
            if (blockFilePos.size() > 0){
                writeBlock(blockFilePos, blockLengths, blockNames, writer, filePosShort2, namesShort2);
            }
            
            offsets = new Long[filePosShort2.size()];
            filePosShort2.toArray(offsets);
            names = new String[namesShort2.size()];
            namesShort2.toArray(names);
            writer.close();
            
            for (int i = 0; i < tempFileNames.size(); i++){
                readers[i].close();
            }
            for (int i = 0; i < tempFileNames.size(); i++){
                deleteFile(tempFileNames.get(i));
            }
            tempFileNames = null;
        //}
        //catch(FileNotFoundException e)
        //{
            //System.err.println("FileNotFoundException caught: " + e.getMessage());
        //}
        //catch(IOException e)
        //{
            //System.err.println("IOException caught: " + e.getMessage());
        //}
    }

    private static void writeBlock(List<Long> blockFilePos, List<Integer> blockLengths, List<String> blockNames, RandomAccessFile writer, List<Long> filePosShort2, List<String> namesShort2) throws IOException
    {
        //try
        //{
            filePosShort2.add(writer.getFilePointer());
            
            for (int l = 0; l < blockFilePos.size(); l++)
                writer.writeLong(blockFilePos.get(l));
        
            for (int i = 0; i < blockLengths.size(); i++)
                writer.writeInt(blockLengths.get(i));
        
            for (int s = 0; s < blockNames.size(); s++)
                writer.writeUTF(blockNames.get(s));
        //}
        //catch(IOException e)
        //{
            //System.err.println("IOException caught: " + e.getMessage());
        //}
        
        namesShort2.add(blockNames.get(0));
        blockFilePos.clear();
        blockLengths.clear();
        blockNames.clear();
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
    
    private static int getSmallestInd(String[] strings)
    {
        int ind = -1;
        String smallest = null;
        for (int i = 0; i < strings.length; i++){
            if (strings[i] == null){
                continue;
            }
            if (smallest == null){
                smallest = strings[i];
                ind = i;
            } else{
                if (strings[i].compareTo(smallest) < 0){
                    smallest = strings[i];
                    ind = i;
                }
            }
        }
        
        return ind;
    }
    
    public void writeMemoryIndexFile() throws FileNotFoundException, IOException
    {
        if (shortIndexPath == null)
            return;
        
        //try
        //{
            RandomAccessFile writer = new RandomAccessFile(shortIndexPath, "rw");
            writer.writeInt(gCounter);
            writer.writeInt(blockSize);
            writer.writeInt(offsets.length);
            for (int i = 0; i < offsets.length; i++){
                writer.writeLong(offsets[i]);
                writer.writeUTF(names[i]);
            }
            
            writer.close();
        //}
        //catch(FileNotFoundException e)
        //{
            //System.err.println("FileNotFoundException caught: " + e.getMessage());
        //}
        //catch(IOException e)
        //{
            //System.err.println("IOException caught: " + e.getMessage());
        //}
    }
    
    public int getSize()
    {
        return gCounter;
    }
    
    private long getFilePos(int index) throws IOException
    {
        int block = index/blockSize;
        int off = index%blockSize;
        
        indexReader.seek(offsets[block] + 8*off);
        return indexReader.readLong();
    }

    public int getLength(int index) throws IOException
    {
        int block = index/blockSize;
        int off = index%blockSize;
        
        indexReader.seek(offsets[block] + 8*blockSize + 4*off);
        return indexReader.readInt();
    }
    
    public String getName(int index) throws IOException
    {
        int block = index/blockSize;
        int off = index%blockSize;
        
        indexReader.seek(offsets[block] + 12*blockSize);
        for (int i = 0; i < off; i++){
            indexReader.readUTF();
        }
        
        return indexReader.readUTF();
    }

    public int getIndex(String name) throws Exception, IOException
    {
        int a = Arrays.binarySearch(names, name);
            if (a < 0)
                a = -2 - a;
        
        indexReader.seek(offsets[a] + 12*blockSize);
        for (int i = 0; i < blockSize; i++)
        {
            String s = indexReader.readUTF();
            if (s.equals(name))
                return a*blockSize + i;
        }
        throw new Exception("Protein " + name + " does not exist in the fasta file.");
    }
    
    public Protein getProtein(String name) 
    {
        try
        {
            return getProtein(getIndex(name));
        }
        catch (Exception e)
        {
            System.err.println(e.getMessage());
            return null;
        }
    }
    
    public Protein getProtein(int index)
    {
        try
        {
            dataReader.seek(getFilePos(index));
            return new Protein(dataReader);
        }
        catch (IOException e)
        {
            System.err.println(e.getMessage());
            return null;
        }
    }
    
    public char getResidueAtProteinPos(int index, int pos)
    {
        try
        {
            if ( pos < 0 || pos >= getLength(index) )
                return ' ';
            
            dataReader.seek(getFilePos(index));
            dataReader.readUTF(); // name
            dataReader.readUTF(); // description
            String result = dataReader.readUTF();// We are now at sequence
            return result.charAt(pos);
        }
        catch (IOException e)
        {
            System.err.println(e.getMessage());
            return ' ';
        }
    }
    
    public String getSubsequenceFromProtein(int index, int pos, int length)
    {
        try
        {
            int protLen = getLength(index);
            
            if ( pos < 0 )
                pos = 0;
            if ( pos >=  protLen)
                return "";
            
            dataReader.seek(getFilePos(index));
            dataReader.readUTF(); // name
            dataReader.readUTF(); // description
            String result = dataReader.readUTF();// We are now at sequence
            
            if ( pos + length > protLen )
                length = protLen - pos;
            
            return result.substring(pos, pos+length);
        }
        catch (IOException e)
        {
            System.err.println(e.getMessage());
            return "";
        }
    
    }
    
    public void dispose() throws IOException
    {
		if (dataReader != null){
			dataReader.close();
			dataReader = null;
		}
		if (indexReader != null){
			indexReader.close();
			indexReader = null;
		}
        
		offsets = null;
		names = null;
        
		if (tempFileNames != null){
			tempFileNames.clear();
			tempFileNames = null;
		}
		if (lengths2 != null){
			lengths2.clear();
			lengths2 = null;
		}
		if (names2 != null){
			names2.clear();
			names2 = null;
		}
		if (offsets2 != null){
			offsets2.clear();
			offsets2 = null;
		}
	}
    
    
    // JUST FOR TESTING
    public static void main(String[] args)
    {
        String directory = "/home/vsolis/Documents/GC/Javacode/XLID/data/test_short";
        String[] mainFastaFiles = {"/home/vsolis/Documents/GC/Javacode/XLID/data/human_short.fasta"}; // UniProt_reviewed_human_canonical_201508.fasta
        String conFastaFile = "/home/vsolis/Documents/GC/Javacode/XLID/data/contaminants_short.fasta"; // Contaminants.fasta
        String specialAAs = "KR";
        boolean includeContaminants = true;
        
        try
        {
            ProteinsIndexation myIndexation = new ProteinsIndexation(directory, mainFastaFiles, conFastaFile, specialAAs, includeContaminants);
            myIndexation.writeMemoryIndexFile();
            ProteinsIndexation myIndexation2 = new ProteinsIndexation(directory);
            System.out.println("Total number of proteins: " + myIndexation2.getSize());
            
            String nameAtIdx = myIndexation2.getName(10);
            System.out.println(nameAtIdx);
            System.out.println(myIndexation2.getIndex(nameAtIdx));
            System.out.println(myIndexation2.getLength(10));
            System.out.println(myIndexation2.getFilePos(10));
            int idxOfName = myIndexation2.getIndex("sp|P31947|1433S_HUMAN");
            System.out.println(idxOfName);
            System.out.println(myIndexation2.getName(idxOfName));
            System.out.println(myIndexation2.getLength(idxOfName));
            System.out.println(myIndexation2.getFilePos(idxOfName));
            
            for (int i = 0; i < myIndexation2.getSize(); i++)
                System.out.println(myIndexation2.getName(i));
        }
        catch(Exception e)
        {
            System.out.println("There was an error while loading/reading/writing files.");
            System.err.println(e.toString());
            System.exit(1);
        }
    }

}
