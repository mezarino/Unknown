package utils;
/*
*
*
*/

import java.util.Comparator;
import java.util.Arrays;

public class ArrayIndexComparator<T extends Comparable<T>> implements Comparator<Integer>
{   
    private final T[] mainArray;
    
    public ArrayIndexComparator(T[] array)
    {
        mainArray = array;
    }
    
    public Integer[] createIndexArray()
    {
        Integer[] indexes = new Integer[mainArray.length];
        
        for (int i = 0; i < mainArray.length; i++)
        {
            indexes[i] = i;
        }
        
        return indexes;
    }
    
    public int compare(Integer index1, Integer index2)
    {
        return mainArray[index1].compareTo(mainArray[index2]);
    }
    
    // Just for debugging
    public static void main(String[] args)
    {
        String[] a = {"pep", "set", "ara", "xav", "bot"};
        ArrayIndexComparator<String> comparator = new ArrayIndexComparator<String>(a);
        Integer[] idxs = comparator.createIndexArray();
        Arrays.sort(idxs, comparator);
        for (int i = 0; i < idxs.length; i++)
            System.out.println("old: " + a[i] + "  new: " + a[idxs[i]]);
        
        Double[] b = {12.4, 100.0, 34.7, 0.98, 72.0};
        ArrayIndexComparator<Double> comparator2 = new ArrayIndexComparator<Double>(b);
        Integer[] idxs2 = comparator2.createIndexArray();
        Arrays.sort(idxs2, comparator2);
        for (int i = 0; i < idxs2.length; i++)
            System.out.println("old: " + b[i] + "  new: " + b[idxs2[i]]);
    }
}


