package utils;
/*
*
*
*/
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.TreeMap;
import java.lang.reflect.Array;


public class ArrayFunctions
{
    public ArrayFunctions()
    {}
    
    public static <T> T[] getSubarray(T[] array, int[] indices)
    {
		T[] result = Arrays.copyOf(array, indices.length);
		for (int i = 0; i < result.length; i++)
			result[i] = array[indices[i]];
        
		return result;
    }
    
    public static <T> T[] getSubarray(T[] array, Integer[] indices)
    {
		T[] result = Arrays.copyOf(array, indices.length);
		for (int i = 0; i < result.length; i++)
			result[i] = array[indices[i]];
        
		return result;
    }
    
    public static int[] getIntegersArray(int from, int to)
    {
        int[] result = new int[to-from];
        for (int i = 0; i < result.length; i++)
            result[i] = from + i;
        
        return result;   
    }
    
    public static <T extends Comparable<T>> Integer[] getSortedIndexes(T[] array)
    {
        ArrayIndexComparator<T> comparator = new ArrayIndexComparator<T>(array);
        Integer[] idxs = comparator.createIndexArray();
        Arrays.sort(idxs, comparator);
        
        return idxs;
    }
    
    public static String concatenateItems(String[] array)
    {
         String result = "";
        
        for (int i = 0; i < array.length; i++)
        {
            result = result.concat(array[i]);
        }
        
        return result;
    }
    
    public static String concatenateItems(List<String> list)
    {
        String result = "";
        
        for (int i = 0; i < list.size(); i++)
        {
            result = result.concat(list.get(i));
        }
        
        return result;
    }
    
    public static <T extends Comparable<T>> T getSmallestString(T[] array)
    {
        T first = null;
        
        for (int i = 0; i < array.length; i++)
        {
            if (array[i] == null)
            {
                continue;
            }
            if (first == null)
            {
                first = array[i];
                continue;
            }
            if (first.compareTo(array[i]) > 0)
            {
                first = array[i];
            }
        }
        
        return first;
    }
    
    public static int[] getIndicesOf(Object object, Object[] objects)
    {
		int[] result = new int[objects.length];
        int size = 0;
        
		for (int i = 0; i < objects.length; i++){
			if (objects[i] != null && objects[i].equals(object))
            {
				result[size] = i;
                size++;
			}
		}
        
		return Arrays.copyOf(result, size);
	}
    
    public static <T extends Comparable<T>> int floorIndex(T[] array, T value)
    {
		int n = array.length;
        
		if (n == 0)
			return -1;
		else if (value.compareTo(array[n - 1]) > 0)
			return n - 1;
		else if (value.compareTo(array[0]) < 0)
			return -1;
        
		int a = Arrays.binarySearch(array, value);
        
		if (a >= 0)
        {
			while (a < array.length - 1 && array[a + 1].equals(array[a]))
				a++;
		
			return a;
		}
        
		return -2 - a;
    }

	public static <T extends Comparable<T>> int ceilIndex(T[] array, T value)
    {
		int n = array.length;
        
        if (n == 0)
			return -1;
		else if (value.compareTo(array[n - 1]) > 0)
			return -1;
		else if (value.compareTo(array[0]) < 0)
			return 0;
        
		int a = Arrays.binarySearch(array, value);
        
		if (a >= 0)
        {
			while (a > 0 && array[a - 1].equals(array[a]))
				a--;
			
			return a;
		}
        
		return -1 - a;
	}
    
    public static TreeMap<Character, Integer> getCharComposition(String string)
    {
        TreeMap<Character, Integer> result = new TreeMap<Character, Integer>();
        
        for (int i = 0; i < string.length(); i++)
        {
            Character a = string.charAt(i);
            result.put(a, result.get(a) + 1);
        }
        
        return result;
    }
    
    // Combinations (of characters)
    private static String[] combine(String[] strings, Character[] s)
    {
        List<String> result = new ArrayList<String>(s.length);
        int wl = strings[0].length(); //current words' length

        for (int i = 0; i < strings.length; i++)
        {
            for (int j = Arrays.binarySearch(s, strings[i].charAt(wl-1)); j < s.length; j++)
            {
                result.add(strings[i] + s[j]);
            }
        }

        return result.toArray(new String[result.size()]);
    }

    public static String[] getCombinations(int l, Character[] s)
    {
        if (l == 0)
        {
            return null;
        }

        String[] result= new String[s.length];
        for (int i = 0; i < s.length; i++)
            result[i] = Character.toString(s[i]);

        if (l == 1)
        {
            return result;
        }
        else
        {
            for (int i = 1; i < l; i++)
                result = combine(result, s);

            return result;
        }
    }
    
    
    
    // Just for debugging
    public static void main(String[] args)
    {
        String[] a = {"pep", "set", "ara", "xav", "bot"};
        int[] idxs = {2, 4, 1};
        ArrayFunctions myFunctions = new ArrayFunctions();
        String[] asub = myFunctions.getSubarray(a, idxs);
        
        for (int i = 0; i < asub.length; i++)
            System.out.println(asub[i]);
        
        int[] myInts = ArrayFunctions.getIntegersArray(4, 10);
        for (int i = 0; i < myInts.length; i++)
            System.out.println(myInts[i]);
        
        Double[] b = {12.4, 100.0, 34.7, 0.98, 72.0};
        Integer[] idxs2 = ArrayFunctions.getSortedIndexes(b);
        for (int i = 0; i < idxs2.length; i++)
            System.out.println("old: " + b[i] + "  new: " + b[idxs2[i]] + "  old-idx: " + idxs2[i]);
        
        String[] c = {"pep", "set", "ara", "xav", "bot", "ara", "sepa", "copa"};
        String smallest = myFunctions.getSmallestString(c);
        System.out.println(smallest);
        int[] indices = myFunctions.getIndicesOf(smallest, c);
        for (int i = 0; i < indices.length; i++)
            System.out.println(indices[i]);
        System.out.println(a.equals(c));
        
        int[] idxs3 = myFunctions.getIndicesOf(3, new Integer[]{23,3,43,3,53});
        for (int i = 0; i < idxs3.length; i++)
            System.out.println(idxs3[i]);
    }
}


