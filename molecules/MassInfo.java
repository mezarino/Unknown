package molecules;

import java.util.TreeMap;
import java.util.ArrayList;
import java.util.HashSet;

public class MassInfo
{
	// frequently used mass values
	public static final double neutron = 1.003355; // for heavy isotope mass calculation
	public static final double hydrogen = 1.007825;
	public static final double proton = 1.007276;
	public static final double water = 18.010565;
	public static final double ammonia = 17.026549;
	public static final double aIonLoss = 27.994915;
	public static final double metOxLoss = 63.999978;

	// atom mass table
	static TreeMap<Character, Double> atomMasses = new TreeMap<Character, Double>()
	{{
		put(new Character('C'), new Double(12.0));
		put(new Character('H'), new Double(1.007825));
		put(new Character('O'), new Double(15.994915));
		put(new Character('N'), new Double(14.003074));
		put(new Character('S'), new Double(31.972072));
		put(new Character('P'), new Double(30.973763));
	}};

	// amino acid mass table
	static TreeMap<Character, Double> aminoacidMasses = new TreeMap<Character, Double>()
	{{
		put(new Character('G'), new Double(57.02146));
		put(new Character('A'), new Double(71.03711));
		put(new Character('S'), new Double(87.03203));
		put(new Character('P'), new Double(97.05276));
		put(new Character('V'), new Double(99.06841));
		put(new Character('T'), new Double(101.04768));
		put(new Character('C'), new Double(103.00919));
		put(new Character('L'), new Double(113.08406));
		put(new Character('I'), new Double(113.08406));
		put(new Character('N'), new Double(114.04293));
		put(new Character('D'), new Double(115.02694));
		put(new Character('Q'), new Double(128.05858));
		put(new Character('K'), new Double(128.09496));
		put(new Character('E'), new Double(129.04259));
		put(new Character('M'), new Double(131.04049));
		put(new Character('H'), new Double(137.05891));
		put(new Character('F'), new Double(147.06841));
		put(new Character('R'), new Double(156.10111));
		put(new Character('Y'), new Double(163.06333));
		put(new Character('W'), new Double(186.07931));
		put(new Character('['), new Double(0));
		put(new Character(']'), new Double(0));
	}};
    
    // atom weight table
	static TreeMap<Character, Double> atomWeights = new TreeMap<Character, Double>()
	{{
		put(new Character('C'), new Double(12.0107));
		put(new Character('H'), new Double(1.00794));
		put(new Character('O'), new Double(15.9994));
		put(new Character('N'), new Double(14.0067));
		put(new Character('S'), new Double(32.065));
		put(new Character('P'), new Double(30.973761));
	}};
    
    // amino acid weight table
    static TreeMap<Character, Double> aminoacidWeights = new TreeMap<Character, Double>()
	{{
		put(new Character('G'), new Double(57.0513));
		put(new Character('A'), new Double(71.0779));
		put(new Character('S'), new Double(87.0773));
		put(new Character('P'), new Double(97.1152));
		put(new Character('V'), new Double(99.1311));
		put(new Character('T'), new Double(101.1039));
		put(new Character('C'), new Double(103.1429));
		put(new Character('L'), new Double(113.1576));
		put(new Character('I'), new Double(113.1576));
		put(new Character('N'), new Double(114.1026));
		put(new Character('D'), new Double(115.0874));
		put(new Character('Q'), new Double(128.1292));
		put(new Character('K'), new Double(128.1723));
		put(new Character('E'), new Double(129.114));
		put(new Character('M'), new Double(131.1961));
		put(new Character('H'), new Double(137.1393));
		put(new Character('F'), new Double(147.1793));
		put(new Character('R'), new Double(156.1857));
		put(new Character('Y'), new Double(163.1733));
		put(new Character('W'), new Double(186.2099));
		put(new Character('['), new Double(0));
		put(new Character(']'), new Double(0));
	}};
    
    public static final HashSet<Character> aminoAcids = new HashSet<Character>()
    {{
        add(new Character('G')); add(new Character('A'));
        add(new Character('S')); add(new Character('P'));
        add(new Character('V')); add(new Character('T'));
        add(new Character('C')); add(new Character('L'));
        add(new Character('I')); add(new Character('N'));
        add(new Character('D')); add(new Character('Q'));
        add(new Character('K')); add(new Character('E'));
        add(new Character('M')); add(new Character('H'));
        add(new Character('F')); add(new Character('R'));
        add(new Character('Y')); add(new Character('W'));
    }};

	public MassInfo()
	{}

	// return mass of target molecule
	public static double getMass(String molecule, boolean isAminoAcid)
	{
		double mass = 0;

		if (isAminoAcid) // determine which mass table to use
		{
			for (int i = 0; i < molecule.length(); i++)
			{
                Character aa = new Character(molecule.charAt(i));
                if (aminoacidMasses.containsKey(aa))
				    mass += aminoacidMasses.get(new Character(molecule.charAt(i))).doubleValue();
                else
                    System.err.println("Unknown aa " + aa + ". Skipping");
            }
            
			mass += water;
		}

		else
		{
			for (int i = 0; i < molecule.length(); i++)
				mass += atomMasses.get(new Character(molecule.charAt(i))).doubleValue();
		}

		return mass;
	}
    
    // return weight of target molecule
	public static double getWeight(String molecule, boolean isAminoAcid)
	{
		double weight = 0;

		if (isAminoAcid) // determine which weight table to use
		{
			for (int i = 0; i < molecule.length(); i++)
            {
                Character aa = new Character(molecule.charAt(i));
                if (aminoacidWeights.containsKey(aa))
				    weight += aminoacidWeights.get(new Character(molecule.charAt(i))).doubleValue();
                else
                    System.err.println("Unknown aa " + aa + ". Skipping");
            }
            
			weight += water;
		}

		else
		{
			for (int i = 0; i < molecule.length(); i++)
				weight += atomWeights.get(new Character(molecule.charAt(i))).doubleValue();
		}

		return weight;
	}
    
    public static TreeMap<Character, Integer> getCompositionTemplate()
    {
        // amino acid composition template
        TreeMap<Character, Integer> compositionTemplate = new TreeMap<Character, Integer>()
	    {{
		    put(new Character('G'), new Integer(0));
		    put(new Character('A'), new Integer(0));
		    put(new Character('S'), new Integer(0));
		    put(new Character('P'), new Integer(0));
		    put(new Character('V'), new Integer(0));
		    put(new Character('T'), new Integer(0));
		    put(new Character('C'), new Integer(0));
		    put(new Character('L'), new Integer(0));
		    put(new Character('I'), new Integer(0));
		    put(new Character('N'), new Integer(0));
		    put(new Character('D'), new Integer(0));
		    put(new Character('Q'), new Integer(0));
		    put(new Character('K'), new Integer(0));
		    put(new Character('E'), new Integer(0));
		    put(new Character('M'), new Integer(0));
		    put(new Character('H'), new Integer(0));
		    put(new Character('F'), new Integer(0));
		    put(new Character('R'), new Integer(0));
		    put(new Character('Y'), new Integer(0));
		    put(new Character('W'), new Integer(0));;
	    }};
        
        return compositionTemplate;
    }
}


