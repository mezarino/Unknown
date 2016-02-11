package molecules;
/*
*
*
*/

import java.io.*;
import java.util.Map;
import java.util.HashMap;

import molecules.MassInfo;

public class Modification
{
    private String title;
    private String description;
    private double deltaMass;
    private String position;
    private Character[] sites;
    private String symbol;
    private boolean fixed;
    
    public Modification(int code, boolean fixed) throws Exception
    {
        switch (code)
        {
            case 1:
                this.title = "Acetyl (K)";
                this.description = "Acetylation";
                this.symbol = "[Ac]";
                this.deltaMass = 42.0105646863;
                this.position = "any";
                this.sites = new Character[]{'K'};
                this.fixed = fixed;
                break;
            
            case 2:
                this.title = "Acetyl (Protein N-term)";
                this.description = "Acetylation";
                this.symbol = "[Ac]";
                this.deltaMass = 42.0105646863;
                this.position = "protNterm";
                this.sites = new Character[]{'-'};
                this.fixed = fixed;
                break;
            
            case 3:
                this.title = "Carbamidomethyl (C)";
                this.description = "Iodoacetamide derivative";
                this.symbol = "[Cam]";
                this.deltaMass = 57.0214637236;
                this.position = "any";
                this.sites = new Character[]{'C'};
                this.fixed = fixed;
                break;
            
            case 4:
                this.title = "Oxidation (M)";
                this.description = "Oxidation";
                this.symbol = "[Oxi]";
                this.deltaMass = 15.9949146221;
                this.position = "any";
                this.sites = new Character[]{'M'};
                this.fixed = fixed;
                break;
            
            case 5:
                this.title = "Phospho (STY)";
                this.description = "Phosphorylation";
                this.symbol = "[Phos]";
                this.deltaMass = 79.9663304084;
                this.position = "any";
                this.sites = new Character[]{'S', 'T', 'Y'};
                this.fixed = fixed;
                break;
            
            case 6:
                this.title = "GlyGly (K)";
                this.description = "Ubiquitination residue";
                this.symbol = "[Ubi]";
                this.deltaMass = 114.0429274472;
                this.position = "any";
                this.sites = new Character[]{'K'};
                this.fixed = fixed;
                break;
            
            case 7:
                this.title = "Methyl (KR)";
                this.description = "Methylation";
                this.symbol = "[Me]";
                this.deltaMass = 14.015650064199999;
                this.position = "any";
                this.sites = new Character[]{'K', 'R'};
                this.fixed = fixed;
                break;
            
            case 8:
                this.title = "Dimethyl (KR)";
                this.description = "di-Methylation";
                this.symbol = "[Me2]";
                this.deltaMass = 28.031300128399998;
                this.position = "any";
                this.sites = new Character[]{'K', 'R'};
                this.fixed = fixed;
                break;
            
            case 9:
                this.title = "Trimethyl (KR)";
                this.description = "tri-Methylation";
                this.symbol = "[Me3]";
                this.deltaMass = 42.0469501926;
                this.position = "any";
                this.sites = new Character[]{'K', 'R'};
                this.fixed = fixed;
                break;
            
            case 10:
                this.title = "Deamidation (N)";
                this.description = "Asparagine deamidation";
                this.symbol = "[Deam]";
                this.deltaMass = 0.9840155848;
                this.position = "any";
                this.sites = new Character[]{'N'};
                this.fixed = fixed;
                break;
            
            case 11:
                this.title = "Deamidation 18O (N)";
                this.description = "Asparagine heavy-deamidation";
                this.symbol = "[18ODeam]";
                this.deltaMass = 2.9882613627;
                this.position = "any";
                this.sites = new Character[]{'N'};
                this.fixed = fixed;
                break;
            
            case 12:
                this.title = "18O";
                this.description = "C-terminal heavy water label";
                this.symbol = "[18O]";
                this.deltaMass = 4.0084915558;
                this.position = "pepCterm";
                this.sites = new Character[]{'c'};
                this.fixed = fixed;
                break;
            
            case 13:
                this.title = "Acetyl (Peptide N-term)";
                this.description = "Acetylation";
                this.symbol = "[Ac]";
                this.deltaMass = 42.0105646863;
                this.position = "pepNterm";
                this.sites = new Character[]{'n'};
                this.fixed = fixed;
                break;
            
            case 14:
                this.title = "Amidation (Protein C-term)";
                this.description = "C-terminal amide CONH2";
                this.symbol = "[Am]";
                this.deltaMass = -0.98402;
                this.position = "protCterm";
                this.sites = new Character[]{'_'};
                this.fixed = fixed;
                break;
            
            case 15:
                this.title = "Methylamine (Protein C-term)";
                this.description = "C-terminal methylamine CONHCH3";
                this.symbol = "[Am]";
                this.deltaMass = 13.0316300642;
                this.position = "protCterm";
                this.sites = new Character[]{'_'};
                this.fixed = fixed;
                break;
            
            default:
                throw new Exception("Unknown modification. Use second constructor instead.");
        }
    }
    
    public Modification(String title, String description, String symbol, double deltaMass, String position, Character[] sites, boolean fixed)
    {
        this.title = title;
        this.description = description;
        this.symbol = "[" + symbol + "]";
        this.deltaMass = deltaMass;
        this.position = position;
        this.sites = sites;
        this.fixed = fixed;
    }
    
    public String getPosition()
    {
        return position;
    }
    
    public Character[] getSites()
    {
        return sites;
    }
    
    public double getDeltaMass()
    {
        return deltaMass;
    }
    
    public String getSymbol()
    {
        return symbol;
    }
    
    public boolean isFixed()
    {
        return fixed;
    }
    
    public static Map<Character, Double> arrayToDict(Modification[] modsArray) // Note that if any 2 or more modifications occur in the same kind of aa, only one will be considered. 
    {
        Map<Character, Double> result = new HashMap<Character, Double>();
        //result.put('c', 0.0); // Cterm (not confuse with 'C' cysteine)
        //result.put('n', 0.0); // Nterm (not confuse with 'N' arginine)
        
        for (int i = 0; i < modsArray.length; i++)
        {
            if (modsArray[i] == null)
                continue;
            
            Character[] sites = modsArray[i].getSites();
            for (int j = 0; j < sites.length; j++)
                result.put(sites[j], modsArray[i].getDeltaMass());
        }
        
        return result;
    }
}
