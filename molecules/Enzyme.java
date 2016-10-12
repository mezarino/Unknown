package molecules;
/*
*
*
*/

import java.io.*;


public class Enzyme
{
    private String name;
    private String cleavageSitePattern;
    private String cutSide;
    
    public Enzyme(int code)
    {
        switch (code)
        {
            case 1:
                this.name = "Trypsin";
                this.cleavageSitePattern = "(K|R)[^P]";
                this.cutSide = "C-term";
                break;
            
            case 2:
                this.name = "Trypsin/P";
                this.cleavageSitePattern = "(K|R)";
                this.cutSide = "C-term";
                break;
            
            case 3:
                this.name = "LysC";
                this.cleavageSitePattern = "(K)[^P]";
                this.cutSide = "C-term";
                break;
            
            case 4:
                this.name = "LysC/P";
                this.cleavageSitePattern = "K";
                this.cutSide = "C-term";
                break;
            
            case 5:
                this.name = "ArgC";
                this.cleavageSitePattern = "(R)[^P]";
                this.cutSide = "C-term";
                break;
            
            case 6:
                this.name = "ArgC/P";
                this.cleavageSitePattern = "R";
                this.cutSide = "C-term";
                break;
            
            case 7:
                this.name = "GluC";
                this.cleavageSitePattern = "E";
                this.cutSide = "C-term";
                break;
            
            case 8:
                this.name = "Chymotrypsin";
                this.cleavageSitePattern = "(F|L|W|Y)[^P]";
                this.cutSide = "C-term";
                break;
            
            case 9:
                this.name = "PepsinA";
                this.cleavageSitePattern = "(F|L)";
                this.cutSide = "C-term";
                break;
            
            default: // 'TrypsinP'
                this.name = "Trypsin/P";
                this.cleavageSitePattern = "(K|R)";
                this.cutSide = "C-term";
                break;
        }
    }
    
    public Enzyme(String name, String regex, String termini)
    {
        this.name = name;
        this.cleavageSitePattern = regex;
        this.cutSide = termini;
    }
    
    public String getCleavageSitePattern()
    {
        return cleavageSitePattern;
    }

    public String getCutSide()
    {
        return cutSide;
    }
    
    public String getName()
    {
        return name;
    }
}
