package molecules;
/*
*
*
*/

import java.io.*;


public class Enzyme
{
    private String name;
    private String regex;
    private String restrictedAA = "";
    
    public Enzyme(int code)
    {
        switch (code)
        {
            case 1:
                this.name = "Trypsin";
                this.regex = "[^RK]*(K|R)";
                this.restrictedAA = "P";
                break;
            
            case 2:
                this.name = "Trypsin/P";
                this.regex = "[^RK]*(K|R)";
                break;
            
            case 3:
                this.name = "LysC";
                this.regex = "[^K]*K";
                this.restrictedAA = "P";
                break;
            
            case 4:
                this.name = "LysC/P";
                this.regex = "[^K]*K";
                break;
            
            case 5:
                this.name = "ArgC";
                this.regex = "[^R]*R";
                this.restrictedAA = "P";
                break;
            
            case 6:
                this.name = "ArgC/P";
                this.regex = "[^R]*R";
                break;
            
            case 7:
                this.name = "GluC";
                this.regex = "[^E]*E";
                break;
            
            case 8:
                this.name = "Chymotrypsin";
                this.regex = "[^FLWY]*(F|L|W|Y)";
                this.restrictedAA = "P";
                break;
            
            case 9:
                this.name = "PepsinA";
                this.regex = "[^FL]*(F|L)";
                break;
            
            default: // 'TrypsinP'
                this.name = "Trypsin/P";
                this.regex = "[^RK]*(K|R)";
                break;
        }
    }
    
    public Enzyme(String name, String regex, String restrictedAA)
    {
        this.name = name;
        this.regex = regex;
        this.restrictedAA = restrictedAA;
    }
    
    public String getRegex()
    {
        return regex;
    }
    
    public String getRestrictedAA()
    {
        return restrictedAA;
    }
    
    public String getName()
    {
        return name;
    }
}
