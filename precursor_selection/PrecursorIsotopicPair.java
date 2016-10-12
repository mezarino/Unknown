/******************************************************************************
*
*
******************************************************************************/


public class PrecursorIsotopicPair {
    private final Peak precursor1;
    private final Peak precursor2;
    private final double mzshift;
    
    // Constructor
    /* precondition: peaks must have the same charge state
    */
    public PrecursorIsotopicPair(Peak p1, Peak p2, double mzshift) {
        if( p1.getCharge() != p2.getCharge() )
            throw new IllegalArgumentException();
        
        precursor1 = p1;
        precursor2 = p2;
        this.mzshift = mzshift;
    }
    
    public int getCharge() {
        return precursor1.getCharge();
    }
    
    public double getMassShift() {
        return mzshift * precursor1.getCharge();
    }
}
