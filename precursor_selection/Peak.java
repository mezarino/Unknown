public class Peak {
    private double mz;
    private double intensity;
    private short charge;
    
    public Peak(double mz, double inten, int z) {
        this.mz = mz;
        this.intensity = inten;
        this.charge = (short) z;
    }
    
    public double getMz() {
        return mz;
    }
    
    public double getIntensity() {
        return intensity;
    }
    
    public int getCharge() {
        return (int) charge;
    }
}
