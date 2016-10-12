import java.io.*;

public class FindPrecursorIsotopicPairs {
    public static void main(String[] args) {
        // Obtain the number of ms2 scans in mzXML file
        //System.out.println("Calculating number of MS2 scans...");
        int ms2scancounts = 0;
        
        try
		{
            String mzXMLPath = args[0];
            BufferedReader reader = new BufferedReader(new FileReader(mzXMLPath));
            String st = reader.readLine().trim();
            
            while (true)
            {   
                if (st.startsWith("<precursorMz")) {
                    ms2scancounts++;
                } else if (st.startsWith("<index ")) { // reached end of scans branch in mzXML file
                    break;
                }
                
                st = reader.readLine().trim();
            }
        }
        catch (Exception e)
		{ e.toString(); }
        
        //System.out.println("There are " + ms2scancounts + " MS2 scans.");
        
        // Parse precursor peaks into an Array by decreasing rt
        //System.out.println("Parsing mzXML file and extracting precursor peaks...");
        Peak[] precpeaks = new Peak[ms2scancounts];
        double[] precrts = new double[ms2scancounts];
        double[] precscannumbers = new double[ms2scancounts];
        
        try
		{
            String mzXMLPath = args[0];
            BufferedReader reader = new BufferedReader(new FileReader(mzXMLPath));
            
            while (ms2scancounts > 0)
            {
                String st = reader.readLine().trim();
                String[] tempst;

		        int parentScanNumber = -1, scanNumber = -1, precursorCharge = -1;
		        double retentionTime = -1, precursorIntensity = -1, precursorMz = -1;

		        while (!st.startsWith("<precursorMz"))
		        {
                    if (st.startsWith("<index "))  // reached end of scans branch in mzXML file
                    {
                        break;
                    }
                    else if (st.startsWith("<scan "))
                    {
                        st = st.split("=")[1]; //<scan num="1"
                        st = st.substring(1, st.length() - 1); // remove "..."
                        scanNumber = Integer.valueOf(st);
                    }
                    else if (st.startsWith("retentionTime"))
			        {
				        st = st.split("=")[1];
				        st = st.substring(3, st.length() - 2); // remove "PT...S"
				        retentionTime = Double.valueOf(st) / 60.0; // convert to minutes
			        }

			        st = reader.readLine().trim();
		        }

		        // reached precursorMz
                precursorMz = Double.valueOf(st.substring(st.indexOf(">") + 1, st.length() - 14));
		        tempst = st.split(" ");

		        for (int i = 0; i < tempst.length; i++)
		        {
			        if (tempst[i].startsWith("precursorScanNum"))
			        {
				        st = tempst[i].split("=")[1];
				        st = st.substring(1, st.length() - 1); // remove "..."
				        parentScanNumber = Integer.valueOf(st);
			        }
                    else if (tempst[i].startsWith("precursorIntensity"))
			        {
				        st = tempst[i].split("=")[1];
				        st = st.substring(1, st.length() - 1); // remove "..."
				        precursorIntensity = Double.valueOf(st);
			        }
                    else if (tempst[i].startsWith("precursorCharge"))
			        {
				        st = tempst[i].split("=")[1];
				        st = st.substring(1, st.length() - 1); // remove "..."
				        precursorCharge = Integer.valueOf(st);
			        }
		        }
                
                // Push to peaks array
                precpeaks[ms2scancounts-1] = new Peak(precursorMz, precursorIntensity, precursorCharge);
                precrts[ms2scancounts-1] = retentionTime;
                precscannumbers[ms2scancounts-1] = scanNumber;
                ms2scancounts--;
            }
        }
		catch (Exception e)
		{ e.toString(); }
        
        // Find the isotopic pairs;
        //System.out.println("Finding precursor isotopic pairs...");
        double mztolerance = Double.valueOf(args[1]);
        String mztoleranceunit = args[2];
        double expectedmassshift = Double.valueOf(args[3]);
        double rttolerance = Double.valueOf(args[4]); // in minutes
        int mincharge = Integer.valueOf(args[5]);
        int maxcharge = Integer.valueOf(args[6]);
        
        System.out.println("Scan1\tScan2\tShift\tCharge\t\tMz1\tMz2\tIntensity1\tIntensity2\tRt1\tRt2");
        ms2scancounts = precpeaks.length;
        
        for (int i = 0; i < ms2scancounts - 1; i++) {
            double mzerror = mztolerance;
            if (mztoleranceunit.equals("ppm"))
                mzerror = mztolerance * 1e-6 * precpeaks[i].getMz();
            
            for (int j = i+1; j < ms2scancounts; j++) {
                if (precrts[i] - precrts[j] > rttolerance)
                    break;
                else if (precpeaks[i].getCharge() != precpeaks[j].getCharge())
                    continue;
                
                double mzshift = Math.abs(precpeaks[j].getMz() - precpeaks[i].getMz());
                //if ( mzshift < (expectedmassshift + mzerror) / maxcharge) {
                 //   continue;
                //}
                
                //double expcharge = expectedmassshift / mzshift;
                //double chargetolerance = mzerror / mzshift;
                
                if ( Math.abs(expectedmassshift/precpeaks[i].getCharge() - mzshift) <=  mzerror) {
                    System.out.println( precscannumbers[i] + "\t" + precscannumbers[j] + "\t" + mzshift + "\t" + precpeaks[i].getCharge() + "\t" + precpeaks[i].getMz() + "\t" + precpeaks[j].getMz() + "\t" + precpeaks[i].getIntensity() + "\t" + precpeaks[j].getIntensity() + "\t" + precrts[i] + "\t" + precrts[j]);
                }
            }
        } 
    }
}
