/*
*
*
*/

public class MolMeasures
{
    public static double protonMass = 1.0072764666;
	public static double electronMass = 5.4857990943*Math.Pow(10, -4);

	public static readonly double carbonMass = 12.0;
	public static readonly double hydrogenMass = 1.0078250321;
	public static readonly double nitrogenMass = 14.0030740052;
	public static readonly double oxygenMass = 15.9949146221;
	public static readonly double phosphorusMass = 30.97376151;
	public static readonly double sulphurMass = 31.97207069;
	
	public static readonly double cTerminusMass = oxygenMass + hydrogenMass;
	public static readonly double nTerminusMass = hydrogenMass;

	public static readonly double carbonWeight = calcWeight("C");
	public static readonly double hydrogenWeight = calcWeight("H");
	public static readonly double nitrogenWeight = calcWeight("N");
	public static readonly double oxygenWeight = calcWeight("O");
	public static readonly double phosphorusWeight = calcWeight("P");
	public static readonly double sulphurWeight = calcWeight("S");
	public static readonly double cTerminusWeight = oxygenWeight + hydrogenWeight;
	public static readonly double NTerminusWeight = hydrogenWeight;

	public static readonly double BIonMassOffset = -hydrogenMass;
	public static readonly double AIonMassOffset = - carbonMass - oxygenMass - hydrogenMass;
	public static readonly double CIonMassOffset = MassN + 2*hydrogenMass;
	public static readonly double YIonMassOffset = +hydrogenMass;
	public static readonly double XIonMassOffset = carbonMass + oxygenMass - hydrogenMass;
	public static readonly double ZIonMassOffset = - MassN - 2*hydrogenMass;

	public static readonly double MassWater = CalcMonoMass("H2O");
	public static readonly double MassAmmonia = CalcMonoMass("NH3");

}
