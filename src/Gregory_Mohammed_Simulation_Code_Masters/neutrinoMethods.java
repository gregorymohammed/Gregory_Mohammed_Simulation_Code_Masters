package Gregory_Mohammed_Simulation_Code_Masters;

/**
 *
 * @author Gregory Mohammed
 */
public final class neutrinoMethods {

    public double neutrinoE;
    public double epsilonnu;
    public double epsilonnugeo;
    public double temperature;
    public double neutrinoP;
    public double neutrinoF;
    public double luminosity;
    public double kappa;
    public double gravitationalPotential;
    public double qt;
    public double qx;

    public neutrinoMethods(greg_const gc) {
        /**
         * \varepsilon_{\nu} = 15 MeV = 2.403e-18 kgkm^{2}s^{-2}.
         */
        epsilonnu = 2.403e-18;
        epsilonnugeo = 1.78e-52; /// km^{-1}

        /**
         * \varepsilon _{\nu} = kT
         */
        temperature = epsilonnu / gc.boltzmann;
        neutrinoEnergy(gc); 
        neutrinoPressure();
        neutrinoFlux(gc);
        luminosity = 1.78e4; /// This is geometrized to km^{-2}.
        kappa = 1.7e-2; /// This is geometrized.

        gravitationalPotential = 0.0;
        qt = 0.0;
        qx = 0.0;

        //printNeutrinoMethods();
    }

    public void temperature(greg_const gc) {
        temperature = epsilonnu / gc.boltzmann;
    }

    public void neutrinoEnergy(greg_const gc) {
        /**
         * This is geometrized.
         */
        neutrinoE = ((gc.radiationConstantGeo * Math.pow(temperature, 4.0))
                * (7.0 + (7.0 * Math.sqrt(3.0)))) / (20.0 + (8.0 * Math.sqrt(3.0)));
    }

    public void neutrinoPressure() {
        /**
         * This is geometrized.
         */
        neutrinoP = neutrinoE / 3.0;
    }

    public void neutrinoFlux(greg_const gc) {
        /**
         * This is geometrized.
         */
        neutrinoF = ((7.0 * gc.stefanBoltzmannGeo * Math.pow(temperature, 4.0))) / (5.0 + (2.0 * Math.sqrt(3)));
    }

    public void makeGravitationalPotential(double x, double rho0, double massCore) {
        if (x == 0.0 || rho0 == 0.0) {
            System.err.println("ERROR :: neutrinoMethods :: makeGravitationalPotential ::\n"
                    + "x = " + x + " => 1/x = " + 1 / x + "\n"
                    + "rho0 = " + rho0);
            System.exit(1);
        }
        gravitationalPotential = (rho0 * massCore) / Math.pow(x, 2.0);
    }

    public double sourceQt(IO io, double W, double x, double vel, double rho0) {
        /**
         * This is geometrized.
         */
        makeGravitationalPotential(x, rho0, io.massCore);
        
        qt = kappa * (rho0 - gravitationalPotential) * ((W * io.energyTuner * neutrinoE) - (io.fluxTuner * neutrinoF * vel));
        return qt;
    }

    public double sourceQx(IO io, double W, double x, double vel, double rho0) {
        /**
         * This is geometrized.
         */
        makeGravitationalPotential(x, rho0, io.massCore);
        
        qx = kappa * (rho0 - gravitationalPotential) * ((io.energyTuner * (neutrinoE / 3.0) * vel) - (W * io.fluxTuner * neutrinoF));
        return qx;
    }

    public void printNeutrinoMethods() {
        System.err.println("OUTPUT INFO :: neutrinoMethods :: All variable values ::\n"
                + "neutrinoE = " + neutrinoE + " :: neutrinoF = " + neutrinoF + " :: neutrinoP = " + neutrinoP + "\n"
                + "epsilonnu = " + epsilonnu + " :: epsilonnugeo = " + epsilonnugeo + " :: luminosity = " + luminosity + "\n"
                + "kappa = " + kappa + " :: temperature = " + temperature + "\n"
                + "qt = " + qt + " :: qx = " + qx);
    }
}
