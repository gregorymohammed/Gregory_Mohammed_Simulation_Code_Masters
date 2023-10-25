package Gregory_Mohammed_Simulation_Code_Masters;

public class greg_const {
    /**
     * =======================================
     * Just some general constants (SI Units).
     * =======================================
     */
    final double epsilonCorrection = 1.0e-1;
    final double radiationConstant = 7.5657e-13; /// km^{-1}s^{-2}K^{-4}
    final double stefanBoltzmann = 5.6704e-8; /// kgs^{-3}K^{-4}
    final double stefanBoltzmanncgs = 5.67051e-5; ///erg cm^{-2} K^{-4} s^{-1}
    final double boltzmann = 1.381e-29; /// km^{2}kgs^{-2}K^{-1}
    final double boltzmanncgs = 1.380658e-16; /// erg K^{-1}
    final double electronMass = 9.1093e-31; /// kg.
    final double electronMasscgs = 9.1094e-28; /// g.
    /**
     * =======================================
     * Just some general geometrized constants.
     * =======================================
     */
    final double radiationConstantGeo = 1.68e-55; /// km^{-1}K^{-4}
    final double stefanBoltzmannGeo = 1.26e-52; /// km^{-3}K^{-4}
    final double boltzmannGeo = 3.07e-72; /// km^{2}K^{-1}
    final double mfe = 55.845; /// Atomic Mass Unit of Iron (Fe).
    final double me = 2.0260e-57; /// Mass of Electron.
    /**
     * =======================================
     * Constants for astrophysical properties.
     * Masses are in kilograms, lengths in
     * metres,times in seconds.
     * =======================================
     */
    final double K = 1.0;
    final double GAMMA = 5.0 / 3.0;
    final double RELGAMMA = 4.0 / 3.0;
    final double G = 6.67259e-11;
    final double Gcgs = 6.67259e-8; /// cm^{3} g^{-1} s^{-2}
    final double c = 3.0e8; // m/s
    final double ccgs = 3.0e10; // cm/s
    final double PI = 3.141592654;
    /**
     * =====================================
     */
    final double solarMass = 1.989e30;
    final double solarRadius = 6.9599e8;
    final double solarCenDensity = 6.0e8;
    final double solarCenPressure = 2.5e16;
    /**
     * =====================================
     * Constants for the various classes.
     * =====================================
     */
    final double MAX_RADIUS = 150e9;
    final double MAX_TIME = 1.0e9;
    final int MIN_RADIUS_POINTS_ALLOWED = 10000;
    final int MAX_SIZE_FILE = 10;
    final int MAX_CHARACTERS_INFILE = 11;
    final int MAX_CHARACTERS_LISTFILE = 18;
    final int MAX_CHARACTERS_BATCHFILE = 64;
    final int CHAR_LINE = 1024;
    final int MAX_CHARACTERS_OUTFILE = 16;
    /*
     * =====================================
     */
}
