package Gregory_Mohammed_Simulation_Code_Masters;

/**
 *
 * @author gregory.mohammed
 */
public class geometrize {

    greg_const gconst;
    double georho0l;
    double georho0r;
    double geopl;
    double geopr;
    double geoxl;
    double geoxr;
    double geoepsl;
    double geoepsr;
    double geovl;
    double geovr;
    double rho0Un;
    double pUn;
    double epsUn;
    double vUn;
    double xUn;
    double aUn;
    double Dun;
    double Eun;
    double Sun;

    public geometrize() {
        gconst = new greg_const();

        georho0l = 0.0;
        georho0r = 0.0;
        geopl = 0.0;
        geopr = 0.0;
        geoxl = 0.0;
        geoxr = 0.0;
        geoepsl = 0.0;
        geoepsr = 0.0;
        geovl = 0.0;
        geovr = 0.0;

        rho0Un = 0.0;
        pUn = 0.0;
        epsUn = 0.0;
        vUn = 0.0;
        xUn = 0.0;

        aUn = 0.0;
        Dun = 0.0;
        Eun = 0.0;
        Sun = 0.0;
    }

    public void makeGeometric(eqs eqs, IO io) {
        /**
         * This method is only used with data which is in cgs units. It is also
         * ONLY used with KKT data.
         */
        georho0l = io.lrho0nu * (gconst.Gcgs / Math.pow(gconst.ccgs, 2.0));
        georho0r = io.rrho0nu * (gconst.Gcgs / Math.pow(gconst.ccgs, 2.0));

        geopl = io.lpnu * (gconst.Gcgs / Math.pow(gconst.ccgs, 4.0));
        geopr = io.rpnu * (gconst.Gcgs / Math.pow(gconst.ccgs, 4.0));

        /**
         * Distance in km.
         * The distance is in km in order to zoom in on the active region 
         * in the solution. If this is cm then the solution looks smoothed
         * and uninteresting except for tantalizing hints that something is
         * happening in the central region.
         */
        geoxl = io.xLeft / 1.0e5;
        geoxr = io.xRight / 1.0e5;

        geoepsl = eqs.EPS(geopl, georho0l);
        geoepsr = eqs.EPS(geopr, georho0r);

        geovl = (io.vLeft / gconst.ccgs);
        geovr = (io.vRight / gconst.ccgs);
    }

    public void makeUnGeometric(int signal, double rho0, double p, double eps, double v, double x) {
        /**
         * signal is an artefact of old coding. The if-else statement is never used.
         * signal == 1 is all that is used when makeUnGeometric is called.
         * The conversion is done in cgs units. x is reduced to km in order to 
         * zoom in on the active region in the solution.
         */
        if (signal == 0) {
            rho0Un = rho0;
            pUn = p;
            epsUn = eps;
            vUn = v;
            xUn = x;
        } else if (signal == 1) {
            rho0Un = rho0 * (Math.pow(gconst.ccgs, 2.0) / gconst.Gcgs);
            pUn = p * (Math.pow(gconst.ccgs, 4.0) / gconst.Gcgs);
            epsUn = eps * (Math.pow(gconst.ccgs, 4.0) / gconst.Gcgs);
            vUn = (v * gconst.ccgs);
            xUn = x * 1.0e5; /// km.
        } else {
            System.err.println("ERROR :: geometrize : makeUnGeometric :\n"
                    + "signal is neither 0 nor 1.");
            System.exit(1);
        }
    }
}
