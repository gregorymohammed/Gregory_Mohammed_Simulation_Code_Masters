package Gregory_Mohammed_Simulation_Code_Masters;

public class timeData {

    static final int LOTS_OF_TIME = 0;
    static final int OUT_OF_TIME = 1;
    static final int OUT_OF_TIME_STEPS = 2;
    public double x;
    public double c;
    public double p1;
    public double p2;
    public int n;
    public int n_max;
    public double coord;
    public double coord_max;
    public double delta;
    public double delta_new, delta_min, delta_max;

    public timeData() {
        x = 0.0;
        c = 0.0;
        p1 = 0.0;
        p2 = 0.0;
        n = 0;
        n_max = 0;
        coord = 0.0;
        coord_max = 0.0;
        delta = 0.0;
        delta_new = 0.0;
        delta_min = 0.0;
        delta_max = 0.0;
    }

    public void Make(IO io) {
        c = io.courant;
        p1 = io.p1;
        n_max = io.totalSteps;

        if (io.dataSwitch == 0) {
            coord_max = io.totalTimeSod;
        } else if (io.dataSwitch == 1) {
            coord_max = io.totalTimeSod;
        } else if (io.dataSwitch == 2) {
            coord_max = io.totalTimeKKT * 3.0e10; /// speed of light, ccgs = 3.0e10 cm/s.
        } else {
            System.err.println("ERROR :: timeData :: Make ::\n"
                    + "Unresolved set of data and related time. \n"
                    + "Check your data and see if there is a time that corresponds to it.\n"
                    + "Make sure that it is either the Sod data or the Kuroda et al. data.");
            System.exit(1);
        }

        /// Test the input

        if (n_max <= 0) {
            System.err.println("TimeData: ERROR: n = " + n + "\n");
            System.exit(1);
        }
        if (coord_max <= 0.0) {
            System.err.println("TimeData: ERROR: coord_max = " + coord_max + "\n");
            System.exit(1);
        }
        if (c <= 0.0) {
            System.err.println("TimeData: ERROR: step_control.c = " + c + "\n");
            System.exit(1);
        }
        if (p1 <= 0.0) {
            System.err.println("TimeData: ERROR: step_control.p1 = " + p1 + "\n");
            System.exit(1);
        }

        /// Initialize basic times

        n = 0;
        coord = 0.0;

        /// initial time step

        delta = 0.01 * coord_max;
        delta_min = 1.0e-12 * coord_max;
        delta_max = 0.5 * coord_max;
    }

    public int testForMax() {
        if (coord >= coord_max) {
            return OUT_OF_TIME;
        }
        if (n >= n_max) {
            return OUT_OF_TIME_STEPS;
        }
        return LOTS_OF_TIME;
    }

    /**
     * ================================================================= Time
     * step choice: (new time step put in dtnew)
     *
     * Time step is controlled by: <UL> <LI> sound speed courant limit <LI>
     * velocity limit <LI> viscosity limit </UL>
     */
    public void makeDelta(eqs eqs, grid grid) {
        /**
         * Maximum value given as multiple of previous dt and as multiple of
         * coordinate "freefall time". Note: tunit is usually some
         * characteristic crossing time.
         */
        if ((p1 <= 0.0)) {
            System.err.println("timeData: ERROR: incorrect parameters; p1 = " + p1 + "\n");
            System.exit(1);
        }

        delta_new = p1 * coord_max;
        /**
         * Courant limit and tests. Go through each Cell.
         */
        for (int ie = 0; ie < grid.clist.size(); ++ie) {
            double vmag = Math.abs(grid.clist.get(ie).cur._v);
            double dx = grid.nlist.get(ie + 1).x - grid.nlist.get(ie).x;
            // Courant limits

            double cs = eqs.CS(grid.clist.get(ie).cur._eps);

            // "Velocity" for Courant limit includes sound speed, fluid velocity

            double courant_v = cs + vmag;

            if (delta_new * courant_v > c * dx) {
                delta_new = c * dx / courant_v;
                if (Double.isInfinite(delta_new)) {
                    System.out.println("timeData :: makeDelta() :: delta_new = " + delta_new);
                    System.exit(1);
                }
            }
            // Artificial Viscosity Limits

            double q = 0.5 * (grid.nlist.get(ie).q + grid.nlist.get(ie + 1).q);
            double artvis_v = Math.sqrt(q / grid.clist.get(ie).cur._rho0);
            double artvis_c = 0.25 * c;

            if (delta_new * artvis_v > artvis_c * dx) {
                delta_new = artvis_c * dx / artvis_v;
            }
        }

        // Some tests (can be removed for production version)

        if (delta_new <= delta_min) {
            System.err.println("timeData: ERROR: small dt. delta_new = " + delta_new + " delta = " + delta + "\n");
            System.exit(1);
        }

        /**
         * Decrease early time steps to get a good start
         */
        final int n_decrease = 5;

        if (n < n_decrease) {
            delta_new *= (double) (n + 1) / (double) (n_decrease + 1);
        }

        /**
         * Maximum time step allowed
         */
        if (delta_new > 1.5 * delta) {
            delta_new = 1.5 * delta;
        }

        /**
         * update the time step
         */
        delta = delta_new;
    }

    public void Update() {
        ++n;
        coord += delta;
    }
}
