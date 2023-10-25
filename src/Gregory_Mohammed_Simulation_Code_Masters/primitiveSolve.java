package Gregory_Mohammed_Simulation_Code_Masters;

public class primitiveSolve {
    /// ============================================================

    public int N_ITER_MAX;
    public double PREC;
    public double V_LIGHT_MIN;
    /// ============================================================
    public double f;
    public double fprime;
    public double pcurrent;
    public double dpprev;
    public double dp;
    /// ============================================================
    public double vv;
    public double gam1;
    public double a;
    /// ============================================================
    public double _v;
    public double _rho0;
    public double _p;
    public double _eps;

    public primitiveSolve() {
        N_ITER_MAX = 10;
        PREC = 1.0E-10;
        V_LIGHT_MIN = 0.99;

        pcurrent = 0.0;
        dpprev = 1.0e2;
        dp = 0.0;

        vv = 0.0;
        gam1 = 0.0;
        a = 0.0;

        _v = 0.0;
        _rho0 = 0.0;
        _p = 0.0;
        _eps = 0.0;
    }

    /**
     * =======================================================================
     * Computes f and fprime for use by the Newton-Raphson solver
     * pcurrent  Current value of estimate for p
     * S  Momentum variable
     * D  Density variable
     * E  Energy variable
     * gamma  equation of state
     * =======================================================================
     * fprime OUTPUT: value of derivative of function for Newton-Raphson
     * =======================================================================
     * 0 all well
     * 1 v>=1
     * =======================================================================
     */
    public int fvalue(eqs eqs, final double pcurrent, final double S, final double D, final double E) {

        final double v = S / (E + pcurrent + D);

        if (v >= 1.0) {
            System.out.println("primitiveSolve :: fvalue :: v >= 1.0: values are: \n" + "S = " + S + "\nD = " + D + "\nE = " + E + "\npcurrent = " + pcurrent + "\nv = " + v + "\n");
        }

        vv = v * v;

        if (vv >= V_LIGHT_MIN * V_LIGHT_MIN) {
            return 1;
        }

        a = 1.0 / Math.sqrt(1.0 - vv);

        final double rho0 = D / a;
        final double eps = (E + pcurrent * (1.0 - a * a) + D * (1.0 - a)) / (D * a);

        final double p = eqs.P(rho0, eps);

        f = p - pcurrent;
        fprime = eqs.gam1 * vv * (E + pcurrent + D * (1.0 - a)) / (E + pcurrent + D) - 1.0;

        return 0;
    }

    /**
     * =======================================================================
     * Solver for primitive variables from conserved variables
     * S Momentum variable
     * D Density variable
     * E Energy variable
     * gamma Equation of state (p = (gamma-1) epsilon rho_0
     * pguess starting guess
     * v Output 3-velocity
     * rho0 Output rest density
     * eps Output specific energy density
     * p   Output pressure
     * =======================================================================
     * 0 ok (convergence)
     * 1 did not converge
     * 2 diverging
     * 3 bad initial data
     * =======================================================================
     */
    public int primitiveSolveFunction(eqs eqs, final double S, final double D, final double E, 
            final double pguess, double v, double rho0, double eps, double p) {
        basic bass = new basic();

        if (D <= 0.0) {
            v = 0.0;
            rho0 = 0.0;
            eps = 0.0;
            p = 0.0;

            _v = v;
            _rho0 = rho0;
            _p = p;
            _eps = eps;

            return 0;
        }
        if (E <= 0.0) {
            v = S / D;
            if (Math.abs(v) > V_LIGHT_MIN) {
                v = bass.signum(v) * V_LIGHT_MIN;
            }
            a = 1.0 / Math.sqrt(1.0 - v * v);
            rho0 = D / a;
            eps = 0.0;
            p = 0.0;

            _v = v;
            _rho0 = rho0;
            _p = p;
            _eps = eps;

            return 0;
        }
        if (pguess < 0.0) {
            System.err.println("primitiveSolve :: ERROR :: pguess = " + pguess + "\n");
            return 3;
        }

        pcurrent = pguess;

        for (int iter = 1; iter <= N_ITER_MAX; ++iter) {
            int ierr = fvalue(eqs, pcurrent, S, D, E);
            if (ierr != 0) {
                System.err.println("primitiveSolve :: ERROR :: function returns error :: probably v > 1 (light)\n");
                return 4;
            }

            dp = -f / fprime;
            pcurrent += dp;

            if (pcurrent <= 0.0) {
                //System.err.println("primitiveSolve :: ERROR :: primitiveSolveFunction :: pcurrent = " + pcurrent + " <= 0.0.\nThis is being repaired ...\n");
                v = S / D;
                if (Math.abs(v) > V_LIGHT_MIN) {
                    v = bass.signum(v) * V_LIGHT_MIN;
                }
                a = 1.0 / Math.sqrt(1.0 - v * v);
                rho0 = D / a;
                eps = 0.0;
                p = 0.0;

                _v = v;
                _rho0 = rho0;
                _p = p;
                _eps = eps;

                return 0;
            }

            if (dp <= PREC * pcurrent) {
                /**
                 * ===============================================================
                 * all well, complete the computation and return
                 * ===============================================================
                 */
                v = S / (E + pcurrent + D);
                if (Math.abs(v) > V_LIGHT_MIN) {
                    v = bass.signum(v) * V_LIGHT_MIN;
                }
                a = 1.0 / Math.sqrt(1.0 - v * v);
                rho0 = D / a;
                eps = (E + pcurrent * (1.0 - a * a) + D * (1.0 - a)) / (D * a);
                p = pcurrent;

                _v = v;
                _rho0 = rho0;
                _p = p;
                _eps = eps;

                return 0;
            }
            if (dp > dpprev) {
                System.err.println("primitiveSolve :: ERROR :: Diverging,...\n" + "f = " + f + "\nfprime = " + fprime + "\ncurrent p value = " + pcurrent + "\ndp(previous) = " + dpprev + "\ndp(current) = " + dp + "\n");
                return 2;
            }

            dpprev = dp;
        }
        System.err.println("primitiveSolve :: ERROR :: Did not converge:\n" + "f = " + f + "\nfprime = " + fprime + "\ncurrent p value = " + pcurrent + "\ndp(previous) = " + dpprev + "\ndp(current) = " + dp + "\n");
        return 1;
    }
}
