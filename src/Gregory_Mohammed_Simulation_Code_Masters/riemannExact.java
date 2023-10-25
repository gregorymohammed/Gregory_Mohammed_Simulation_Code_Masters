/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

import java.util.*;

/**
 *
 * @author gregory
 */
public class riemannExact {

    ///==========================================================
    /**
     * Elapsed time required
     */
    double T;
    /**
     * Position of discontinuity
     */
    double x_discontinuity;
    /**
     * Positions characteristics (left-going or right-going)
     */
    double x1;
    double x2;
    double x3;
    double x4;
    double x5;
    double ps;
    double vels;
    waveVariables lwave;
    waveVariables rwave;
    ///==========================================================

    public void Print() {
        System.out.print(".x_discontinuity=" + x_discontinuity + " .x1=" + x1 + " .x2=" + x2 + " .x3=" + x3 + " .x4=" + x4 + " .x5=" + x5
                + " .vels=" + vels + " .T=" + T);
    }

    public enum LeftRight {

        LEFT, RIGHT
    }

    public double max(double a, double b) {
        return ((a) > (b) ? (a) : (b));
    }

    public double signum(double a) {
        return ((a) >= 0.0 ? (+1) : (-1));
    }

    public double DSIGN(double a, double b) {
        return Math.abs(a) * signum(b);
    }
    ///==========================================================

    public riemannExact() {
        T = -10.0;
        x_discontinuity = 0.0;
        x1 = 0.0;
        x2 = 0.0;
        x3 = 0.0;
        x4 = 0.0;
        x5 = 0.0;
        ps = 0.0;
        vels = -10.0;  //PJM  c=1 so this default value should be right out

        lwave = new waveVariables();
        rwave = new waveVariables();
    }

    /**
     * Solve Calculates the exact solution given left and right states.
     *
     * @param[in] T_in
     * @param[in] x_discontinuity_in
     * @param[in] eqs
     * @param[in] lstate
     * @param[in] rstate
     */
    public void Solve(double T_in, double x_discontinuity_in, eqs eqs, stateVariables lstate, stateVariables rstate) throws Exception {
        if (T_in <= 0.0) {
            System.out.println("riemannExact:Solve: ERROR: input paramater T_in=" + T_in + " <= 0");
            System.exit(-1);
        }
        //System.out.println( "riemannExact:Solve: INFO: eqs.gamma=" + eqs.gamma );
        T = T_in;
        x_discontinuity = x_discontinuity_in;

        /// Solve for min/max pressures over states

        int istraddle = 0;
        final int N_STRADDLE_MAX = 10;

        double pmin = (lstate.p + rstate.p) / 2.0;
        double pmax = pmin;

        double check = 1.0;

        while (check > 0.0) {   /// "5" loop in the fortran
            ++istraddle;
            if (istraddle > N_STRADDLE_MAX) {
                Exception straddleError = new Exception("\nCannot find a straddle.");

                /// Throw an exception for use in the flux calculation higher up.
                throw straddleError;
            }

            pmin = 0.5 * max(pmin, 0.0);
            pmax = 2.0 * pmax;

            double dvel1 = GetDVel(eqs, pmin, lstate, rstate, lwave, rwave);
            double dvel2 = GetDVel(eqs, pmax, lstate, rstate, lwave, rwave);

            //System.out.println("dvel1 = " + dvel1 + "\ndvel2 = " + dvel2 + "\n");

            check = dvel1 * dvel2;
        }
        /**
         * --------------------------- PRESSURE AND FLOW VELOCITY IN THE
         * INTERMEDIATE STATES ---------------------------
         */
        ps = GetP(eqs, pmin, pmax, lstate, rstate, lwave, rwave);
        vels = 0.5 * (lwave.vel + rwave.vel);

        /**
         * ----------- POSITIONS OF THE WAVES -----------
         */
        if (lstate.p >= ps) {
            x1 = x_discontinuity
                    + T * (lstate.vel - lstate.cs) / (1.0 - lstate.vel * lstate.cs);
            x2 = x_discontinuity
                    + T * (vels - lwave.cs) / (1.0 - vels * lwave.cs);
        } else {
            x1 = x_discontinuity + lwave.vshock * T;
            x2 = x1;
        }

        x3 = x_discontinuity + vels * T;

        if (rstate.p >= ps) {
            x4 = x_discontinuity
                    + T * (vels + rwave.cs) / (1.0 + vels * rwave.cs);
            x5 = x_discontinuity
                    + T * (rstate.vel + rstate.cs) / (1.0 + rstate.vel * rstate.cs);
        } else {
            x4 = x_discontinuity + rwave.vshock * T;
            x5 = x4;
        }
    }
///=========================================================

    /**
     * ------- NAME: G E T P -------
     *
     * PURPOSE: FIND THE PRESSURE IN THE INTERMEDIATE STATE OF A RIEMANN PROBLEM
     * IN RELATIVISTIC HYDRODYNAMICS
     *
     *
     * COMMENTS: THIS ROUTINE USES A COMBINATION OF INTERVAL BISECTION AND
     * INVERSE QUADRATIC INTERPOLATION TO FIND THE ROOT IN A SPECIFIED INTERVAL.
     * IT IS ASSUMED THAT DVEL(PMIN) AND DVEL(PMAX) HAVE OPPOSITE SIGNS WITHOUT
     * A CHECK. ADAPTED FROM "COMPUTER METHODS FOR MATHEMATICAL COMPUTATION", BY
     * G. E. FORSYTHE, M. A. MALCOLM, AND C. B. MOLER, PRENTICE-HALL, ENGLEWOOD
     * CLIFFS N.J.
     *
     * /retval pressure (PS)
     */
    public double GetP(eqs eqs, double pmin, double pmax, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {

        double EPS = 1.0;
        double tmp1 = 2.0;

        while (tmp1 > 1.0) {
            EPS = EPS / 2.0;
            tmp1 = 1.0 + EPS;
        }

        final double TOL = 10.0 * EPS;
        final int N_MAX_ITERATIONS = 100;

        //printAnything.message("RiemannExact::GetP: INFO: Entry: machine precision: EPS  = " + EPS + "\n");
        /**
         * INITIALIZATION -looks like a,b and function values f(a) and f(b) are
         * initialized.
         */
        double a = pmin;
        double b = pmax;

        double fa = GetDVel(eqs, a, lstate, rstate, lwave, rwave);
        double fb = GetDVel(eqs, b, lstate, rstate, lwave, rwave);

        /// Check to see if there is a straddle

        if (fa * fb > 0.0) {
            System.err.println("RiemannExact::GetP: ERROR: fa and fb have same sign\n");
            System.exit(1);
        }

        //printAnything.message("RiemannExact::GetP: INFO: on entry\n a=" + a + " fa=" + fa + "\n" + "  b=" + b + " fb=" + fb + "\n");

        /// Initialize the "c" point (for inverse quadratic interpolation)

        double c = a;
        double fc = fa;
        double diff = b - a;
        double old_diff = diff;

        /// Iterate

        for (int i = 1; i <= N_MAX_ITERATIONS; ++i) {

            /// Rearrange so |f(b)| <= |f(c)|
            /// I don't really understand this.  "a" is used as a temp variable??

            if (Math.abs(fc) < Math.abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            /// CONVERGENCE TEST
            /// I don't know why there are factors of 2.0 and 0.5 here?

            double RELATIVE_TOL_PRESSURE = 2.0 * EPS * Math.abs(b) + 0.5 * TOL;
            double xm = 0.5 * (c - b);

            if (Math.abs(xm) <= RELATIVE_TOL_PRESSURE) {
                return b;
            }

            if (fb == 0.0) {
                return b;
            }

            /// IS BISECTION NECESSARY?

            if ((Math.abs(old_diff) < RELATIVE_TOL_PRESSURE) || (Math.abs(fa) <= Math.abs(fb))) {
                diff = xm;
                old_diff = diff;
            } else {

                /// IS QUADRATIC INTERPOLATION POSSIBLE?

                double p, q;
                if (a == c) {

                    /// LINEAR INTERPOLATION

                    double fbfa = fb / fa;
                    p = 2.0 * xm * fbfa;
                    q = 1.0 - fbfa;
                } else {

                    /// INVERSE QUADRATIC INTERPOLATION

                    if (fc == 0) {
                        return c;
                    }
                    if (fa == 0) {
                        return a;
                    }

                    double fafc = fa / fc;
                    double fbfc = fb / fc;
                    double fbfa = fb / fa;

                    p = fbfa * (2.0 * xm * fafc * (fafc - fbfc) - (b - a) * (fbfc - 1.0));
                    q = (fafc - 1.0) * (fbfc - 1.0) * (fbfa - 1.0);
                }

                /// ADJUST SIGNS

                if (p > 0.0) {
                    q = -q;
                }
                p = Math.abs(p);

                /// IS INTERPOLATION ACCEPTABLE?

                if (((2.0 * p) >= (3.0 * xm * q - Math.abs(RELATIVE_TOL_PRESSURE * q)))
                        || (p >= Math.abs(0.5 * old_diff * q))) { //no, back to bisection
                    /// Setup for bisection
                    diff = xm;
                    old_diff = diff;
                } else {
                    /// Setup for using the interpolation
                    old_diff = diff;
                    diff = p / q;
                }
            }

            /// COMPLETE THE STEP
            /// a is replaced by b, and b is given the new (improved) value

            a = b;
            fa = fb;

            if (Math.abs(diff) > RELATIVE_TOL_PRESSURE) { /// Just ensure DIFF>roundoff
                b = b + diff;
            } else {                                  /// add a pre-defined small change
                b = b + DSIGN(RELATIVE_TOL_PRESSURE, xm);
            }
            fb = GetDVel(eqs, b, lstate, rstate, lwave, rwave);

            if (fb * fc > 0.0) {
                c = a;
                fc = fa;
                diff = b - a;
                old_diff = diff;
            }
        } /// Iteration loop

        System.err.println("RiemannExact::GetP: ERROR: Did not converge after N_MAX_ITERATIONS = " + N_MAX_ITERATIONS + "\n");
        System.err.println("a = " + a + " fa = " + fa + "\n"
                + "b = " + b + " fb = " + fb + "\n"
                + "diff = " + diff + " old_diff = " + old_diff + "\n");
        System.err.println("lstate:");
        lstate.Print();
        System.err.println();
        System.err.println("rstate:");
        rstate.Print();
        System.err.println();
        System.exit(1);
        return 1;
    }
//=======================================================

    /**
     * GetDVel
     *
     * NAME: G E T D V E L ----------
     *
     * PURPOSE: COMPUTE THE DIFFERENCE IN FLOW SPEED BETWEEN LEFT AND RIGHT
     * INTERMEDIATE STATES FOR GIVEN LEFT AND RIGHT STATES AND PRESSURE
     *
     * COMMENTS NONE
     *
     * \retval difference (dvel)
     */
    public double GetDVel(eqs eqs, double p, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        /**
         * ----- Waves -----
         */
        LeftRight LEFT = LeftRight.LEFT;
        LeftRight RIGHT = LeftRight.RIGHT;

        int ierr_left = GetVel(eqs.gamma, lstate, lwave, p, LEFT);
        if (ierr_left != 0) {
            System.out.println();
            System.out.println("GetDVel: ERROR: GetVel (left side) returns " + ierr_left + " at p=" + p);
            System.out.print("lstate: ");
            lstate.Print();
            System.out.println();
            System.out.print("lwave: ");
            lwave.Print();
            System.out.println();
            System.exit(1);
        }
        int ierr_right = GetVel(eqs.gamma, rstate, rwave, p, RIGHT);
        if (ierr_right != 0) {
            System.out.println();
            System.out.println("GetDVel: ERROR: GetVel (right side) returns " + ierr_right + " at p=" + p);
            System.out.print("rstate: ");
            rstate.Print();
            System.out.println();
            System.out.print("rwave: ");
            rwave.Print();
            System.out.println();
            System.exit(1);
        }

        double dvel = lwave.vel - rwave.vel;
        return dvel;
    }
///====================================================================

    /**
     * Find solution at a point
     *
     * This is a general approach for any point. Might be faster if there is a
     * priori knowledge of the region the point is in.
     *
     * \param[in] x Position that values are to be calculated. \param[out] node
     * the values on this node
     */
    public void GetExactValues(final double x, eqs eos, exactNode node,
            stateVariables lstate, stateVariables rstate) {
        LeftRight LEFT = LeftRight.LEFT;
        LeftRight RIGHT = LeftRight.RIGHT;

        if (x <= x1) {
            node.p = lstate.p;
            node.rho0 = lstate.rho0;
            node.vel = lstate.vel;
            node.u = lstate.u;

        } else if (x <= x2) {

            double xi = (x - x_discontinuity) / T;
            RareF(xi, eos.gamma, lstate, node, LEFT);
        } else if (x <= x3) {

            node.p = ps;
            node.rho0 = lwave.rho0;
            node.vel = vels;
            node.u = lwave.u;
        } else if (x <= x4) {

            node.p = ps;
            node.rho0 = rwave.rho0;
            node.vel = vels;
            node.u = rwave.u;
        } else if (x <= x5) {

            double xi = (x - x_discontinuity) / T;
            RareF(xi, eos.gamma, rstate, node, RIGHT);
        } else {

            node.p = rstate.p;
            node.rho0 = rstate.rho0;
            node.vel = rstate.vel;
            node.u = rstate.u;
        }
        node.MakewhDSE();
    }

//====================================================================
    /**
     * --------- NAME: G E T V E L ---------
     *
     * PURPOSE: COMPUTE THE FLOW VELOCITY BEHIND A RAREFACTION OR SHOCK IN TERMS
     * OF THE POST-WAVE PRESSURE FOR A GIVEN STATE AHEAD THE WAVE IN A
     * RELATIVISTIC FLOW
     *
     *
     * COMMENTS: THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND
     * MUELLER, J. FLUID MECH., (1994)
     *
     *
     * In original code the first set of variables (with A postfix) are State
     * variables. The second set (no postfix) are wave vars.
     */
    public int GetVel(final double gamma, stateVariables state, waveVariables wave, final double p, final LeftRight LR) {
        final double gam1 = gamma - 1.0;

        //PJM: Test to see if state.p ~= p.  This function fails in this situation (0/0 issues)

        final double EPS = 1.0e-13;  // Just a quick estimate of machine eps

        final double p_diff = p - state.p;
        if (Math.abs(p_diff) < EPS) {
            wave.h = state.h;
            wave.u = state.u;
            wave.rho0 = state.rho0;
            wave.cs = state.cs;
            wave.vshock = 0.0;
            wave.vel = state.vel;
            return 0;
        }

        /**
         * --------------- LEFT OR RIGHT PROPAGATING WAVE ---------------
         */
        LeftRight LEFT = LeftRight.LEFT;
        LeftRight RIGHT = LeftRight.RIGHT;

        double SIGN = 0.0;

        if (LR == LEFT) {
            SIGN = -1.0;
        }
        if (LR == RIGHT) {
            SIGN = 1.0;
        }

        if (p > state.p) {
            /**
             * --- SHOCK ---
             */
            double a = 1.0 + gam1 * (state.p - p) / (gamma * p);
            double b = 1.0 - a;
            double c = state.h * (state.p - p) / state.rho0 - state.h * state.h;
            /**
             * ---------------- CHECK FOR UNPHYSICAL ENTHALPIES ----------------
             */
            //if (c * 4.0 * a > b * b) {
            if (c > (b * b) / (4.0 * a)) {
                System.err.println("RiemannExact::GetVel: ERROR: in shock section: unphysical specific enthalpy in intermediate state\n");
                System.out.println("wave:");
                wave.Print();
                System.out.println();
                return -1;
                //System.exit(1);
            }
            /**
             * ----------------------------- SPECIFIC ENTHALPY IN THE POST-WAVE
             * STATE (FROM THE EQUATION OF STATE AND THE TAUB ADIABAT, EQ.(74),
             * MM94) -----------------------------
             */
            wave.h = (-b + Math.sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

            if (wave.h <= 1.0) {
                System.err.println("RiemannExact::GetVel: ERROR: in shock section: wave.h = " + wave.h + " (<1.0)\nunphysical specific enthalpy in intermediate state\n");
                System.out.println("state:");
                state.Print();
                System.out.println();
                System.out.println("wave:");
                wave.Print();
                System.out.println();
                return -1;
                //System.exit(1);
            }
            /**
             * --------------- DENSITY IN THE POST-WAVE STATE (FROM EQ.(73),
             * MM94) ---------------
             */
            wave.rho0 = gamma * p / (gam1 * (wave.h - 1.0));
            /**
             * ------------------------ SPECIFIC INTERNAL ENERGY IN THE
             * POST-WAVE STATE (FROM THE EQUATION OF STATE)
             * ------------------------
             */
            wave.u = p / (gam1 * wave.rho0);
            /**
             * -------------------------- MASS FLUX ACROSS THE WAVE (FROM THE
             * RANKINE-HUGONIOT RELATIONS, EQ.(71), MM94)
             * --------------------------
             */
            final double denom = state.h / state.rho0 - wave.h / wave.rho0;
            if (Math.abs(denom) < 1.0e-20) {
                System.err.println("RiemannExact::GetVel: ERROR: in shock section denominator=0\nUsually due to state ~= wave so 0/0 issues.\n");
                System.out.println("   denom=" + denom + "  p=" + p + "  state.p=" + state.p);
                System.out.println("state:");
                state.Print();
                System.out.println();
                System.out.println("wave:");
                wave.Print();
                System.out.println();
                return -1;
                //System.exit(1);
            }
            double mflux_sqrt = (p - state.p) / denom;
            if (mflux_sqrt < 0.0) {
                System.err.println("RiemannExact::GetVel: ERROR: in shock section sqrt argument < 0");
                System.out.println("  denom=" + denom + "  p=" + p + " state.p=" + state.p);
                System.out.println("state:");
                state.Print();
                System.out.println();
                System.out.println("wave:");
                wave.Print();
                System.out.println();
                return -1;
                //System.exit(1);
            }
            double mflux = SIGN * Math.sqrt(mflux_sqrt);
            /**
             * ---------- SHOCK VELOCITY (FROM EQ.(86), MM94 ----------
             */
            double mflux2 = mflux * mflux;

            double tmpsa = mflux2 + Math.pow((state.rho0 * state.w), 2.0);
            double tmpsb = -state.vel * state.rho0 * state.rho0 * state.w * state.w;

            wave.vshock = (-tmpsb + SIGN * mflux2 * Math.sqrt(1.0 + state.rho0 * state.rho0 / mflux2))
                    / tmpsa;
            double wshock = 1.0 / Math.sqrt(1.0 - wave.vshock * wave.vshock);
            /**
             * ------------------- FLOW VELOCITY IN THE POST-SHOCK STATE (FROM
             * EQ.(67), MM94) -------------------
             */
            double tmpvela = wshock * (p - state.p) / mflux + state.h * state.w * state.vel;
            double tmpvelb = state.h * state.w
                    + (p - state.p) * (wshock * state.vel / mflux + 1.0 / (state.rho0 * state.w));

            if (tmpvelb == 0.0) {
                System.err.println("RiemannExact::GetVel: ERROR in shock section.  tmpvelb = " + tmpvelb + "\n");
                System.out.println("state:");
                state.Print();
                System.out.println();
                System.out.println("wave:");
                wave.Print();
                System.out.println();
                return -1;
                //System.exit(1);
            }

            wave.vel = tmpvela / tmpvelb;
            /**
             * --------------------- LOCAL SOUND SPEED IN THE POST-SHOCK STATE
             * (FROM THE EQUATION OF STATE) ---------------------
             */
            wave.cs = Math.sqrt(gamma * p / (wave.rho0 * wave.h));

        } else {
            /**
             * ------ RAREFACTION ------
             *
             * --------------------------- POLYTROPIC CONSTANT OF THE GAS ACROSS
             * THE RAREFACTION ---------------------------
             */
            if (state.p <= 0.0) {
                System.err.println("RiemannExact::GetVel: ERROR in rarefaction section: state.p = " + state.p + " <= 0\n");
                System.out.println("state:");
                state.Print();
                System.out.println();
                System.out.println("wave:");
                wave.Print();
                System.out.println();
                return -1;
                //System.exit(1);
            }

            double k = state.p / Math.pow(state.rho0, gamma);
            /**
             * --------------- DENSITY BEHIND THE RAREFACTION ---------------
             */
            wave.rho0 = Math.pow((p / k), (1.0 / gamma));
            //printAnything.message("RiemannExact::GetVel: INFO: Rarefaction section: p = " + p + " state.p = " + state.p + " wave.rho0 = " + wave.rho0 + "\n");
            /**
             * ------------------------ SPECIFIC INTERNAL ENERGY BEHIND THE
             * RAREFACTION (FROM THE EQUATION OF STATE) ------------------------
             */
            wave.u = p / (gam1 * wave.rho0);
            /**
             * -------------------- LOCAL SOUND SPEED BEHIND THE RAREFACTION
             * (FROM THE EQUATION OF STATE) --------------------
             */
            wave.cs = Math.sqrt(gamma * p / (wave.rho0 + gamma * p / gam1));
            /**
             * ------------------ FLOW VELOCITY BEHIND THE RAREFACTION
             * ------------------
             */
            double sqrt_gam1 = Math.sqrt(gam1);

            double tmp3 = ((sqrt_gam1 + state.cs) / (sqrt_gam1 - state.cs))
                    * (sqrt_gam1 - wave.cs) / (sqrt_gam1 + wave.cs);

            double tmp_exp = -SIGN * 2.0 / sqrt_gam1;

            double tmp_a = ((1.0 + state.vel) / (1.0 - state.vel))
                    * Math.pow(tmp3, tmp_exp);

            wave.vel = (tmp_a - 1.0) / (tmp_a + 1.0);
            //printAnything.message("RiemannExact::GetVel: INFO: rarefaction section: wave.vel = " + wave.vel + "\n");
        }
        return 0;
    }
///===================================================================

    /**
     * -------- NAME: R A R E F --------
     *
     * PURPOSE: COMPUTE THE FLOW STATE IN A RAREFACTION FOR GIVEN PRE-WAVE STATE
     *
     *
     * COMMENTS: THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND
     * MUELLER, J. FLUID MECH., (1994)
     *
     * SUBROUTINE RAREF( XI, RHOA, PA, UA, CSA, VELA, S, RHO, P, U, VEL )
     *
     * RHOA, PA, UA, CSA, VELA are input StateVars (COMMON/STATES/)
     * S,RHO,P,U,VEL are separate values in the vectors
     *
     * PJM Note: they get the 'A' ending mixed up!!! Have to be very careful
     * during translation of their code.
     *
     */
    public void RareF(final double xi, final double gamma, stateVariables state, exactNode node, final LeftRight LR) {
        final double gam1 = gamma - 1.0;
        /**
         * --------------- LEFT OR RIGHT PROPAGATING WAVE ---------------
         */
        LeftRight LEFT = LeftRight.LEFT;
        LeftRight RIGHT = LeftRight.RIGHT;

        double SIGN = 0.0;

        if (LR == LEFT) {
            SIGN = 1.0;
        } else if (LR == RIGHT) {
            SIGN = -1.0;
        } else {
            System.err.println("RiemannExact::Raref: ERROR: LR is neither LEFT nor RIGHT\n");
            System.exit(1);
        }

        double b = Math.sqrt(gamma - 1.0);
        double c = (b + state.cs) / (b - state.cs);
        double d = -SIGN * b / 2.0;
        double k = (1.0 + xi) / (1.0 - xi);
        double l = c * Math.pow(k, d);
        double v = Math.pow((1.0 - state.vel) / (1.0 + state.vel), d);

        double ocs2 = state.cs;
        double cs2 = 0.0;

        /**
         * This commented out section is an erroneous adaptation from
         * the original PJM C code.
         */
//        label25:
//        {
//            double fcs2 = l * v * Math.pow((1.0 + SIGN * ocs2), d) * (ocs2 - b)
//                    + Math.pow((1.0 - SIGN * ocs2), d) * (ocs2 + b);
//
//            double dfdcs2 = l * v * Math.pow((1.0 + SIGN * ocs2), d)
//                    * (1.0 + SIGN * d * (ocs2 - b) / (1.0 + SIGN * ocs2))
//                    + Math.pow((1.0 - SIGN * ocs2), d)
//                    * (1.0 - SIGN * d * (ocs2 + b) / (1.0 - SIGN * ocs2));
//
//            cs2 = ocs2 - fcs2 / dfdcs2;
//            if (Math.abs(cs2 - ocs2) / ocs2 > 5.0e-7) {
//                ocs2 = cs2;
//                break label25;
//            }
//        }
        
        /**
         * This is PJM's corrections.
         */
        double diff = 1.0;
        while( diff > 1.0e-7 )             // INSTEAD OF THE GOTO
        {
            double fcs2 = l * v * Math.pow((1.0 + SIGN * ocs2), d) * (ocs2 - b)
                    + Math.pow((1.0 - SIGN * ocs2), d) * (ocs2 + b);

            double dfdcs2 = l * v * Math.pow((1.0 + SIGN * ocs2), d)
                    * (1.0 + SIGN * d * (ocs2 - b) / (1.0 + SIGN * ocs2))
                    + Math.pow((1.0 - SIGN * ocs2), d)
                    * (1.0 - SIGN * d * (ocs2 + b) / (1.0 - SIGN * ocs2));

            cs2 = ocs2 - fcs2 / dfdcs2;
            diff = Math.abs(cs2 - ocs2) / ocs2;
            ocs2 = cs2;
        }

        node.vel = (xi + SIGN * cs2) / (1.0 + SIGN * xi * cs2);

        double tmp1 = (cs2 * cs2 * (gam1 - state.cs * state.cs))
                / (state.cs * state.cs * (gam1 - cs2 * cs2));

        node.rho0 = state.rho0 * Math.pow(tmp1, (1.0 / gam1));

        node.p = cs2 * cs2 * gam1 * node.rho0 / ((gam1 - cs2 * cs2) * gamma);

        node.u = node.p / (gam1 * node.rho0);
    }
///=========================================================

    public void Base_Integrals(eqs eqs, final int n, final double x_left, final double x_right, evolutionVariables ivars, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        base_integration(eqs, n, x_left, x_right, ivars, lstate, rstate, lwave, rwave);
    }

    public void Integrals_Past(evolutionVariables ivars, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        base_integration_past(x1, x5, ivars, lstate, rstate, lwave, rwave);
    }

    public void Left_Integrals_Past(evolutionVariables ivars, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        base_integration_past(x1, x_discontinuity, ivars, lstate, rstate, lwave, rwave);
    }

    public void Right_Integrals_Past(evolutionVariables ivars, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        base_integration_past(x_discontinuity, x5, ivars, lstate, rstate, lwave, rwave);
    }
///=========================================================

    /**
     * Integrals between x_left and x_right
     *
     * Since it is quite an effort to compute the values all integrals are done
     * in one sweep.
     *
     * Uses a simple Trapezoidal rule for now. Might be good to improve this.
     *
     * \param[in] n number of integration points (n-1 = number of integration
     * intervals) \param[in] x_left left boundary \param[in] x_right right
     * boundary \param[out] integral Values of the integrals (S,D,E)
     */
    public void base_integration(eqs eos, final int n, final double x_left, final double x_right, evolutionVariables integral, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        final double range = x_right - x_left;
        final double deltax = range / (double) (n - 1);

        evolutionVariables work_sum = new evolutionVariables();
        /// Interior values

        exactNode node = new exactNode();
        for (int i = 1; i < (n - 1); ++i) {
            double x = x_left + deltax * (double) (i);
            GetExactValues(x, eos, node, lstate, rstate);
            work_sum.S += node.S;
            work_sum.D += node.D;
            work_sum.E += node.E;
        }
        work_sum.multiplyByScalar(2.0);

        /// Boundary values

        GetExactValues(x_left, eos, node, lstate, rstate);
        work_sum.S += node.S;
        work_sum.D += node.D;
        work_sum.E += node.E;

        GetExactValues(x_right, eos, node, lstate, rstate);
        work_sum.S += node.S;
        work_sum.D += node.D;
        work_sum.E += node.E;

        work_sum.multiplyByScalar((range / (double) (2 * (n - 1))));
        integral = work_sum;
        return;
    }
///-----------------------------------------------------------------------

    /**
     * Integrals between x_1 and x_discontinuity
     *
     * \param[in] n number of integration points \param[out] integral Values of
     * the integrals (S,D,E)
     */
    public void Integrals(eqs eqs, final int n, evolutionVariables integral, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        final double x_left = x1;
        final double x_right = x5;
        base_integration(eqs, n, x_left, x_right, integral, lstate, rstate, lwave, rwave);
    }

    public void Left_Integrals(eqs eqs, final int n, evolutionVariables integral, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        final double x_left = x1;
        final double x_right = x_discontinuity;
        base_integration(eqs, n, x_left, x_right, integral, lstate, rstate, lwave, rwave);
        return;
    }

    public void Right_Integrals(eqs eqs, final int n, evolutionVariables integral, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        final double x_left = x_discontinuity;
        final double x_right = x5;
        base_integration(eqs, n, x_left, x_right, integral, lstate, rstate, lwave, rwave);
        return;
    }
///====================================================================
/// Integration on past (initial) time slice

    public void base_integration_past(final double x_left, final double x_right, evolutionVariables integral, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        if (x_left >= x_right) {
            System.err.println("RiemannExact::base_integration_past: ERROR: x_left = " + x_left + " >= " + x_right + "\n");
            System.exit(1);
        }

        if (x_right <= x_discontinuity) {
            final double dx = x_right - x_left;
            evolutionVariables evars = new evolutionVariables(lstate);
            integral.S = evars.S * dx;
            integral.D = evars.D * dx;
            integral.E = evars.E * dx;
            return;
        } else if (x_left >= x_discontinuity) {
            final double dx = x_right - x_left;
            evolutionVariables evars = new evolutionVariables(rstate);
            integral.S = evars.S * dx;
            integral.D = evars.D * dx;
            integral.E = evars.E * dx;
            return;
        } else {
            final double dx_left = x_discontinuity - x_left;
            final double dx_right = x_right - x_discontinuity;
            evolutionVariables eleft = new evolutionVariables(lstate);
            evolutionVariables eright = new evolutionVariables(rstate);
            integral.S = eleft.S * dx_left + eright.S * dx_right;
            integral.D = eleft.D * dx_left + eright.D * dx_right;
            integral.E = eleft.E * dx_left + eright.E * dx_right;
            return;
        }
    }
///====================================================================
/// Make fluxes

    public void MakeFlux(eqs eos, final double x, evolutionVariables flux, stateVariables lstate, stateVariables rstate, waveVariables lwave, waveVariables rwave) {
        exactNode exact = new exactNode();
        GetExactValues(x, eos, exact, lstate, rstate);
        flux.S = exact.S * exact.vel + exact.p;
        flux.D = exact.D * exact.vel;
        flux.E = exact.E * exact.vel + exact.p * exact.vel;
    }
///====================================================================
/// Make evolution variables from primitive variables

    public void makeEvoVars(evolutionVariables ev, stateVariables s) {
        ev.S = s.rho0 * s.h * s.w * s.w * s.vel;
        ev.D = s.rho0 / Math.sqrt(1.0 - s.vel * s.vel);
        ev.E = s.rho0 * s.h * s.w * s.w - s.p - ev.D;  // D must already be made
    }
///====================================================================

    /**
     * Compute solution on a simple 1d mesh
     *
     * \retval ExactNode* pointer to array of ExactNodes containing the mesh
     */
    public ArrayList<exactNode> MakeMeshValues(eqs eos, final int n, final double x_left, 
            final double x_right, stateVariables lstate, stateVariables rstate) {
        if (n <= 5) {
            System.err.println("RiemannExact::MakeMeshValues: ERROR: b too small\n n = " + n + "\n");
            System.exit(1);
        }

        if (x_right <= x_left) {
            System.err.println("RiemannExact::MakeMeshValues: ERROR: x_right <= x_left\n x_left = " + x_left + " x_right = " + x_right + "\n");
            System.exit(1);
        }
        /**
         * This array list is public within an exactNode. Even though it looks
         * like it constructs an array list of size "n", it does not. The "for"
         * loop does this.
         */
        ArrayList<exactNode> nlist = new ArrayList<exactNode>(n);
        for (int i = 0; i < n; i++) {
            exactNode a = new exactNode();
            nlist.add(a);
        }
        /// Values at nodes for graphing purposes
        final double dx = (x_right - x_left) / ((double) (n - 1));

        for (int i = 0; i < n; i++) {
            double x = x_left + ((double) (i)) * dx;
            GetExactValues(x, eos, nlist.get(i), lstate, rstate);
        }

        return nlist;
    }
}