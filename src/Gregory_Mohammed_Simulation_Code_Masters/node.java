package Gregory_Mohammed_Simulation_Code_Masters;

import java.util.ArrayList;

public class node {

    /**
     * ==============================================
     * This class contains the variables and methods
     * needed to generate and keep track of the nodes
     * and any process related to them.
     * ==============================================
     */
    /// ==============================================
    flux oflux;
    flux nflux;
    exactNode enode;
    /// ==============================================
    /**
     * ----------------------------------------------
     * An id to keep tabs on a node. Also its position.
     * ----------------------------------------------
     */
    public int id;
    public double x;
    public double q;
    public int solveOk;
    public double aminus;
    public double aplus;

    /**
     * ==============================================
     * Constructor
     * ==============================================
     */
    public node() {
        id = 0;
        q = 0.0;
        /// ==============================================
        oflux = new flux();
        nflux = new flux();
        enode = new exactNode();
        ///=========================================================
        solveOk = 1;
        aminus = 0.0;
        aplus = 0.0;
    }
    
    public cell cleft(int i, ArrayList<cell> clist) {
        cell cleft = clist.get(i - 1);
        return cleft;
    }

    public cell cright(int i, ArrayList<cell> clist) {
        cell cright = clist.get(i);
        return cright;
    }

    /**
     * ==============================================
     * For the fluxes ...
     * ==============================================
     */
    /**
     * ----------------------------------------------
     * This method will be used to make the fluxes at
     * the current node. Remember, the Finite Volume
     * Method is cell-centric, so it makes sense for
     * the cell to inherit the methods and implement
     * them.
     * ----------------------------------------------
     */
    public void makeFlux(double leftx, double rightx, int sign, artVisControl artvis, eqs eqs,
            baseVariables left, baseVariables right, timeData timeData, flux flux) {
        /**
         * Primitive variables are simply averaged.
         * Generally these do not need the high-resolution approach.
         */
        double v = 0.5 * (left._v + right._v);
        double p = 0.5 * (left._p + right._p);
        double eps = 0.5 * (left._eps + right._eps);
        /**
         * High resolution evolution variables
         */
        final double deltat = timeData.delta;

        double D, S, E;

        if ((v >= 0.0)) {
            S = left._S + left._sigmaS * (x - leftx + sign * v * deltat / 2.0);
            D = left._D + left._sigmaD * (x - leftx + sign * v * deltat / 2.0);
            E = left._E + left._sigmaE * (x - leftx + sign * v * deltat / 2.0);
        } else { // v<0
            S = right._S + right._sigmaS * (x - rightx + sign * v * deltat / 2.0);
            D = right._D + right._sigmaD * (x - rightx + sign * v * deltat / 2.0);
            E = right._E + right._sigmaE * (x - rightx + sign * v * deltat / 2.0);
        }

        /**
         * now calculate the artificial viscosity contribution
         */
        final double dv = right._v - left._v;
        double ba = 1.0 / Math.sqrt(1.0 - v * v);
        final double cs = eqs.CS(eps);

        if (dv < 0.0) {
            final double absdv = Math.abs(dv);
            q = (D + eps * D + p * ba) * absdv * (artvis._k1 * cs + artvis._k2 * absdv);
        } else {
            q = 0.0;
        }

        /**
         * Compute fluxes
         */
        flux.make(S, D, E, v, p, q);
    }

    /**
     * ----------------------------------------------
     * This method predicts the flux on the node, before
     * the cell values can be calculated.
     * ----------------------------------------------
     */
    public void makeFluxPredictor(int i, artVisControl artvis, eqs eqs, ArrayList<cell> clist, timeData ptime) {
        try {
            final int ieLeft = i - 1;
            final int ieRight = i;
            if (ieLeft == -1) {
                q = 0.0;
                oflux.make(clist.get(i).old._S, clist.get(i).old._D, clist.get(i).old._E,
                        clist.get(i).old._v, clist.get(i).old._p, q);
                return;
            }
            if (ieRight == clist.size()) {
                q = 0.0;
                oflux.make(clist.get(i - 1).old._S, clist.get(i - 1).old._D, clist.get(i - 1).old._E,
                        clist.get(i - 1).old._v, clist.get(i - 1).old._p, q);
                return;
            }
            final int sign = -1;  // forward differencing (predictor)

            makeFlux(clist.get(i - 1)._x, clist.get(i)._x, sign, artvis, eqs, clist.get(i - 1).old, clist.get(i).old, ptime, oflux);

        } catch (IndexOutOfBoundsException ioobe) {
            System.err.println("node :: makeFluxPredictor :: message from IndexOutOfBoundsException :: ioobe\n" + ioobe);
            System.exit(1);
        }
    }

    /**
     * ----------------------------------------------
     * This method corrects the flux on the node, that
     * is, it gives the correct flux and does away with
     * the predictions.
     * ----------------------------------------------
     */
    public void makeFluxCorrector(int i, artVisControl artvis, eqs eqs, ArrayList<cell> clist, timeData ptime) {
        try {
            final int ieLeft = i - 1;
            final int ieRight = i;
            if (ieLeft == -1) {
                q = 0.0;
                nflux.make(clist.get(i).cur._S, clist.get(i).cur._D, clist.get(i).cur._E,
                        clist.get(i).cur._v, clist.get(i).cur._p, q);
                return;
            }
            if (ieRight == clist.size()) {
                q = 0.0;
                nflux.make(clist.get(i - 1).cur._S, clist.get(i - 1).cur._D, clist.get(i - 1).cur._E,
                        clist.get(i - 1).cur._v, clist.get(i - 1).cur._p, q);
                return;
            }

            // Compute fluxes

            //const int sign = +1;   // backward differencing (corrector)
            final int sign = +1;   // backward differencing (corrector)

            flux curFlux = new flux();
            makeFlux(clist.get(i - 1)._x, clist.get(i)._x, sign, artvis, eqs, clist.get(i - 1).cur, clist.get(i).cur, ptime, curFlux);

            nflux._S = 0.5 * (curFlux._S + oflux._S);
            nflux._D = 0.5 * (curFlux._D + oflux._D);
            nflux._E = 0.5 * (curFlux._E + oflux._E);
        } catch (IndexOutOfBoundsException ioobe) {
            System.err.println("node :: makeFluxCorrector :: message from IndexOutOfBoundsException :: ioobe\n" + ioobe);
            System.exit(1);
        }
    }

    /**
     * 
     * @param left
     * @param right
     * @param leftx
     * @param rightx
     * @param eqs
     * @param ptime
     * @param artvis
     * @param sign
     * @param flux 
     * 
     * This code (MakeFluxUpOrRie) was written by PJM in C++.
     * Modified for Java by GM.
     */
    public void MakeFluxUpOrRie(ArrayList<cell> clist, ArrayList<node> nlist, baseVariables left, baseVariables right, 
            timeData ptime, artVisControl artvis, int sign, flux flux, int i, eqs eqs, IO io) {

        // Averages are used for various tests and upwind fallbacks, so compute them all here. 

        final double v = 0.5 * (left._v + right._v);
        final double p = 0.5 * (left._p + right._p);
        final double rho0 = 0.5 * (left._rho0 + right._rho0);
        final double eps = 0.5 * (left._eps + right._eps);
        final double S = 0.5 * (left._S + right._S);
        final double D = 0.5 * (left._D + right._D);
        final double E = 0.5 * (left._E + right._E);

        final double a = 1.0 / Math.sqrt(1.0 - v * v);
        final double cs = eqs.CS(eps);
        final double dv = right._v - left._v;

        // Assuming rho0 is never zero this gives a reasonable characteristic p for convergence/error testing


        final double p_characteristic = io.k * Math.pow(rho0, eqs.gamma);
        final double eps_characteristic = p_characteristic / ((eqs.gamma - 1.0) * rho0);
        final double w_characteristic = 1.0 + (eqs.gamma * eps_characteristic);
        final double E_characteristic = (rho0*w_characteristic*Math.pow(a, 2.0)) - p_characteristic - (a*rho0);
        final double S_characteristic = rho0*w_characteristic*Math.pow(a, 2.0); /// Assumes v = 1.

        // HIGH RESOLUTION FITTING
        // -----------------------
        // BaseVariables is actually used for the Element class (current and old) but
        // also has a built-in primitive solver.  So it's handy to use it here even
        // though other element-specific members are unused.

        baseVariables highres_left = new baseVariables();
        baseVariables highres_right = new baseVariables();
        MakeHighRes(clist, nlist, i, io, highres_left, highres_right);



        //========================================================================
        // Simple v-based Upwind.  Used as a fallback if Riemann throws an exception

        double Sup, Dup, Eup, vup, pup;
        if (v >= 0) {
            Sup = highres_left._S;
            Dup = highres_left._D;
            Eup = highres_left._E;
            vup = highres_left._v;
            pup = highres_left._p;
        } else {
            Sup = highres_right._S;
            Dup = highres_right._D;
            Eup = highres_right._E;
            vup = highres_right._v;
            pup = highres_right._p;
        }

        // Artificial viscosity.  Necessary if straight upwind is used.

        if (dv < 0.0) {
            final double absdv = Math.abs(dv);
            q = (Dup + eps * Dup + p * a) * absdv * (artvis._k1 * cs + artvis._k2 * absdv);
        } else {
            q = 0.0;
        }

        if (io.method == 0) {
            flux.make(Sup, Dup, Eup, vup, pup, q);
            flux.make(Sup, Dup, Eup, v, p, q);
            return;
        }
        //=====================================================================================
        // If Right ~= Left then just use upwind and don't bother with Riemann

        // **** NOTE: should consider whether to use primitive variables here.

        // Check to see if these are almost the same.
        // If they are then just use the upwind.
        // Note that v is compared to 1 (light speed) and h compared to 1 (1+..)
        // rho0 is never zero so compare to the mean

        final double diff_S = highres_left._S - highres_right._S;
        final double diff_D = highres_left._D - highres_right._D;
        final double diff_E = highres_left._E - highres_right._E;

        final double diff = Math.abs(diff_S / S_characteristic) + Math.abs(diff_D / D) + Math.abs(diff_E / E_characteristic);

        final double EPS = 1.0e-8;
        if ((diff < EPS) || (p < EPS * p_characteristic) || (eps < EPS * eps_characteristic)) {
            flux.make(Sup, Dup, Eup, v, p, q);
            return;
        }

        //======================================================================================
        // Ok, let's try Riemann.

        // construct left and right states for Riemann solver

        stateVariables state_left = new stateVariables(eqs, highres_left._rho0, highres_left._p, highres_left._v);

        stateVariables state_right = new stateVariables(eqs, highres_right._rho0, highres_right._p, highres_right._v);

        // Compute Riemann solution
        // -note that position is not significant, just the difference between states.

        final double x_discontinuity = 0.0;

        // Evolution time is NOT SIGNIFICANT!
        // -Only interested in the (constant) solution along the x=x_node ray.
        // -So evolve for any non-zero time, and return the variables on that ray
        // -I've tested this by running with different time steps.  No difference.

        final double evolution_time = ptime.delta;

        // Compute all variable values on interface from Riemann solution
        // These should be constant along the interface ray

        try {

            riemannExact rExact = new riemannExact();

            rExact.Solve(evolution_time, x_discontinuity, eqs, state_left, state_right);


            rExact.GetExactValues(x_discontinuity, eqs, enode, state_left, state_right);

            flux.make(enode.S, enode.D, enode.E, enode.vel, enode.p, q);


        } catch (Exception straddleError) {

            //System.err.println("Node :: MakeFluxUpOrRie: WARNING: id = " + id);
            //System.out.println("\nIt's raining " + straddleError + "!");

            flux.make(Sup, Dup, Eup, v, p, q);
            return;

        }
    }

    ///=============================================================================
    public void makeRiemannFluxPredictor(int i, artVisControl artvis, ArrayList<cell> clist,
            ArrayList<node> nlist, eqs eqs, timeData ptime, IO io) {
        ///=============================================================================
        /**
        ///=====================================================
        Provide boundary handling.
        ///=====================================================
         */
        final int ieLeft = i - 1;
        final int ieRight = i;
        
        if (ieLeft == -1) {
            q = 0.0;
            oflux.make(clist.get(ieRight).old._S, clist.get(ieRight).old._D,
                    clist.get(ieRight).old._E, clist.get(ieRight).old._v, clist.get(ieRight).old._p, q);

            return;
        }
        if (ieRight == clist.size()) {  
            q = 0.0;
            oflux.make(clist.get(ieLeft).old._S, clist.get(ieLeft).old._D,
                    clist.get(ieLeft).old._E, clist.get(ieLeft).old._v, clist.get(ieLeft).old._p, q);

            return;
        }
        
        final cell left = clist.get(i - 1);
        final cell right = clist.get(i);
        final int sign = -1;
        
        MakeFluxUpOrRie(clist, nlist, left.old, right.old, ptime, artvis, sign, oflux, i, eqs, io);
    }

    ///=============================================================================
    public void makeRiemannFluxCorrector(int i, artVisControl artvis, ArrayList<cell> clist,
            ArrayList<node> nlist, eqs eqs, timeData ptime, IO io) {
        ///=============================================================================
        /**
        ///=====================================================
        Provide boundary handling.
        ///=====================================================
         */
        final int ieLeft = i - 1;
        final int ieRight = i;

        if (ieLeft == -1) {
            q = 0.0;
            nflux.make(clist.get(ieRight).old._S, clist.get(ieRight).old._D, clist.get(ieRight).old._E, clist.get(ieRight).old._v, clist.get(ieRight).old._p, q);

            return;
        }
        if (ieRight == clist.size()) {  //PJM:  should be ..-1???  But this works here.
            q = 0.0;
            nflux.make(clist.get(ieLeft).old._S, clist.get(ieLeft).old._D, clist.get(ieLeft).old._E, clist.get(ieLeft).old._v, clist.get(ieLeft).old._p, q);

            return;
        }
        
        final cell left = clist.get(i - 1);
        final cell right = clist.get(i);
        final int sign = +1;
        
        MakeFluxUpOrRie(clist, nlist, left.cur, right.cur, ptime, artvis, sign, nflux, i, eqs, io);
    }

    //============================================================================
    /** Make the high resolution approximations on the right and left of the current Node
    
    \param[in] reconstruction_type   Switch for linear or cubic reconstruction
    \param[out] highres_left    High resolution base variables on left
    \param[out] highres_right   High resolution base variables on right 
     */
    public void MakeHighRes(ArrayList<cell> clist, ArrayList<node> nlist, int i, IO io, 
            baseVariables highres_left, baseVariables highres_right) {
        /**
         * Reconstruction Types:
         * Constant = 0
         * Linear = 1
         * CubicHermite = 2
         */

        if (io.reconstructionType == 0) {
            highres_left.MakeHighResConstant(cleft(i, clist), x);
            highres_right.MakeHighResConstant(cright(i, clist), x);
            return;
        }

        if (io.reconstructionType == 1) {
            highres_left.MakeHighResLinear(cleft(i, clist), x);
            highres_right.MakeHighResLinear(cright(i, clist), x);
            return;
        }

        if (io.reconstructionType != 2) {
            System.err.println("Node :: MakeHighRes: ERROR: unknown reconstruction type.\n"
                    + "Reconstruction_type = " + io.reconstructionType);
            System.exit(1);
        }

        // Elements

        cell cell1 = new cell();
        cell cell2 = new cell();  // will be ordered, elem1 to the left of elem2

        // Left

        cell2 = cleft(i, clist);
        
        cell1 = cell2.nleft(i, nlist).cleft(i, clist); // Left element pointing to lefter node pointing to lefter element

        if (cell1 == null) {
            highres_left.MakeHighResLinear(cell2, x);
        } else {
            highres_left.MakeHighResCubicHermite(cell1, cell2, x);
        }

        // Right

        cell1 = cright(i, clist);
        
        cell2 = cell1.nright(i, nlist).cright(i, clist); 

        if (cell2 == null) {
            highres_right.MakeHighResLinear(cell1, x);
        } else {
            highres_right.MakeHighResCubicHermite(cell1, cell2, x);
        }
    }
    
    /**
     * =====================================================
     * This is the function which loads the left and right
     * state variables for the Riemann solver to work on.
     * =====================================================
     */
    public void loadAllOldStateVars(eqs eqs, int i, ArrayList<cell> clist, stateVariables lstate, stateVariables rstate) {
        final int ieLeft = i - 1;
        final int ieRight = i;

        if (ieLeft == -1) {
            ///=====================================================
            /**
            At node 0, let the Riemann solver know what to do for
            left = 0.
             */
            ///=====================================================
            lstate.rho0 = clist.get(i).old._rho0;
            lstate.p = clist.get(i).old._p;
            lstate.vel = clist.get(i).old._v;
            lstate.u = clist.get(i).old._eps;

            lstate.MakeCS(eqs);
            lstate.MakeH();
            lstate.w = clist.get(i).old._a;
            lstate.icell = -1;  // left side
            ///=====================================================
            /// Load the rstate ...
            ///=====================================================
            rstate.rho0 = clist.get(i).old._rho0;
            rstate.p = clist.get(i).old._p;
            rstate.vel = clist.get(i).old._v;
            rstate.u = clist.get(i).old._eps;

            rstate.MakeCS(eqs);
            rstate.MakeH();
            rstate.w = clist.get(i).old._a;
            rstate.icell = ieRight;
            ///=====================================================
            return;
        }
        if (ieRight == clist.size()) {
            ///=====================================================
            /// Load the lstate ...
            ///=====================================================
            lstate.rho0 = clist.get(i - 1).old._rho0;
            lstate.p = clist.get(i - 1).old._p;
            lstate.vel = clist.get(i - 1).old._v;
            lstate.u = clist.get(i - 1).old._eps;

            lstate.MakeCS(eqs);
            lstate.MakeH();
            lstate.w = clist.get(i - 1).old._a;
            lstate.icell = ieLeft;
            ///=====================================================
            /**
            At node N, let the Riemann solver know what to do for
            right = 0.
             */
            ///=====================================================
            /// Load the rstate ...
            ///=====================================================
            rstate.rho0 = clist.get(i - 1).old._rho0;
            rstate.p = clist.get(i - 1).old._p;
            rstate.vel = clist.get(i - 1).old._v;
            rstate.u = clist.get(i - 1).old._eps;

            rstate.MakeCS(eqs);
            rstate.MakeH();
            rstate.w = clist.get(i - 1).old._a;
            rstate.icell = -ieRight;
            ///=====================================================
            return;
        }
        ///=====================================================
        /// Load the lstate ...
        ///=====================================================
        baseVariables old = clist.get(i - 1).old;

        double halfdeltax = 0.5 * (clist.get(i)._x - clist.get(i - 1)._x);

        double hrrho0l = old._rho0 + (halfdeltax * old._sigmarho0);
        double hrpl = old._p + (halfdeltax * old._sigmap);
        double hrepsl = old._eps + (halfdeltax * old._sigmaeps);
        double hrvl = old._v + (halfdeltax * old._sigmav);

        lstate.rho0 = hrrho0l;
        lstate.p = hrpl;
        lstate.vel = hrvl;
        lstate.u = hrepsl;

        lstate.MakeCS(eqs);
        lstate.MakeH();
        lstate.w = 1.0 / Math.sqrt(1.0 - Math.pow(hrvl, 2.0));
        lstate.icell = ieLeft;
        ///=====================================================
        /// Load the rstate ...
        ///=====================================================
        old = clist.get(i).old;

        double hrrho0r = old._rho0 - (halfdeltax * old._sigmarho0);
        double hrpr = old._p - (halfdeltax * old._sigmap);
        double hrepsr = old._eps - (halfdeltax * old._sigmaeps);
        double hrvr = old._v - (halfdeltax * old._sigmav);

        rstate.rho0 = hrrho0r;
        rstate.p = hrpr;
        rstate.vel = hrvr;
        rstate.u = hrepsr;

        rstate.MakeCS(eqs);
        rstate.MakeH();
        rstate.w = 1.0 / Math.sqrt(1.0 - Math.pow(hrvr, 2.0));
        rstate.icell = ieRight;
        ///=====================================================
    }

    public void loadAllCurStateVars(eqs eqs, int i, ArrayList<cell> clist, stateVariables lstate, stateVariables rstate) {
        final int ieLeft = i - 1;
        final int ieRight = i;

        if (ieLeft == -1) {
            ///=====================================================
            /**
            At node 0, let the Riemann solver know what to do for
            left = 0.
             */
            ///=====================================================
            lstate.rho0 = clist.get(i).old._rho0;
            lstate.p = clist.get(i).old._p;
            lstate.vel = clist.get(i).old._v;
            lstate.u = clist.get(i).old._eps;

            lstate.MakeCS(eqs);
            lstate.MakeH();
            lstate.w = clist.get(i).old._a;
            lstate.icell = -1;  // left side
            ///=====================================================
            /// Load the rstate ...
            ///=====================================================
            rstate.rho0 = clist.get(i).old._rho0;
            rstate.p = clist.get(i).old._p;
            rstate.vel = clist.get(i).old._v;
            rstate.u = clist.get(i).old._eps;

            rstate.MakeCS(eqs);
            rstate.MakeH();
            rstate.w = clist.get(i).old._a;
            rstate.icell = ieRight;
            ///=====================================================
            return;
        }
        if (ieRight == clist.size()) {
            ///=====================================================
            /// Load the lstate ...
            ///=====================================================
            lstate.rho0 = clist.get(i - 1).old._rho0;
            lstate.p = clist.get(i - 1).old._p;
            lstate.vel = clist.get(i - 1).old._v;
            lstate.u = clist.get(i - 1).old._eps;

            lstate.MakeCS(eqs);
            lstate.MakeH();
            lstate.w = clist.get(i - 1).old._a;
            lstate.icell = ieLeft;
            ///=====================================================
            /**
            At node N, let the Riemann solver know what to do for
            right = 0.
             */
            ///=====================================================
            /// Load the rstate ...
            ///=====================================================
            rstate.rho0 = clist.get(i - 1).old._rho0;
            rstate.p = clist.get(i - 1).old._p;
            rstate.vel = clist.get(i - 1).old._v;
            rstate.u = clist.get(i - 1).old._eps;

            rstate.MakeCS(eqs);
            rstate.MakeH();
            rstate.w = clist.get(i - 1).old._a;
            rstate.icell = -ieRight;
            ///=====================================================
            return;
        }
        ///=====================================================
        /// Load the lstate ...
        ///=====================================================
        double halfdeltax = 0.5 * (clist.get(i)._x - clist.get(i - 1)._x);

        baseVariables cur = clist.get(i - 1).cur;

        double hrrho0l = cur._rho0 + (halfdeltax * cur._sigmarho0);
        double hrpl = cur._p + (halfdeltax * cur._sigmap);
        double hrepsl = cur._eps + (halfdeltax * cur._sigmaeps);
        double hrvl = cur._v + (halfdeltax * cur._sigmav);

        lstate.rho0 = hrrho0l;
        lstate.p = hrpl;
        lstate.vel = hrvl;
        lstate.u = hrepsl;

        lstate.MakeCS(eqs);
        lstate.MakeH();
        lstate.w = 1.0 / Math.sqrt(1.0 - Math.pow(hrvl, 2.0));
        lstate.icell = ieLeft;
        ///=====================================================
        /// Load the rstate ...
        ///=====================================================
        cur = clist.get(i).cur;

        double hrrho0r = cur._rho0 - (halfdeltax * cur._sigmarho0);
        double hrpr = cur._p - (halfdeltax * cur._sigmap);
        double hrepsr = cur._eps - (halfdeltax * cur._sigmaeps);
        double hrvr = cur._v - (halfdeltax * cur._sigmav);

        rstate.rho0 = hrrho0r;
        rstate.p = hrpr;
        rstate.vel = hrvr;
        rstate.u = hrepsr;

        rstate.MakeCS(eqs);
        rstate.MakeH();
        rstate.w = 1.0 / Math.sqrt(1.0 - Math.pow(hrvr, 2.0));
        rstate.icell = ieRight;
        ///=====================================================
    }
}
