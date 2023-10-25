package Gregory_Mohammed_Simulation_Code_Masters;

import java.util.ArrayList;

public class cell {

    /**
     * ============================================== 
     * This class contains the
     * variables and methods needed to generate and keep track of the cells and
     * any process related to them.
     * ==============================================
     */
    /**
     * ============================================== 
     * Private variables.
     * ==============================================
     */
    baseVariables old;
    baseVariables cur;
    /**
     * ---------------------------------------------- 
     * An id to keep tabs on a
     * cell. 
     * ----------------------------------------------
     */
    public int _id;
    /**
     * ============================================== 
     * private variables.
     * ==============================================
     */
    public double _x;

    /**
     * ============================================== 
     * Accessors.
     * ==============================================
     */
    public double xmid(int i, ArrayList<node> nlist) {
        return (0.5 * (nlist.get(i).x + nlist.get(i + 1).x));
    }

    public double Volume(int i, ArrayList<node> nlist) {
        return (nlist.get(i + 1).x - nlist.get(i).x);
    }

    public double H(int i, ArrayList<node> nlist) {
        return (nlist.get(i + 1).x - nlist.get(i).x);
    }

    public double RestMass() {
        return (cur._D / cur._a);
    }
    
    public double shockTubeLength(ArrayList<node> nlist) {
        return (nlist.get(nlist.size()-1).x - nlist.get(0).x);
    }

    public double attenuate(double currentCellLength, ArrayList<node> nlist) {
        double tmp = shockTubeLength(nlist);
        
        return (currentCellLength / tmp);
    }
    
    public node nleft(int i, ArrayList<node> nlist) {
        node nleft = nlist.get(i);
        return nleft;
    }

    public node nright(int i, ArrayList<node> nlist) {
        node nright = nlist.get(i + 1);
        return nright;
    }

    /**
     * ============================================== 
     * Constructor
     * ==============================================
     */
    public cell() {
        _id = -1;
        _x = 0;

        old = new baseVariables();
        cur = new baseVariables();
    }

    /**
     * ============================================== 
     * The various methods.
     * ==============================================
     */
    /**
     * ============================================== 
     * Makes. Note: p must be
     * made by a call to the equation of state class
     * ==============================================
     */
    public void Make(double v_in, double rho0_in, double eps_in) {
        cur._v = v_in;
        cur._rho0 = rho0_in;
        cur._eps = eps_in;
    }

    public void makeP(eqs eqs) {
        cur._p = eqs.P(cur._rho0, cur._eps);
    }

    public void makeA() {
        cur._a = 1.0 / Math.sqrt(1.0 - cur._v * cur._v);
    }

    /**
     * ============================================== 
     * Upwind functions.
     * ==============================================
     */
    public void makePrimitiveOld(eqs eqs) {
        primitiveSolve psolve = new primitiveSolve();

        int iError = old.makePrimitive(eqs, psolve);
        /**
         * ========================================== 
         * Test for bad values in the
         * old set. 
         * ==========================================
         */
        if (iError != 0) {
            System.err.println("cell :: ERROR :: makePrimitiveOld failed to solve for primitive variables.");
            System.exit(1);
        }
    }

    public void makePrimitiveCur(eqs eqs) {
        primitiveSolve psolve = new primitiveSolve();

        int iError = cur.makePrimitive(eqs, psolve);
        /**
         * ========================================== 
         * Test for bad values in the
         * new set. 
         * ==========================================
         */
        if (iError != 0) {
            System.err.println("cell :: ERROR :: makePrimitiveCur failed to solve for primitive variables.");
            System.exit(1);
        }
    }

    public void makeConserved(eqs eqs) {
        makeA();
        final double w = 1.0 + eqs.gamma * cur._eps;
        final double tmp1 = cur._rho0 * w * cur._a * cur._a;
        cur._S = tmp1 * cur._v;
        cur._D = cur._a * cur._rho0;
        cur._E = tmp1 - cur._p - cur._D;
        cur._p = eqs.P(cur._rho0, cur._eps);
        /**
         * ========================================== 
         * Test for negative energy.
         * ==========================================
         */
        if (cur._E < 0.0) {
            System.err.println("cell :: ERROR :: makeConserved generates negative energy.\n"
                    + "eqs.gamma = " + eqs.gamma + " :: cur._eps = " + cur._eps + " :: cur._p = " + cur._p + " :: w = " + w + "\n"
                    + "cur._rho0 = " + cur._rho0 + " :: cur._a = " + cur._a + " :: cur._v = " + cur._v + " :: tmp1 = " + tmp1 + "\n"
                    + "cur._S = " + cur._S + " :: cur._D = " + cur._D + " :: cur._E = " + cur._E);
            System.exit(1);
        }
    }

    /**
     * ============================================== 
     * Evolve
     * ==============================================
     */
    public void evolveSimplePredictor(int i, neutrinoMethods neutrinoM, ArrayList<node> nlist, timeData ptime, IO io) {
        double dtdx = ptime.delta / H(i, nlist);

        if (io.neutrinos == 0) {
            cur._S = old._S - (dtdx * (nlist.get(i + 1).oflux._S - nlist.get(i).oflux._S));
            cur._D = old._D - (dtdx * (nlist.get(i + 1).oflux._D - nlist.get(i).oflux._D));
            cur._E = old._E - (dtdx * (nlist.get(i + 1).oflux._E - nlist.get(i).oflux._E));
        } else {
            double Qt = neutrinoM.sourceQt(io, old._a, _x, old._v, old._rho0);
            double Qx = neutrinoM.sourceQx(io, old._a, _x, old._v, old._rho0);
            
            cur._S = old._S + (Qx * ptime.delta) - (dtdx * (nlist.get(i + 1).oflux._S - nlist.get(i).oflux._S));
            cur._D = old._D - (dtdx * (nlist.get(i + 1).oflux._D - nlist.get(i).oflux._D));
            cur._E = old._E - (Qt * ptime.delta) - (dtdx * (nlist.get(i + 1).oflux._E - nlist.get(i).oflux._E));
        }
    }

    public void evolveSimpleCorrector(int i, neutrinoMethods neutrinoM, ArrayList<node> nlist, timeData ptime, IO io) {
        double dtdx = ptime.delta / H(i, nlist);

        if (io.neutrinos == 0) {
            cur._S = old._S - (dtdx * (((nlist.get(i + 1).nflux._S + nlist.get(i + 1).oflux._S) / 2.0)
                    - ((nlist.get(i).nflux._S + nlist.get(i).oflux._S) / 2.0)));
            cur._D = old._D - (dtdx * (((nlist.get(i + 1).nflux._D + nlist.get(i + 1).oflux._D) / 2.0)
                    - ((nlist.get(i).nflux._D + nlist.get(i).oflux._D) / 2.0)));
            cur._E = old._E - (dtdx * (((nlist.get(i + 1).nflux._E + nlist.get(i + 1).oflux._E) / 2.0)
                    - ((nlist.get(i).nflux._E + nlist.get(i).oflux._E) / 2.0)));
        } else {
            double Qt = neutrinoM.sourceQt(io, cur._a, _x, cur._v, cur._rho0);
            double Qx = neutrinoM.sourceQx(io, cur._a, _x, cur._v, cur._rho0);

            cur._S = old._S + (Qx * ptime.delta) - (dtdx * (((nlist.get(i + 1).nflux._S + nlist.get(i + 1).oflux._S) / 2.0)
                    - ((nlist.get(i).nflux._S + nlist.get(i).oflux._S) / 2.0)));
            cur._D = old._D - (dtdx * (((nlist.get(i + 1).nflux._D + nlist.get(i + 1).oflux._D) / 2.0)
                    - ((nlist.get(i).nflux._D + nlist.get(i).oflux._D) / 2.0)));
            cur._E = old._E - (Qt * ptime.delta) - (dtdx * (((nlist.get(i + 1).nflux._E + nlist.get(i + 1).oflux._E) / 2.0)
                    - ((nlist.get(i).nflux._E + nlist.get(i).oflux._E) / 2.0)));
        }
    }

    /**
     * ============================================== 
     * Process for exchanging the
     * old for the new. 
     * ==============================================
     */
    public void cellExchange() {
        old._S = cur._S;
        old._D = cur._D;
        old._E = cur._E;

        old._p = cur._p;
        old._eps = cur._eps;
        old._a = cur._a;
        old._rho0 = cur._rho0;
        old._v = cur._v;

        old._sigmaD = cur._sigmaD;
        old._sigmaE = cur._sigmaE;
        old._sigmaS = cur._sigmaS;

        old._sigmarho0 = cur._sigmarho0;
        old._sigmap = cur._sigmap;
        old._sigmaeps = cur._sigmaeps;
        old._sigmav = cur._sigmav;
    }
}
