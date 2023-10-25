package Gregory_Mohammed_Simulation_Code_Masters;

import java.util.*;

/**
 * @author Gregory Mohammed
 */
public class grid {

    /**
     * ===================================================================
     */
    ArrayList<node> nlist;
    ArrayList<cell> clist;
    /**
     * ===================================================================
     */
    public double dx_mean;
    public int corrector_level;
    public int n_corrector;
    public int run_id;
    double x_left;
    double x_right;
    double v_left;
    double v_right;
    double rho0_left;
    double rho0_right;
    double eps_left;
    double eps_right;
    double artvisk1;
    double artvisk2;
    public int UpwindorGodunov;

    /**
     * ===================================================================
     *
     * @param n_elements
     * @param x_left
     * @param x_right
     * @return
     * ===================================================================
     */
    public grid() {
        nlist = new ArrayList<node>();
        clist = new ArrayList<cell>();

        dx_mean = 0.0;
        corrector_level = -1;
        n_corrector = 1;
        run_id = 0;
        x_left = -0.5;
        x_right = 0.5;

        v_left = 0.0;
        v_right = 0.0;
        rho0_left = 1.0;
        rho0_right = 0.5;
        eps_left = 0.2;
        eps_right = 0.1;

        artvisk1 = 0.0;
        artvisk2 = 0.0;

        UpwindorGodunov = 0;
    }

    /**
     * ===================================================================
     * Switch to make the initial data as a smoothed shock or a real shock.
     * ===================================================================
     */
    public void MakeInitial(eqs eqs, IO io, greg_const gc, geometrize geo, neutrinoMethods neutrinoM) {

        switch (io.initShock) {
            case 0:
                makeSmoothJump(eqs, io);
                break;
            case 1:
                switch (io.dataSwitch) {
                    case 0:
                        makeShockTube(eqs, io);
                        break;
                    case 1:
                        makeShockTube(eqs, io);
                        break;
                    case 2:
                        makeKKTShockTube(eqs, io, geo);
                        break;
                    default:
                        System.err.println("Grid::MakeInitial: ERROR: unknown data type.");
                        System.exit(1);
                }
                break;
            default:
                System.err.println("Grid::MakeInitial: ERROR: unknown init_type \n");
                System.exit(1);
        }

        gridExchange();
    }

    /**
     * =================================================================== Make
     * the geometry of the grid.
     * ===================================================================
     */
    public int makeGeometry(int n_elements, double x_left, double x_right) {
        if ((Math.abs(x_left - x_right) <= 1.0e-4)) {
            System.err.println("grid::makeGeometry: ERROR: initial left/right data | x_left = " + x_left + " x_right = " + x_right + "\n");
            System.exit(1);
        }
        if (n_elements < 4) {
            System.err.println("grid::makeGeometry: ERROR: n_elements is too small | n_elements = " + n_elements + "\n");
            System.exit(1);
        }

        int n_nodes = n_elements + 1;

        dx_mean = (x_right - x_left) / (double) (n_elements);

        maken(n_nodes);
        makec(n_elements);
        /// ===========================================================================
        /**
         * Positions Make the nodes first and then the elements can be
         * calculated from the nodal positions.
         */
        for (int in = 0; in < n_nodes; ++in) {
            nlist.get(in).x = x_left + dx_mean * (double) (in);
        }

        for (int ie = 0; ie < n_elements; ++ie) {
            clist.get(ie)._x = 0.5 * (nlist.get(ie).x + nlist.get(ie + 1).x);
        }

        return 0;
    }

    /**
     * =================================================================== Make
     * the real shock tube.
     * ===================================================================
     */
    /**
     * =================================================================== Make
     * initial data for shock tube. This is a single call which does everything
     * ===================================================================
     */
    public int makeShockTube(eqs eqs, IO io) {
        /**
         * These values are hard-coded for the Sod shock tube. So, this method
         * cannot be used for any other data. This allows the param.dat file to
         * be flexible for "playing" with data to observe results for different
         * physical situations (or unphysical ones too).
         */
        x_left = io.lx;
        x_right = io.rx;

        v_left = io.lv;
        v_right = io.rv;

        rho0_left = io.lrho0;
        rho0_right = io.rrho0;

        /**
         * 17 March, 2012: Because I edited PJM's data file, param.dat, to have
         * the left and right pressures and densities, now I have to put the
         * calculation for energy here, so that I don't have to change anything
         * else.
         */
        eps_left = eqs.EPS(io.lp, rho0_left); //io.lp / ((io.gamma - 1.0) * rho0_left);
        eps_right = eqs.EPS(io.rp, rho0_right); //io.rp / ((io.gamma - 1.0) * rho0_right);

        /**
         * Make the grid and geometry.
         */
        int ierr = makeGeometry(io.resolution, x_left, x_right);
        if (ierr != 0) {
            System.err.println("grid::makeShockTube: ERROR: makeGeometry returns " + ierr + "\n");
            System.exit(1);
        }

        /**
         * Make the primitive variables
         */
        makeShockTube(eqs, v_left, rho0_left, eps_left,
                v_right, rho0_right, eps_right);

        /**
         * Make the conserved quantities
         */
        makeConserved(eqs);

        return 0;
    }

    /**
     * =================================================================== Make
     * Riemann Shock Tube initial data. Note that only the elements have data.
     * The nodes compute data as required from the adjoining elements.
     * ===================================================================
     */
    public int makeShockTube(eqs eqs, double v_left, double rho0_left, double eps_left,
            double v_right, double rho0_right, double eps_right) {
        /**
         * should check parameters
         */
        if ((rho0_left < 0.0)
                || (rho0_right < 0.0)
                || (eps_left < 0.0)
                || (eps_right < 0.0)) {
            System.err.println("ERROR :: grid :: makeShockTube ::<n"
                    + "Incorrect initial data:\n"
                    + "rho0_left = " + rho0_left + " :: rho0_right = " + rho0_right + "\n"
                    + "eps_left  = " + eps_left + " :: eps_right  = " + eps_right);
            System.exit(1);
        }

        int ne_mid = clist.size() / 2;
        for (int ie = 0; ie < ne_mid; ++ie) {
            clist.get(ie).cur._v = v_left;
            clist.get(ie).cur._rho0 = rho0_left;
            clist.get(ie).cur._eps = eps_left;
        }
        for (int ie = ne_mid; ie < clist.size(); ++ie) {
            clist.get(ie).cur._v = v_right;
            clist.get(ie).cur._rho0 = rho0_right;
            clist.get(ie).cur._eps = eps_right;
        }

        makeP(eqs);
        makeConserved(eqs);
        makeA();

        return 0;
    }

    /**
     * =================================================================== Make
     * the smoothed shock tube.
     * ===================================================================
     */
    /**
     * =================================================================== Make
     * initial data for a smooth jump
     * ===================================================================
     */
    public int makeSmoothJump(eqs eqs, IO io) {
        /**
         * Set default jump parameters and then get the rest.
         */
        x_left = io.lx;
        x_right = io.rx;

        v_left = io.lv;
        v_right = io.rv;

        rho0_left = io.lrho0;
        rho0_right = io.rrho0;

        /**
         * 17 March, 2012: Because I edited PJM's data file, param.dat, to have
         * the left and right pressures and densities, now I have to put the
         * calculation for energy here, so that I don't have to change anything
         * else.
         */
        eps_left = io.lp / ((eqs.gamma - 1.0) * rho0_left);
        eps_right = io.rp / ((eqs.gamma - 1.0) * rho0_right);
        /**
         * Make the grid and geometry.
         */
        int ierr = makeGeometry(io.resolution, x_left, x_right);
        if (ierr != 0) {
            System.err.println("grid::makeSmoothJump: ERROR: makeGeometry returns " + ierr + "\n");
            System.exit(1);
        }

        /**
         * Make the data
         */
        makeSmoothJump(eqs, v_left, rho0_left, eps_left,
                v_right, rho0_right, eps_right);

        /**
         * Make the conserved quantities
         */
        makeConserved(eqs);

        return 0;
    }

    /**
     * =================================================================== Make
     * smooth jump Note that only the elements have data. The nodes compute data
     * as required from the adjoining elements.
     * ===================================================================
     */
    public int makeSmoothJump(eqs eqs, double v_left, double rho0_left, double eps_left,
            double v_right, double rho0_right, double eps_right) {
        /**
         * should check parameters
         */
        if ((rho0_left < 0.0)
                || (rho0_right < 0.0)
                || (eps_left < 0.0)
                || (eps_right < 0.0)) {
            System.err.println("grid::makeSmoothJump: ERROR: incorrect initial data\n" + " rho0_left = " + rho0_left + " rho0_right = " + rho0_right + "\n" + "  eps_left  = " + eps_left + " eps_right  = " + eps_right + "\n");
            System.exit(1);
        }

        /**
         * Set the jump data
         */
        int ne_mid = clist.size() / 2;
        for (int ie = 0; ie < ne_mid; ++ie) {
            clist.get(ie).cur._v = v_left;
            clist.get(ie).cur._rho0 = rho0_left;
            clist.get(ie).cur._eps = eps_left;
        }
        for (int ie = ne_mid; ie < clist.size(); ++ie) {
            clist.get(ie).cur._v = v_right;
            clist.get(ie).cur._rho0 = rho0_right;
            clist.get(ie).cur._eps = eps_right;
        }

        /**
         * smooth it a few times
         */
        final int N_SMOOTH = clist.size() / 4;

        for (int ismooth = 1; ismooth <= N_SMOOTH; ++ismooth) {
            for (int ie = 1; ie < (clist.size() - 1); ++ie) {
                clist.get(ie).cur._v = 0.25
                        * (clist.get(ie - 1).cur._v + 2.0 * clist.get(ie).cur._v + clist.get(ie + 1).cur._v);
                clist.get(ie).cur._rho0 = 0.25
                        * (clist.get(ie - 1).cur._rho0 + 2.0 * clist.get(ie).cur._rho0 + clist.get(ie + 1).cur._rho0);
                clist.get(ie).cur._eps = 0.25
                        * (clist.get(ie - 1).cur._eps + 2.0 * clist.get(ie).cur._eps + clist.get(ie + 1).cur._eps);
            }
        }

        makeP(eqs);
        makeConserved(eqs);
        makeA();

        return 0;
    }

    /**
     * ===================================================================
     */
    /**
     * ===================================================================
     * Methods of the nodeList class from version 2.0.0.
     * ===================================================================
     */
    public void maken(final int n_nodes) {
        try {
            if (n_nodes < 4) {
                System.err.println("nodeList: ERROR: n_nodes = " + n_nodes);
                System.exit(1);
            }
            nlist.clear();
            for (int i = 0; i < n_nodes; i++) {
                node nnode = new node();
                nnode.id = i;
                nlist.add(i, nnode);
            }
        } catch (IndexOutOfBoundsException ioobe) {
            System.err.println("grid :: maken() :: IndexOutOfBoundsException :: ioobe\n" + ioobe);
        }
    }

    public void makeFluxPredictor(artVisControl artvis, eqs eqs, timeData ptime) {
        for (int i = 0; i < nlist.size(); ++i) {
            nlist.get(i).makeFluxPredictor(i, artvis, eqs, clist, ptime);
        }
    }

    public void makeFluxCorrector(artVisControl artvis, eqs eqs, timeData ptime) {
        for (int i = 0; i < nlist.size(); ++i) {
            nlist.get(i).makeFluxCorrector(i, artvis, eqs, clist, ptime);
        }
    }

    /**
     * ===================================================== This is the
     * function which makes the Riemann solutions on the nodes in the node list
     * when the Godunov method is selected.
     * =====================================================
     */
    public void makeRiemannFluxPredictor(timeData ptime, artVisControl artvis, eqs eqs, IO io) {
        for (int j = 0; j < nlist.size(); j++) {
            nlist.get(j).makeRiemannFluxPredictor(j, artvis, clist, nlist, eqs, ptime, io);
        }
    }

    public void makeRiemannFluxCorrector(timeData ptime, artVisControl artvis, eqs eqs, IO io) {
        for (int j = 0; j < nlist.size(); j++) {
            nlist.get(j).makeRiemannFluxCorrector(j, artvis, clist, nlist, eqs, ptime, io);
        }
    }

    /**
     * ===================================================================
     */
    /**
     * ===================================================================
     * Methods of the cellList class from version 2.0.0.
     * ===================================================================
     */
    public void makec(final int n_elements) {
        try {
            if (n_elements < 4) {
                System.err.println("elementList: ERROR: n_elements = " + n_elements + "\n");
                System.exit(1);
            }
            clist.clear();
            /**
             * Cells inside the limits of the size of clist.
             */
            for (int i = 0; i < n_elements; i++) {
                cell ncell = new cell();
                ncell._id = i;
                clist.add(i, ncell);
            }
        } catch (IndexOutOfBoundsException ioobe) {
            System.err.println("grid :: makec() :: IndexOutOfBoundsException :: ioobe\n" + ioobe);
        }
    }

    public void makeConserved(eqs eqs) {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).makeConserved(eqs);
        }
    }

    public void makeP(eqs eqs) {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).makeP(eqs);
        }
    }

    public void makeA() {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).makeA();
        }
    }

    public void makePrimitiveOld(eqs eqs) {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).makePrimitiveOld(eqs);
        }
    }

    public void makePrimitiveCur(eqs eqs) {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).makePrimitiveCur(eqs);
        }
    }

    /**
     * ======================================================== Predictor:
     * Evolve the element quantities using the fluxes at the interfaces.
     * ========================================================
     */
    public void predict(neutrinoMethods neutrinoM, timeData ptime, IO io) {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).evolveSimplePredictor(i, neutrinoM, nlist, ptime, io);
        }
    }

    /**
     * ======================================================== Corrector:
     * Evolve the element quantities using the fluxes at the interfaces.
     * ========================================================
     */
    public void correct(neutrinoMethods neutrinoM, timeData ptime, IO io) {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).evolveSimpleCorrector(i, neutrinoM, nlist, ptime, io);
        }
    }

    public void gridExchange() {
        for (int i = 0; i < clist.size(); ++i) {
            clist.get(i).cellExchange();
        }
    }

    public void makeSigma(IO io) {
        for (int i = 0; i < clist.size(); ++i) {
            makeSigma(i, io);
        }
    }

    public void MakeSigmaMinMod(int i) {

        basic bass = new basic();

        final int iLeft = i - 1;
        final int iCenter = i;
        final int iRight = i + 1;

        if (iLeft == -1) {
            clist.get(iCenter).cur.zeroSigma();
            return;
        }
        if (iRight == clist.size()) {
            clist.get(iCenter).cur.zeroSigma();
            return;
        }

        final cell left = clist.get(iLeft);
        final cell center = clist.get(iCenter);
        final cell right = clist.get(iRight);

        final double dx = right._x - left._x;
        final double dxleft = center._x - left._x;
        final double dxright = right._x - center._x;

        double dS_mid = (right.cur._S - left.cur._S) / dx;
        double dS_left = (center.cur._S - left.cur._S) / dxleft;
        double dS_right = (right.cur._S - center.cur._S) / dxright;

        double dS_mag = bass.min(Math.abs(dS_left), Math.abs(dS_mid), Math.abs(dS_right));
        center.cur._sigmaS = dS_mag * bass.signum(dS_mid);
        if (dS_left * dS_right <= 0.0) {
            center.cur._sigmaS = 0.0;
        }

        double dD_mid = (right.cur._D - left.cur._D) / dx;
        double dD_left = (center.cur._D - left.cur._D) / dxleft;
        double dD_right = (right.cur._D - center.cur._D) / dxright;

        double dD_mag = bass.min(Math.abs(dD_left), Math.abs(dD_mid), Math.abs(dD_right));
        center.cur._sigmaD = dD_mag * bass.signum(dD_mid);
        if (dD_left * dD_right <= 0.0) {
            center.cur._sigmaD = 0.0;
        }

        double dE_mid = (right.cur._E - left.cur._E) / dx;
        double dE_left = (center.cur._E - left.cur._E) / dxleft;
        double dE_right = (right.cur._E - center.cur._E) / dxright;

        double dE_mag = bass.min(Math.abs(dE_left), Math.abs(dE_mid), Math.abs(dE_right));
        center.cur._sigmaE = dE_mag * bass.signum(dE_mid);
        if (dE_left * dE_right <= 0.0) {
            center.cur._sigmaE = 0.0;
        }

        double drho0_mid = (right.cur._rho0 - left.cur._rho0) / dx;
        double drho0_left = (center.cur._rho0 - left.cur._rho0) / dxleft;
        double drho0_right = (right.cur._rho0 - center.cur._rho0) / dxright;

        double drho0_mag = bass.min(Math.abs(drho0_left), Math.abs(drho0_mid), Math.abs(drho0_right));
        center.cur._sigmarho0 = drho0_mag * bass.signum(drho0_mid);
        if (drho0_left * drho0_right <= 0.0) {
            center.cur._sigmarho0 = 0.0;
        }

        double dp_mid = (right.cur._p - left.cur._p) / dx;
        double dp_left = (center.cur._p - left.cur._p) / dxleft;
        double dp_right = (right.cur._p - center.cur._p) / dxright;

        double dp_mag = bass.min(Math.abs(dp_left), Math.abs(dp_mid), Math.abs(dp_right));
        center.cur._sigmap = dp_mag * bass.signum(dp_mid);
        if (dp_left * dp_right <= 0.0) {
            center.cur._sigmap = 0.0;
        }

        double deps_mid = (right.cur._eps - left.cur._eps) / dx;
        double deps_left = (center.cur._eps - left.cur._eps) / dxleft;
        double deps_right = (right.cur._eps - center.cur._eps) / dxright;

        double deps_mag = bass.min(Math.abs(deps_left), Math.abs(deps_mid), Math.abs(deps_right));
        center.cur._sigmaeps = deps_mag * bass.signum(deps_mid);
        if (deps_left * deps_right <= 0.0) {
            center.cur._sigmaeps = 0.0;
        }

        double dv_mid = (right.cur._v - left.cur._v) / dx;
        double dv_left = (center.cur._v - left.cur._v) / dxleft;
        double dv_right = (right.cur._v - center.cur._v) / dxright;

        double dv_mag = bass.min(Math.abs(dv_left), Math.abs(dv_mid), Math.abs(dv_right));
        center.cur._sigmav = dv_mag * bass.signum(dv_mid);
        if (dv_left * dv_right <= 0.0) {
            center.cur._sigmav = 0.0;
        }
    }

    public void MakeSigmaMC(int i) {

        basic bass = new basic();

        final int iLeft = i - 1;
        final int iCenter = i;
        final int iRight = i + 1;

        if (iLeft == -1) {
            clist.get(iCenter).cur.zeroSigma();
            return;
        }
        if (iRight == clist.size()) {
            clist.get(iCenter).cur.zeroSigma();
            return;
        }

        final cell left = clist.get(iLeft);
        final cell center = clist.get(iCenter);
        final cell right = clist.get(iRight);

        final double dx = right._x - left._x;
        final double dxleft = center._x - left._x;
        final double dxright = right._x - center._x;

        double dS_mid = (right.cur._S - left.cur._S) / dx;
        double dS_left = (center.cur._S - left.cur._S) / dxleft;
        double dS_right = (right.cur._S - center.cur._S) / dxright;

        double dS_mag = bass.min(Math.abs(2.0 * dS_left), Math.abs(dS_mid), Math.abs(2.0 * dS_right));
        center.cur._sigmaS = dS_mag * bass.signum(dS_mid);
        if (dS_left * dS_right <= 0.0) {
            center.cur._sigmaS = 0.0;
        }

        double dD_mid = (right.cur._D - left.cur._D) / dx;
        double dD_left = (center.cur._D - left.cur._D) / dxleft;
        double dD_right = (right.cur._D - center.cur._D) / dxright;

        double dD_mag = bass.min(Math.abs(2.0 * dD_left), Math.abs(dD_mid), Math.abs(2.0 * dD_right));
        center.cur._sigmaD = dD_mag * bass.signum(dD_mid);
        if (dD_left * dD_right <= 0.0) {
            center.cur._sigmaD = 0.0;
        }

        double dE_mid = (right.cur._E - left.cur._E) / dx;
        double dE_left = (center.cur._E - left.cur._E) / dxleft;
        double dE_right = (right.cur._E - center.cur._E) / dxright;

        double dE_mag = bass.min(Math.abs(2.0 * dE_left), Math.abs(dE_mid), Math.abs(2.0 * dE_right));
        center.cur._sigmaE = dE_mag * bass.signum(dE_mid);
        if (dE_left * dE_right <= 0.0) {
            center.cur._sigmaE = 0.0;
        }

        double drho0_mid = (right.cur._rho0 - left.cur._rho0) / dx;
        double drho0_left = (center.cur._rho0 - left.cur._rho0) / dxleft;
        double drho0_right = (right.cur._rho0 - center.cur._rho0) / dxright;

        double drho0_mag = bass.min(Math.abs(2.0 * drho0_left), Math.abs(drho0_mid), Math.abs(2.0 * drho0_right));
        center.cur._sigmarho0 = drho0_mag * bass.signum(drho0_mid);
        if (drho0_left * drho0_right <= 0.0) {
            center.cur._sigmarho0 = 0.0;
        }

        double dp_mid = (right.cur._p - left.cur._p) / dx;
        double dp_left = (center.cur._p - left.cur._p) / dxleft;
        double dp_right = (right.cur._p - center.cur._p) / dxright;

        double dp_mag = bass.min(Math.abs(2.0 * dp_left), Math.abs(dp_mid), Math.abs(2.0 * dp_right));
        center.cur._sigmap = dp_mag * bass.signum(dp_mid);
        if (dp_left * dp_right <= 0.0) {
            center.cur._sigmap = 0.0;
        }

        double deps_mid = (right.cur._eps - left.cur._eps) / dx;
        double deps_left = (center.cur._eps - left.cur._eps) / dxleft;
        double deps_right = (right.cur._eps - center.cur._eps) / dxright;

        double deps_mag = bass.min(Math.abs(2.0 * deps_left), Math.abs(deps_mid), Math.abs(2.0 * deps_right));
        center.cur._sigmaeps = deps_mag * bass.signum(deps_mid);
        if (deps_left * deps_right <= 0.0) {
            center.cur._sigmaeps = 0.0;
        }

        double dv_mid = (right.cur._v - left.cur._v) / dx;
        double dv_left = (center.cur._v - left.cur._v) / dxleft;
        double dv_right = (right.cur._v - center.cur._v) / dxright;

        double dv_mag = bass.min(Math.abs(2.0 * dv_left), Math.abs(dv_mid), Math.abs(2.0 * dv_right));
        center.cur._sigmav = dv_mag * bass.signum(dv_mid);
        if (dv_left * dv_right <= 0.0) {
            center.cur._sigmav = 0.0;
        }
    }

    public void makeSigma(int i, IO io) {
        if (io.highResType == 0) {
            clist.get(i).cur.zeroSigma();
            return;
        } else if (io.highResType == 1) {
            MakeSigmaMinMod(i);
            return;
        } else if (io.highResType == 2) {
            MakeSigmaMC(i);
            return;
        } else {
            System.err.println("Did not understand your choice of high resolution method.\n"
                    + "Please choose one of:\n"
                    + "0 = None\n"
                    + "1 = MinMod\n"
                    + "2 = MC\n"
                    + "Thank you.");
            System.exit(1);
        }
    }

    /**
     * =================================================================== Make
     * the shock tube which uses Kuroda et al.'s data for neutrino radiation.
     * ===================================================================
     */
    public int makeKKTShockTube(eqs eqs, IO io, geometrize geo) {
        geo.makeGeometric(eqs, io);

        x_left = geo.geoxl;
        x_right = geo.geoxr;

        v_left = geo.geovl;
        v_right = geo.geovr;

        rho0_left = geo.georho0l;
        rho0_right = geo.georho0r;

        eps_left = geo.geoepsl;
        eps_right = geo.geoepsr;

        /**
         * Make the grid and geometry.
         */
        int ierr = makeGeometry(io.resolution, x_left, x_right);
        if (ierr != 0) {
            System.err.println("grid::makeShockTube: ERROR: makeGeometry returns " + ierr + "\n");
            System.exit(1);
        }

        /**
         * Make the primitive variables
         */
        makeKKTShockTube(eqs, io, x_left, v_left, rho0_left, eps_left,
                x_right, v_right, rho0_right, eps_right);
        /**
         * Make the conserved quantities
         */
        makeConserved(eqs);

        return 0;
    }

    public void makeKKTShockTube(eqs eqs, IO io, double x_left, double v_left, double rho0_left, double eps_left,
            double x_right, double v_right, double rho0_right, double eps_right) {
        /**
         * should check parameters
         */
        if ((rho0_left < 0.0)
                || (rho0_right < 0.0)
                || (eps_left < 0.0)
                || (eps_right < 0.0)) {
            System.err.println("Grid :: makeKKTShockTube: ERROR: Incorrect initial data.\n"
                    + "rho0_left = " + rho0_left + " rho0_right = " + rho0_right + "\n"
                    + "eps_left  = " + eps_left + " eps_right  = " + eps_right);
            System.exit(1);
        }

        if (io.vprofile == 2) {
            /**
             * Make the density initial data, which is a negative ramp.
             */
            
            for (int cells = 0; cells < clist.size(); cells++) {
                clist.get(cells).cur._rho0 = 1.0e14 * (rho0_left * ((x_right - clist.get(cells)._x) / (x_right - x_left))
                        + rho0_right * ((clist.get(cells)._x - x_left) / (x_right - x_left)));

                if (clist.get(cells).cur._rho0 <= 0.0) {
                    System.err.println("WARNING :: grid :: makeKKTShockTube :\n"
                            + "clist.get(cells).cur._rho0 = " + clist.get(cells).cur._rho0
                            + " at cell " + cells + " where: \n"
                            + "rho0_left = " + rho0_left + " :: rho0_right = " + rho0_right
                            + "\nx_left = " + x_left + " :: x_right = " + x_right
                            + "\nclist.get(cells)._x = " + clist.get(cells)._x + "\n"
                            + "Making positive ...");
                    System.exit(1);
                }
            }

            /**
             * Make the energy initial data, which is a negative ramp.
             */
            
            for (int ic = 0; ic < clist.size(); ic++) {
                clist.get(ic).cur._eps = 1.0e15 * (eps_left * ((x_right - clist.get(ic)._x) / (x_right - x_left))
                        + eps_right * ((clist.get(ic)._x - x_left) / (x_right - x_left)));

                if (clist.get(ic).cur._eps <= 0.0) {
                    System.err.println("ERROR :: grid :: makeKKTShockTube :\n"
                            + "clist.get(ic).cur._eps = " + clist.get(ic).cur._eps
                            + " at cell " + ic + " where: \n"
                            + "eps_left = " + eps_left + " :: eps_right = " + eps_right
                            + "\nx_left = " + x_left + " :: x_right = " + x_right
                            + "\nclist.get(ic)._x = " + clist.get(ic)._x + "\n"
                            + "Making positive ...");
                    clist.get(ic).cur._eps = 0.0;
                }
            }

            /**
             * Make the velocity initial data, which is v = v_{0r} *
             * ((r/x_{r})^{2} - 2.0 * (r/x_{r}) + 1). The left velocity and the
             * right velocity are 0 ms^{-1}. This cannot (obviously) be used.
             */
            int mid = clist.size() / 2;
            double xmid = (x_right + x_left) / 2.0;

            for (int icel = 0; icel < mid; icel++) {
                clist.get(icel).cur._v = v_left * (Math.pow(((x_left - clist.get(icel)._x) / (x_left - xmid)), 2.0));
            }

            for (int icel = mid; icel < clist.size(); icel++) {
                clist.get(icel).cur._v = v_right * (Math.pow(((clist.get(icel)._x - x_right) / (xmid - x_right)), 2.0));
            }
        } else if (io.vprofile == 1) {
            int ne_mid = clist.size() / 2;
            for (int ie = 0; ie < ne_mid; ++ie) {
                clist.get(ie).cur._v = v_left;
                clist.get(ie).cur._rho0 = rho0_left;
                clist.get(ie).cur._eps = eps_left;
            }
            for (int ie = ne_mid; ie < clist.size(); ++ie) {
                clist.get(ie).cur._v = v_right;
                clist.get(ie).cur._rho0 = rho0_right;
                clist.get(ie).cur._eps = eps_right;
            }
        }

        makeP(eqs);
        makeConserved(eqs);
        makeA();
    }
}