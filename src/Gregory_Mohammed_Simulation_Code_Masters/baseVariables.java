package Gregory_Mohammed_Simulation_Code_Masters;

public class baseVariables {

    /**
     * ============================================== 
     * This class contains the
     * variables which are fundamental to the node and cell classes.
     * ==============================================
     */
    public double _rho0;
    public double _p;
    public double _eps;
    public double _v;
    public double _a;
    public double _D;
    public double _E;
    public double _S;
    public double _sigmaD;
    public double _sigmaE;
    public double _sigmaS;
    public double _sigmarho0;
    public double _sigmap;
    public double _sigmaeps;
    public double _sigmav;

    /**
     * ============================================== 
     * Constructor.
     * ==============================================
     */
    public baseVariables() {
        _rho0 = 0.0;
        _p = 0.0;
        _eps = 0.0;
        _v = 0.0;
        _a = 0.0;
        _D = 0.0;
        _E = 0.0;
        _S = 0.0;
    }

    /**
     * ============================================== 
     * Methods.
     * ==============================================
     */
    public int makePrimitive(eqs eqs, primitiveSolve psolve) {
        /**
         * ========================================== 
         * Test for sensible values.
         * Exit if they are bad. 
         * ==========================================
         */
        int iError = psolve.primitiveSolveFunction(eqs, _S, _D, _E, _p, _v, _rho0, _eps, _p);

        if (iError != 0) {
            System.err.println("baseVariables :: ERROR :: makePrimitive failed to solve for primitive variables.\n");
            System.err.println("_v = " + _v + "\n_rho0 = " + _rho0 + "\n_eps = " + _eps + "\n_p = " + _p + "\nS = " + _S + "\nD = " + _D + "\nE = " + _E + "\n");
            return 1;
        }
        /**
         * ========================================== 
         * If they are good, set a.
         * ==========================================
         */
        _a = (1.0 / Math.sqrt(1.0 - _v * _v));

        copyPSolveIn(psolve);

        return 0;
    }

    public void copyPSolveIn(primitiveSolve ps) {
        _rho0 = ps._rho0;
        _p = ps._p;
        _eps = ps._eps;
        _v = ps._v;
    }

    public void zeroSigma() {
        //System.out.println("Zero'ing sigma ...");

        _sigmaD = 0.0;
        _sigmaE = 0.0;
        _sigmaS = 0.0;
        _sigmarho0 = 0.0;
        _sigmap = 0.0;
        _sigmaeps = 0.0;
        _sigmav = 0.0;
    }

    //============================================================================
    /** Make constant extrapolations from base position to highres position.
    
    \param[in] elem Element from which to start
    \param[in] x Position at which to calculate the high-resolution extrapolation
    \param[out] highres High resolution base variables
     */
    public void MakeHighResConstant(cell cell, double x) {
        _S = cell.cur._S;
        _D = cell.cur._D;
        _E = cell.cur._E;
        _rho0 = cell.cur._rho0;
        _p = cell.cur._p;
        _v = cell.cur._v;
    }

    //============================================================================
    /** Make linear extrapolations from base position to highres position.
    
    \param[in] elem Element from which to start
    \param[in] x Position at which to calculate the high-resolution extrapolation
    \param[out] highres High resolution base variables
     */
    public void MakeHighResLinear(cell cell, double x) {
        final double x_base = cell._x;

        final double dx = x - x_base;

        _S = cell.cur._S + cell.cur._sigmaS * dx;
        _D = cell.cur._D + cell.cur._sigmaD * dx;
        _E = cell.cur._E + cell.cur._sigmaE * dx;
        _rho0 = cell.cur._rho0 + cell.cur._sigmarho0 * dx;
        _p = cell.cur._p + cell.cur._sigmap * dx;
        _v = cell.cur._v + cell.cur._sigmav * dx;
    }

//=============================================================================
    /** Cubic Hermite interpolation basis functions
    
    \param[in] nim Linear basis function N_{i-1} on the left
    \param[in] ni  Linear basis function N_{i} on the right
    \param[in] h length left to right edge
    \param[out] cubic  Cubic Hermite basis functions
     */
    public void CubicHermite(double nim, double ni, double h, double[] cubic) {
        cubic[0] = nim * nim * (3.0 - nim * 2.0);
        cubic[2] = 1.0 - cubic[0];

        final double a = h * nim * ni;

        cubic[3] = -a * ni;
        cubic[1] = a * nim;
    }

//============================================================================
    /** Make Cubic Hermite extrapolations from the elements elem1 and elem2 to x
    
    Note: elem1 must be to the left of elem2
    
    \param[in] elem1 Element 1
    \param[in] elem2 Element 2
    \param[in] x  Position at which to calculate the high-resolution (cubic hermite) extrapolation
    \param[out] highres_left High resolution base variables extrapolated to x
     */
    public void MakeHighResCubicHermite(cell cell1, cell cell2, double x) {
        if (cell1._x >= cell2._x) {
            System.err.append("BaseVariables :: MakeHighResCubicHermite: ERROR: cell1 and cell2 are not in order.");
            System.exit(1);
        }

        final double h = cell2._x - cell2._x;
        final double ni = (x - cell1._x) / h;
        //const double ni = (elem2->x - x) / h;

        final double nim = 1.0 - ni;

        double[] cubic = new double[4];

        CubicHermite(nim, ni, h, cubic);

        _S = cell1.cur._S * cubic[0] + cell1.cur._sigmaS * cubic[1]
                + cell2.cur._S * cubic[2] + cell2.cur._sigmaS * cubic[3];

        _D = cell1.cur._D * cubic[0] + cell1.cur._sigmaD * cubic[1]
                + cell2.cur._D * cubic[2] + cell2.cur._sigmaD * cubic[3];

        _E = cell1.cur._E * cubic[0] + cell1.cur._sigmaE * cubic[1]
                + cell2.cur._E * cubic[2] + cell2.cur._sigmaE * cubic[3];

        _p = cell1.cur._p * cubic[0] + cell1.cur._sigmap * cubic[1]
                + cell2.cur._p * cubic[2] + cell2.cur._sigmap * cubic[3];

        _rho0 = cell1.cur._rho0 * cubic[0] + cell1.cur._sigmarho0 * cubic[1]
                + cell2.cur._rho0 * cubic[2] + cell2.cur._sigmarho0 * cubic[3];

        _v = cell1.cur._v * cubic[0] + cell1.cur._sigmav * cubic[1]
                + cell2.cur._v * cubic[2] + cell2.cur._sigmav * cubic[3];
    }
}
