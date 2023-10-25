package Gregory_Mohammed_Simulation_Code_Masters;

public class flux {

    /**
     * ==============================================
     * The class which is responsible for anything to
     * do with fluxes.
     * ==============================================
     */
    public double _D;
    public double _E;
    public double _S;

    /**
     * ==============================================
     * Constructor.
     * ==============================================
     */
    public flux() {
        _D = 0;
        _E = 0;
        _S = 0;
    }
    
    void Print(){
        System.out.print( "._D=" + _D + " ._E=" + _E + " ._S=" + _S );
    }

    /**
     * ==============================================
     * Methods.
     * ==============================================
     */
    public void make(baseVariables b, double q) {
        _D = b._D * b._v;
        _E = (b._E * b._v) + ((b._p + q) * b._v);
        _S = (b._S * b._v) + b._p + q;
    }

    public void make(double aS, double aD, double aE, double av, double ap, double q) {
        /**
         * =========================================
         * Bug in PJM's code, which manifested in the Godunov code.
         * Because the Riemann solver only works on the evolution of the shock,
         * S, D, E need only be calculated then. Thus, v needs to be present on
         * all values. This only happens for S = [E + p + D]v (MM, p38).
         * =========================================
         */
        _D = aD * av;
        _E = (aE * av) + ((ap + q) * av);
        _S = (aS * av) + ap + q;
        //_S = (_E + ap + _D) * av;
    }
}
