package Gregory_Mohammed_Simulation_Code_Masters;

/**
 *
 * @author gregory
 */
public class exactNode {

    int N;
    double x;
    double rho0; /// rest density
    double p;
    double u;   /// epsilon in my notation
    double vel;
    double h;   /// enthalpy w in my notation
    double w;   /// w=1/sqrt(1-v^2).  a in my notation
    double S, D, E; /// These are the conserved quantities in my notation

    public exactNode () {}

    public exactNode(final int n) {
        N = n;
    }

    public void MakewhDSE() {
        Makew();
        Makeh();
        MakeD();
        MakeS();
        MakeE();
    } // order is important.
    
    void Print()
    {
        System.out.print( " .N=" + N + " .x=" + x + " .rho0=" + rho0 + " .p=" + p
                + " .u=" + u + " .vel=" + vel + ". h=" + h + " .w=" + w
                + " .S=" + S + " .D=" + D + " .E=" + E );
    }

    public void Makew() {
        w = 1.0 / Math.sqrt(1.0 - vel * vel);
    }

    public void Makeh() {
        h = 1.0 + u + p / rho0;
    }

    public void MakeD() {
        D = rho0 / Math.sqrt(1.0 - vel * vel);
    }

    public void MakeS() {
        S = rho0 * h * w * w * vel;
    }

    public void MakeE() {
        E = rho0 * h * w * w - p - D;
    }  // D must already be made
}
