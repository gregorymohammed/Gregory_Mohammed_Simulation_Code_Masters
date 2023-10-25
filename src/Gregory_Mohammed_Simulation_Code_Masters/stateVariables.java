/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

/**
 * @author Gregory Mohammed
 */
public final class stateVariables {

    /** rest density (PJM: rho_0) */
    public double rho0;
    /** pressure */
    public double p;
    /** specific internal energy (PJM: epsilon) */
    public double u;
    /** specific enthalpy */
    public double h;
    /** sound speed */
    public double cs;
    /** flow velocity */
    public double vel;
    /** Lorentz factor 1/sqrt(1-v^2) */
    public double w;
    /** Cell from which this state was corrected.  Used only for debugging output. */
    public int icell;
    
    public stateVariables() {
        rho0 = 0.0;
        p = 0.0;
        u = 0.0;
        h = 0.0;
        cs = 0.0;
        vel = 0.0;
        w = 0.0;
        icell = -10;
    }

    public stateVariables(final stateVariables st) {
        rho0 = st.rho0;
        p = st.p;
        u = st.u;
        h = st.h;
        cs = st.cs;
        vel = st.vel;
        w = st.w;
        icell = st.icell;
    }
    
    public stateVariables(eqs eqs, double rho0_in, double p_in, double v_in) {
        rho0 = rho0_in;
        p = p_in;
        vel = v_in;
        u = eqs.EPS(p, rho0);
        cs = eqs.CS(u);
        MakeH();
        MakeW();
    }

    public void MakeU(eqs eqs) {
        u = eqs.EPS(p, rho0);
    }

    public void MakeH() {
        h = 1.0 + u + p / rho0;
    }

    public void MakeCS(eqs eqs) {
        cs = eqs.CS(u);
    }

    public void MakeW() {
        w = 1.0 / Math.sqrt(1.0 - vel * vel);
    }
    
    public void MakeP(eqs eqs) {
        p = (eqs.gamma - 1.0) * rho0 * u;
    }
    
    public void Print() {
        System.out.println("OUTPUT :: stateVariables :: \n"
                + "icell = " + icell + " :: rho0 = " + rho0 + "\n"
                + "p = " + p + " :: u = " + u + "\n"
                + "h = " + h + " :: cs = " + cs + "\n"
                + "vel = " + vel + " :: w = " + w);
    }
}
