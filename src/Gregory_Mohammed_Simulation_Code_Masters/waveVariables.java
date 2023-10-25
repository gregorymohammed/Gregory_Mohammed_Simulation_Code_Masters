/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Gregory_Mohammed_Simulation_Code_Masters;

/**
 * @author Gregory Mohammed
 */
public class waveVariables {

    public double rho0;
    public double u;
    public double h;
    public double cs;
    public double vel;
    public double vshock;

    public waveVariables() {
        rho0 = 0.0;
        u = 0.0;
        h = 0.0;
        cs = 0.0;
        vel = 0.0;
        vshock = 0.0;
    }

    public void MakeU(eqs eqs, double p) {
        eqs.EPS(p, rho0);
    }
    
    public void Print() {
        System.out.print( "waveVariables: .rho0=" + rho0 + " .u=" + u + " .h=" + h
                + " .cs=" + cs + " .vel=" + vel + " .vshock=" + vshock);               
    }
}
