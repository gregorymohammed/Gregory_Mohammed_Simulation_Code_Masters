/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

/**
 *
 * @author gregory
 */
public class evolutionVariables {

    double S, D, E;

    public evolutionVariables() {
        S = 0.0;
        D = 0.0;
        E = 0.0;
    }

    public evolutionVariables(stateVariables st) {
    }

    public evolutionVariables(final double invalue) {
        S = invalue;
        D = invalue;
        E = invalue;
    }

    public evolutionVariables(final evolutionVariables ev) {
        S = ev.S;
        D = ev.D;
        E = ev.E;
    }

    public void setSelf(final evolutionVariables ev) {
        S = ev.S;
        D = ev.D;
        E = ev.E;
    }

    public void setToValue(final double invalue) {
        S = invalue;
        D = invalue;
        E = invalue;
    }

    public void addToSelf(final evolutionVariables ev) {
        S += ev.S;
        D += ev.D;
        E += ev.E;
    }

    public void addScalar(final double scalar) {
        S += scalar;
        D += scalar;
        E += scalar;
    }

    public void multiplyByScalar(final double scalar) {
        S *= scalar;
        D *= scalar;
        E *= scalar;
    }
}
