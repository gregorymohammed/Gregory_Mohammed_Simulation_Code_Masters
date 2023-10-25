package Gregory_Mohammed_Simulation_Code_Masters;

import java.util.*;

/**
 *
 * @author gregory
 */
public class nodalRiemannSolver {

    public double x_left;
    public double x_right;
    public int n;
    public double T;
    public int goodSolve;

    public nodalRiemannSolver() {
        x_left = 0.0;
        x_right = 0.0;
        n = 100;
        T = 1.0;
        goodSolve = 0;
    }

    public void makeTheRiemannSolution(int j, eqs eqs, ArrayList<cell> clist, exactNode node,
            stateVariables lstate, stateVariables rstate, timeData ptime, riemannExact rexact) throws Exception {

        if (j == 0) {
            return;
        }

        if (j == clist.size()) {    //PJM: j=node number, clist is cell list, which is 1 less than nlist
            return;
        }
        /**
         * =====================================================
         * REMEMBER: The position, X(), of the NODE, is NOT the discontinuity!
         * Define the size of the Riemann shock tube.
         * =====================================================
         */
        /**PJM:  Not correct, the discontinuity is the position of the node.
         * However the position should not be significant, so it would be easiest just to set = 0
        x_left = clist.get(j - 1)._x;
        x_right = clist.get(j)._x;
        
        double x_discontinuity = 0.5 * (x_left + x_right);
         */
        double x_discontinuity = 0.0;    // Local solution at x=0.0.  Easiest approach.

        T = ptime.delta;

        rexact.Solve(T, x_discontinuity, eqs, lstate, rstate);

        if (rexact.T <= 0.0) {
            System.out.println("nodalRiemannSolver: ERROR: rexact.T=" + rexact.T
                    + " <=0. Means rexact has only default constructed values");
            System.out.print("  rexact: ");
            rexact.Print();
            System.out.println();
            System.out.print("  lstate: ");
            lstate.Print();
            System.out.println();
            System.out.print("  rstate: ");
            rstate.Print();
            System.out.println();
            System.exit(-1);
        }

        // Note: the "node" (an exactNode) is computed by GetExactValues 

        rexact.GetExactValues(x_discontinuity, eqs, node, lstate, rstate);
    }
}
