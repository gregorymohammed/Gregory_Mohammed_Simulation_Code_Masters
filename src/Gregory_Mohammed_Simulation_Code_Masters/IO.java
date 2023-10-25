/*
 * This class handles file reading and writing.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Gregory Mohammed
 */
public final class IO {

    /**
     * ===================================== These 2 variables determine:
     * initShock: 0 = smooth initial shock 1 = real shock method: 0 = Upwind 1 =
     * Godunov resolution = number of cells The values are set in the param.dat
     * file. =====================================
     */
    greg_const gconst;
    int initShock;
    int showInitialSlice;
    int showExact;
    int method;
    int neutrinos;
    double fluxTuner;
    double energyTuner;
    int highResType;
    int reconstructionType;
    int resolution;
    /**
     * =====================================
     */
    int dataSwitch;
    int plotGeo;
    double xLeft;
    double xRight;
    double vLeft;
    double vRight;
    double lrho0nu;
    double rrho0nu;
    double lpnu;
    double rpnu;
    /**
     * =====================================
     */
    double gamma;
    int totalSteps;
    double totalTimeSod;
    double totalTimeKKT;
    double lp;
    double lrho0;
    double lv;
    double rp;
    double rrho0;
    double rv;
    double lx;
    double rx;
    double epsdmp;
    int dmpinterval;
    int corrector;
    double artvisk1;
    double artvisk2;
    double k;
    double courant;
    double p1;
    int vprofile;
    double massCore;
    stateVariables forRiemannLeft;
    stateVariables forRiemannRight;

    public IO(eqs eqs, stateVariables lstate, stateVariables rstate) {
        gconst = new greg_const();
        forRiemannLeft = new stateVariables();
        forRiemannRight = new stateVariables();

        showInitialSlice = 0;
        showExact = 0;
        initShock = 0;
        method = 0;
        neutrinos = 0;
        fluxTuner = 0.0;
        energyTuner = 0.0;
        highResType = 0;
        reconstructionType = 0;
        resolution = 100;

        dataSwitch = 0;
        plotGeo = 0;
        xLeft = 0.0;
        xRight = 0.0;
        vLeft = 0.0;
        vRight = 0.0;
        lrho0nu = 0.0;
        rrho0nu = 0.0;
        lpnu = 0.0;
        rpnu = 0.0;

        gamma = 0.0;
        totalSteps = 0;
        totalTimeSod = 0.0;
        totalTimeKKT = 0.0;
        lp = 0.0;
        lrho0 = 0.0;
        lv = 0.0;
        rp = 0.0;
        rrho0 = 0.0;
        rv = 0.0;
        lx = 0.0;
        rx = 0.0;

        epsdmp = 0.0;
        dmpinterval = 0;
        corrector = 0;
        artvisk1 = 0.0;
        artvisk2 = 0.0;
        k = 0.0;
        courant = 0.0;
        p1 = 0.0;
        vprofile = 0;
        massCore = 0;
        try {
            readData(eqs, lstate, rstate);
        } catch (IOException ex) {
            System.err.println("ERROR :: In IO, readData cannot read the param.dat file.");
            Logger.getLogger(IO.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * ====================================================================
     * readData reads the data from the datafile "sodValues.dat".
     * ====================================================================
     */
    public void readData(eqs eqs, stateVariables lstate, stateVariables rstate) throws IOException {
        try {
            FileInputStream fileinput = new FileInputStream("param.dat");
            Scanner input = new Scanner(fileinput);

            while (input.hasNextLine()) {

                String line = input.nextLine();
                if (!line.startsWith("#")) {
                    if (line.startsWith("initShock = ")) {
                        initShock = Integer.parseInt(line.substring(12));
                    } else if (line.startsWith("showInitialSlice = ")) {
                        showInitialSlice = Integer.parseInt(line.substring(19));
                    } else if (line.startsWith("showExact = ")) {
                        showExact = Integer.parseInt(line.substring(12));
                    } else if (line.startsWith("method = ")) {
                        method = Integer.parseInt(line.substring(9));
                    } else if (line.startsWith("neutrinos = ")) {
                        neutrinos = Integer.parseInt(line.substring(12));
                    } else if (line.startsWith("fluxTuner = ")) {
                        fluxTuner = Double.parseDouble(line.substring(12));
                    } else if (line.startsWith("energyTuner = ")) {
                        energyTuner = Double.parseDouble(line.substring(14));
                    } else if (line.startsWith("highResType = ")) {
                        highResType = Integer.parseInt(line.substring(14));
                    } else if (line.startsWith("reconstructionType = ")) {
                        reconstructionType = Integer.parseInt(line.substring(21));
                    } else if (line.startsWith("resolution = ")) {
                        resolution = Integer.parseInt(line.substring(13));
                    } else if (line.startsWith("gamma = ")) {
                        gamma = Double.parseDouble(line.substring(8));
                    } else if (line.startsWith("k = ")) {
                        k = Double.parseDouble(line.substring(4));
                    } else if (line.startsWith("totalSteps = ")) {
                        totalSteps = Integer.parseInt(line.substring(13));
                    } else if (line.startsWith("totalTimeSod = ")) {
                        totalTimeSod = Double.parseDouble(line.substring(15));
                    } else if (line.startsWith("totalTimeKKT = ")) {
                        totalTimeKKT = Double.parseDouble(line.substring(15));
                    } else if (line.startsWith("epsdmp = ")) {
                        epsdmp = Double.parseDouble(line.substring(9));
                    } else if (line.startsWith("dmpinterval = ")) {
                        dmpinterval = Integer.parseInt(line.substring(14));
                    } else if (line.startsWith("corrector = ")) {
                        corrector = Integer.parseInt(line.substring(12));
                    } else if (line.startsWith("artvis.k1 = ")) {
                        artvisk1 = Double.parseDouble(line.substring(12));
                    } else if (line.startsWith("artvis.k2 = ")) {
                        artvisk2 = Double.parseDouble(line.substring(12));
                    } else if (line.startsWith("courant = ")) {
                        courant = Double.parseDouble(line.substring(10));
                    } else if (line.startsWith("p1 = ")) {
                        p1 = Double.parseDouble(line.substring(5));
                    } else if (line.startsWith("dataSwitch = ")) {
                        dataSwitch = Integer.parseInt(line.substring(13));
                    } else if (line.startsWith("plotGeo = ")) {
                        plotGeo = Integer.parseInt(line.substring(10));
                    } else if (line.startsWith("x_left = ")) {
                        lx = Double.parseDouble(line.substring(9));
                    } else if (line.startsWith("x_right = ")) {
                        rx = Double.parseDouble(line.substring(10));
                    } else if (line.startsWith("v_left = ")) {
                        lv = Double.parseDouble(line.substring(9));
                    } else if (line.startsWith("v_right = ")) {
                        rv = Double.parseDouble(line.substring(10));
                    } else if (line.startsWith("rho0_left = ")) {
                        lrho0 = Double.parseDouble(line.substring(12));
                    } else if (line.startsWith("rho0_right = ")) {
                        rrho0 = Double.parseDouble(line.substring(13));
                    } else if (line.startsWith("p_left = ")) {
                        lp = Double.parseDouble(line.substring(9));
                    } else if (line.startsWith("p_right = ")) {
                        rp = Double.parseDouble(line.substring(10));
                    } else if (line.startsWith("xLeft = ")) {
                        xLeft = Double.parseDouble(line.substring(8));
                    } else if (line.startsWith("xRight = ")) {
                        xRight = Double.parseDouble(line.substring(9));
                    } else if (line.startsWith("vLeft = ")) {
                        vLeft = Double.parseDouble(line.substring(8));
                    } else if (line.startsWith("vRight = ")) {
                        vRight = Double.parseDouble(line.substring(9));
                    } else if (line.startsWith("lrho0nu = ")) {
                        lrho0nu = Double.parseDouble(line.substring(10));
                    } else if (line.startsWith("rrho0nu = ")) {
                        rrho0nu = Double.parseDouble(line.substring(10));
                    } else if (line.startsWith("lpnu = ")) {
                        lpnu = Double.parseDouble(line.substring(7));
                    } else if (line.startsWith("rpnu = ")) {
                        rpnu = Double.parseDouble(line.substring(7));
                    } else if (line.startsWith("vprofile = ")) {
                        vprofile = Integer.parseInt(line.substring(11));
                    } else if (line.startsWith("massCore = ")) {
                        massCore = Double.parseDouble(line.substring(11));
                    }
                }
            }
            input.close();

            lstate.p = lp;
            lstate.rho0 = lrho0;
            lstate.vel = lv;

            rstate.p = rp;
            rstate.rho0 = rrho0;
            rstate.vel = rv;

            lstate.MakeU(eqs);
            lstate.MakeH();
            lstate.MakeCS(eqs);
            lstate.MakeW();

            rstate.MakeU(eqs);
            rstate.MakeH();
            rstate.MakeCS(eqs);
            rstate.MakeW();

            /**
             * This chunk is to provide left and right states for the exact
             * Riemann solution at t = 0.25. This is ONLY for the classic Sod
             * shock tube, which is why those values are hard-coded.
             */
            forRiemannLeft.p = 1.0;
            forRiemannLeft.rho0 = 1.0;
            forRiemannLeft.vel = 0.0;

            forRiemannRight.p = 0.1;
            forRiemannRight.rho0 = 0.125;
            forRiemannRight.vel = 0.0;

            forRiemannLeft.MakeU(eqs);
            forRiemannLeft.MakeH();
            forRiemannLeft.MakeCS(eqs);
            forRiemannLeft.MakeW();

            forRiemannRight.MakeU(eqs);
            forRiemannRight.MakeH();
            forRiemannRight.MakeCS(eqs);
            forRiemannRight.MakeW();

        } catch (IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void printIO(stateVariables lstate, stateVariables rstate) {
        /**
         * This prints all of IO's variable values.
         */
        System.err.println("===========================================================================================================\n"
                + "INFO OUTPUT :: IO :: printIO :\n"
                + "===========================================================================================================\n"
                + "showInitialSlice = " + showInitialSlice + " :: initShock = " + initShock + " :: method = " + method + "\n"
                + "neutrinos = " + neutrinos + " :: fluxTuner = " + fluxTuner + " :: energyTuner = " + energyTuner + "\n"
                + "highResType = " + highResType + " :: reconstructionType = " + reconstructionType + "\n"
                + "resolution = " + resolution + " :: dataSwitch = " + dataSwitch + " :: plotGeo = " + plotGeo + "\n"
                + "===========================================================================================================\n"
                + "xLeft = " + xLeft + " :: xRight = " + xRight + " :: vLeft = " + vLeft + " :: vRight = " + vRight + "\n"
                + "lrho0nu = " + lrho0nu + " :: rrho0nu = " + rrho0nu + " :: lpnu = " + lpnu + " :: rpnu = " + rpnu + "\n"
                + "gamma = " + gamma + " :: totalSteps = " + totalSteps + " :: \n"
                + "totalTimeSod = " + totalTimeSod + " :: totalTimeKKT = " + totalTimeKKT + "\n"
                + "lp = " + lp + " :: lrho0 = " + lrho0 + " :: lv = " + lv + "\n"
                + "rp = " + rp + " :: rrho0 = " + rrho0 + " :: rv = " + rv + "\n"
                + "lx = " + lx + " :: rx = " + rx + "\n"
                + "epsdmp = " + epsdmp + " :: dmpinterval = " + dmpinterval + " :: corrector = " + corrector + "\n"
                + "artvisk1 = " + artvisk1 + " :: artvisk2 = " + artvisk2 + " :: courant = " + courant + "\n"
                + "p1 = " + p1 + " :: vprofile = " + vprofile + "\n"
                + "===========================================================================================================\n"
                + "lstate.rho0 = " + lstate.rho0 + " :: rstate.rho0 = " + rstate.rho0 + "\n"
                + "lstate.p = " + lstate.p + " :: rstate.p = " + rstate.p + "\n"
                + "lstate.vel = " + lstate.vel + " :: rstate.vel = " + rstate.vel + "\n"
                + "lstate.u = " + lstate.u + " :: rstate.u = " + rstate.u + "\n"
                + "===========================================================================================================");
    }
}
