/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author gregory
 */
public final class eqs {

    double gamma;      /// Polytropic index
    double gam1;       /// gamma-1 (set in constructors)
    double kappa;
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

    public eqs() {

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
        massCore = 0.0;
        try {
            readParameters();
        } catch (IOException ex) {
            System.err.println("ERROR :: In eqs, readParameters cannot read the param.dat file.");
            Logger.getLogger(IO.class.getName()).log(Level.SEVERE, null, ex);
        }
        Make(gamma);
    }

    public final void Make(double gamma_in) {
        if (gamma_in < 1.0) {
            System.err.println("eqs :: ERROR :: gamma = " + gamma + "\n");
            System.exit(1);
        }
        gamma = gamma_in;
        gam1 = gamma - 1.0;
    }

    public double P(double rho0, double eps) {
        return (gam1 * rho0 * eps);
    }

    public double EPS(double p, double rho0) {
        return (p / (gam1 * rho0));
    }

    public double CS(double eps) {
        return (Math.sqrt(gamma * gam1 * eps / (1.0 + gamma * eps)));
    }

    /**
     * ====================================================================
     * readParameters reads the contents of "param.dat" as input to the Godunov
     * code, NOT the Riemann solver.
     * ====================================================================
     */
    public void readParameters() throws IOException {
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
        } catch (IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
