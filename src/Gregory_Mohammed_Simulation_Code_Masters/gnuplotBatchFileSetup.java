package Gregory_Mohammed_Simulation_Code_Masters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Gregory Mohammed
 *
 * 5 August, 2012
 *
 * This class is used to generate the gnuplot script which the main program will
 * use to produce plots. The class needs to read all the files in the working
 * folder, then extract the .dat files. It will then use those files for the
 * plots.
 */
public final class gnuplotBatchFileSetup {

    /// Extension to look for.
    String plotpattern;
    String csvpattern;
    /// The Riemann exact solution filename.
    String riemann;
    String sd;
    /// Hold i when riemann is a match with listFile[i].
    int riemannLocation;
    /// The Approximate solution filename.
    String approximate;
    /// Hold i when approximate is a match with listFile[i].
    int approximateLocation;
    /// Directory to work in.
    File dir;
    /// Array to hold .dat filenames.
    File listFile[];
    /// Batch file names.
    String initialVel;
    String initialRho0;
    String initialEps;
    String initialP;
    String rho0;
    String p;
    String u;
    String v;

    public gnuplotBatchFileSetup() {
        initialVel = "initialVel.eps";
        initialRho0 = "initialRho0.eps";
        initialEps = "initialEps.eps";
        initialP = "initialP.eps";
    }

    public gnuplotBatchFileSetup(IO io, String approximate_in) throws FileNotFoundException, IOException {
        plotpattern = ".plt";
        csvpattern = ".csv";
        riemann = "Riemann_Exact_Solution.csv";
        sd = "standard_deviations.csv";
        riemannLocation = 0;
        approximate = approximate_in;
        approximateLocation = 0;
        rho0 = "rho0.eps";
        p = "p.eps";
        u = "u.eps";
        v = "v.eps";

        //changeExtension = new String[58];
        dir = new File(System.getProperty("user.dir"));
        listFile = dir.listFiles();

        FileOutputStream writeAnimation = null;
        writeAnimation = new FileOutputStream("animate_pressure.bat");
        PrintStream gnuplotAnimate = new PrintStream(writeAnimation);

        FileOutputStream writePicture = null;
        writePicture = new FileOutputStream("picture_pressure.bat");
        PrintStream gnuplotFrame = new PrintStream(writePicture);

        /**
         * This constructor loads all the .plt files into the listFile array. It
         * also gets the location of the Riemann_Exact_Solution.plt filename for
         * later use.
         */
        if (listFile != null) {
            for (int i = 0; i < listFile.length; i++) {
                if ((listFile[i].getName().endsWith(plotpattern)) || (listFile[i].getName().endsWith(csvpattern))) {
                    if (listFile[i].getName() == null ? riemann == null : listFile[i].getName().equals(riemann)) {
                        riemannLocation = i;
                    } else if (listFile[i].getName() == null ? approximate == null : listFile[i].getName().equals(approximate)) {
                        approximateLocation = i;
                    }
                }
                /**
                 * This bit is to be used for getting a file with a list of all
                 * the approximate solution files. That is then used with ...
                 * for the purposes of generating a .bat file for gnuplot to
                 * animate the data. .pic files are also generated, which are
                 * used by GIMP to generate an animated .gif.
                 */
                if ((listFile[i].getName().endsWith(csvpattern)) != (listFile[i].getName().equals(riemann))) {
                    if ((listFile[i].getName().endsWith(csvpattern)) != (listFile[i].getName().equals(sd))) {
                        approximateLocation = i;

                        if (i < (listFile.length / 2.0)) {
                            gnuplotAnimate.printf("set title \"Approximate Fluid Pressure.\"\n"
                                    + "set yrange [0:1.7e44]\n"
                                    + "set ylabel \"pressure (x10^-3 Pa)\"\n"
                                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                                    + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\"\n"
                                    + "pause 0.00001\n");

                            gnuplotFrame.printf("set title \"Approximate Fluid Pressure.\"\n"
                                    + "set yrange [0:1.7e44]\n"
                                    + "set ylabel \"pressure (x10^-3 Pa)\"\n"
                                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                                    + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\"\n"
                                    + "set output \"" + listFile[approximateLocation].getName() + ".png\"\n"
                                    + "set terminal pngcairo\n"
                                    + "replot\n");
                        } else {
                            gnuplotAnimate.printf("set title \"Approximate Fluid Pressure.\"\n"
                                    + "set yrange [0:7e43]\n"
                                    + "set ylabel \"pressure (x10^-3 Pa)\"\n"
                                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                                    + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\"\n"
                                    + "pause 0.00001\n");

                            gnuplotFrame.printf("set title \"Approximate Fluid Pressure.\"\n"
                                    + "set yrange [0:7e43]\n"
                                    + "set ylabel \"pressure (x10^-3 Pa)\"\n"
                                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                                    + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\"\n"
                                    + "set output \"" + listFile[approximateLocation].getName() + ".png\"\n"
                                    + "set terminal pngcairo\n"
                                    + "replot\n");
                        }
                    }
                }
            }
            writeAnimation.close();
            writePicture.close();
            gnuplotBatchScript(io);
        } else {
            System.err.println("WARNING :: There are no .csv files!");
            System.exit(1);
        }
    }

    public void gnuplotBatchScript(IO io) {
        if (io.showInitialSlice == 0) {
            if ((listFile[riemannLocation].getName() == null ? riemann == null : listFile[riemannLocation].getName().equals(riemann))
                    && (listFile[approximateLocation].getName() == null ? approximate == null : listFile[approximateLocation].getName().equals(approximate))) {

                FileOutputStream writeToThisrho0BatchFile = null;
                try {
                    writeToThisrho0BatchFile = new FileOutputStream("rho0.plt");

                    PrintStream gnuplotScriptrho0 = new PrintStream(writeToThisrho0BatchFile);

                    if (io.showExact == 1) {
                        gnuplotScriptrho0.printf("set ylabel \"The Rest Density (rho0) in Sod Units\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:8 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:3 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:6 with linespoints title \"Difference\"\n"
                                + "set output \"" + rho0 + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptrho0.printf("set title \"Approximate Rest Density.\"\n"
                                + "set ylabel \"rho0\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:8 title \"Approximate Solution\"\n"
                                + "set output \"" + rho0 + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisrho0BatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThispBatchFile = null;
                try {
                    writeToThispBatchFile = new FileOutputStream("p.plt");
                    PrintStream gnuplotScriptp = new PrintStream(writeToThispBatchFile);
                    if (io.showExact == 1) {
                        gnuplotScriptp.printf("set ylabel \"The Pressure in Sod Units\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:5 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:7 with linespoints title \"Difference\"\n"
                                + "set output \"" + p + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptp.printf("set title \"Approximate Pressure.\"\n"
                                + "set ylabel \"p\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\"\n"
                                + "set output \"" + p + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThispBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisuBatchFile = null;
                try {
                    writeToThisuBatchFile = new FileOutputStream("u.plt");
                    PrintStream gnuplotScriptu = new PrintStream(writeToThisuBatchFile);
                    if (io.showExact == 1) {
                        gnuplotScriptu.printf("set ylabel \"The Internal Energy in Sod Units\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:9 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:4 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:9 with linespoints title \"Difference\"\n"
                                + "set output \"" + u + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptu.printf("set title \"Approximate Internal Energy.\"\n"
                                + "set ylabel \"u\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:9 title \"Approximate Solution\"\n"
                                + "set output \"" + u + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisuBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisvBatchFile = null;
                try {
                    writeToThisvBatchFile = new FileOutputStream("v.plt");
                    PrintStream gnuplotScriptv = new PrintStream(writeToThisvBatchFile);
                    if (io.showExact == 1) {
                        gnuplotScriptv.printf("set ylabel \"The Fluid Velocity\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:7 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:2 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:8 with linespoints title \"Difference\"\n"
                                + "set output \"" + v + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptv.printf("set title \"Approximate Velocity.\"\n"
                                + "set ylabel \"v\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:7 title \"Approximate Solution\"\n"
                                + "set output \"" + v + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisvBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

            } else if ((listFile[approximateLocation].getName() == null ? approximate == null : listFile[approximateLocation].getName().equals(approximate))) {
                /// This part is used for KKT since there is no Riemann solution file(s) to be generated.
                /// So, setting showInitialSlice = 0 and showExact = 0.

                FileOutputStream writeToThisrho0BatchFile = null;
                try {
                    writeToThisrho0BatchFile = new FileOutputStream("rho0.plt");

                    PrintStream gnuplotScriptrho0 = new PrintStream(writeToThisrho0BatchFile);

                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptrho0.printf("set title \"Approximate Rest Density.\"\n"
                                + "set ylabel \"rho0\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:8 title \"Approximate Solution\"\n"
                                + "set output \"" + rho0 + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisrho0BatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThispBatchFile = null;
                try {
                    writeToThispBatchFile = new FileOutputStream("p.plt");
                    PrintStream gnuplotScriptp = new PrintStream(writeToThispBatchFile);
                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptp.printf("set title \"Approximate Pressure.\"\n"
                                + "set ylabel \"p\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\"\n"
                                + "set output \"" + p + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThispBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisuBatchFile = null;
                try {
                    writeToThisuBatchFile = new FileOutputStream("u.plt");
                    PrintStream gnuplotScriptu = new PrintStream(writeToThisuBatchFile);
                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptu.printf("set title \"Approximate Internal Energy.\"\n"
                                + "set ylabel \"u\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:9 title \"Approximate Solution\"\n"
                                + "set output \"" + u + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisuBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisvBatchFile = null;
                try {
                    writeToThisvBatchFile = new FileOutputStream("v.plt");
                    PrintStream gnuplotScriptv = new PrintStream(writeToThisvBatchFile);
                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptv.printf("set title \"Approximate Velocity.\"\n"
                                + "set ylabel \"v\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:7 title \"Approximate Solution\"\n"
                                + "set output \"" + v + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisvBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            } else {
                System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                        + "Can't show .plt files because there are no Approximate solutions generated.");
            }
        } else if (io.showInitialSlice == 1) {
            if ((listFile[riemannLocation].getName() == null ? riemann == null : listFile[riemannLocation].getName().equals(riemann))
                    && (listFile[approximateLocation].getName() == null ? approximate == null : listFile[approximateLocation].getName().equals(approximate))) {

                FileOutputStream writeToThisrho0BatchFile = null;
                try {
                    writeToThisrho0BatchFile = new FileOutputStream("rho0.plt");

                    PrintStream gnuplotScriptrho0 = new PrintStream(writeToThisrho0BatchFile);

                    if (io.showExact == 1) {
                        gnuplotScriptrho0.printf("set ylabel \"The Rest Density (rho0)\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:8 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:3 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:6 with linespoints title \"Difference\"\n"
                                + "set output \"" + rho0 + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptrho0.printf("set title \"Approximate Density.\"\n"
                                + "set ylabel \"rho0\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:8\n"
                                + "set output \"" + rho0 + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisrho0BatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThispBatchFile = null;
                try {
                    writeToThispBatchFile = new FileOutputStream("p.plt");
                    PrintStream gnuplotScriptp = new PrintStream(writeToThispBatchFile);
                    if (io.showExact == 1) {
                        gnuplotScriptp.printf("set ylabel \"p\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:5 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:7 with linespoints title \"Difference\"\n"
                                + "set output \"" + p + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptp.printf("set title \"Approximate Pressure.\"\n"
                                + "set ylabel \"p\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10\n"
                                + "set output \"" + p + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThispBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisuBatchFile = null;
                try {
                    writeToThisuBatchFile = new FileOutputStream("u.plt");
                    PrintStream gnuplotScriptu = new PrintStream(writeToThisuBatchFile);
                    if (io.showExact == 1) {
                        gnuplotScriptu.printf("set ylabel \"u\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:9 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:4 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:9 with linespoints title \"Difference\"\n"
                                + "set output \"" + u + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptu.printf("set title \"Approximate Internal Energy.\"\n"
                                + "set ylabel \"u\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:9\n"
                                + "set output \"" + u + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisuBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisvBatchFile = null;
                try {
                    writeToThisvBatchFile = new FileOutputStream("v.plt");
                    PrintStream gnuplotScriptv = new PrintStream(writeToThisvBatchFile);
                    if (io.showExact == 1) {
                        gnuplotScriptv.printf("set ylabel \"v\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:7 title \"Approximate Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:2 with lines title \"Exact Solution\", \"" + listFile[riemannLocation].getName() + "\" using 1:8 with linespoints title \"Difference\"\n"
                                + "set output \"" + v + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else if (io.showExact == 0) {
                        gnuplotScriptv.printf("set title \"Approximate Velocity.\"\n"
                                + "set ylabel \"v\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:7\n"
                                + "set output \"" + v + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0 to not show or 1 to show.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisvBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

            } else if ((listFile[approximateLocation].getName() == null ? approximate == null : listFile[approximateLocation].getName().equals(approximate))) {
                /// This part is used for KKT since there is no Riemann solution file(s) to be generated.
                /// So, setting showInitialSlice = 0 and showExact = 0.

                FileOutputStream writeToThisrho0BatchFile = null;
                try {
                    writeToThisrho0BatchFile = new FileOutputStream("rho0.plt");

                    PrintStream gnuplotScriptrho0 = new PrintStream(writeToThisrho0BatchFile);

                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptrho0.printf("set title \"Approximate Rest Density.\"\n"
                                + "set ylabel \"rho0\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:8 title \"Approximate Solution\"\n"
                                + "set output \"" + rho0 + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisrho0BatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThispBatchFile = null;
                try {
                    writeToThispBatchFile = new FileOutputStream("p.plt");
                    PrintStream gnuplotScriptp = new PrintStream(writeToThispBatchFile);
                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptp.printf("set title \"Approximate Pressure.\"\n"
                                + "set ylabel \"p\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:10 title \"Approximate Solution\"\n"
                                + "set output \"" + p + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThispBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisuBatchFile = null;
                try {
                    writeToThisuBatchFile = new FileOutputStream("u.plt");
                    PrintStream gnuplotScriptu = new PrintStream(writeToThisuBatchFile);
                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptu.printf("set title \"Approximate Internal Energy.\"\n"
                                + "set ylabel \"u\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:9 title \"Approximate Solution\"\n"
                                + "set output \"" + u + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisuBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                FileOutputStream writeToThisvBatchFile = null;
                try {
                    writeToThisvBatchFile = new FileOutputStream("v.plt");
                    PrintStream gnuplotScriptv = new PrintStream(writeToThisvBatchFile);
                    if (io.showExact == 1) {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. Make 0.");
                    } else if (io.showExact == 0) {
                        gnuplotScriptv.printf("set title \"Approximate Velocity.\"\n"
                                + "set ylabel \"v\"\n"
                                + "set xlabel \"Length of Shock Tube (km)\"\n"
                                + "plot \"" + listFile[approximateLocation].getName() + "\" using 2:7 title \"Approximate Solution\"\n"
                                + "set output \"" + v + "\"\n"
                                + "set terminal eps\n"
                                + "replot");
                    } else {
                        System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                                + "showExact is an unknown value. It is neither 0 nor 1.");
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        writeToThisvBatchFile.close();
                    } catch (IOException ex) {
                        Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            } else {
                System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                        + "Can't show .plt files because there are no Approximate solutions generated..");
            }
        } else {
            System.err.println("ERROR :: gnuplotBatchFileSetup ::\n"
                    + "showInitialSlice is an invalid value. Check your param.dat file.");
            System.exit(1);
        }
    }

    public void gnuplotBatchScript(String initialName) {
        FileOutputStream writeToThisInitialvelBatchFile = null;
        try {
            writeToThisInitialvelBatchFile = new FileOutputStream("initialVel.plt");

            PrintStream gnuplotScriptInitialvel = new PrintStream(writeToThisInitialvelBatchFile);
            gnuplotScriptInitialvel.printf("set title \"The Initial Velocity.\"\n"
                    + "set ylabel \"velocity\"\n"
                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                    + "plot \"" + initialName + "\" using 2:7 title \"Velocity\"\n"
                    + "set output \"" + initialVel + "\"\n"
                    + "set terminal eps\n"
                    + "replot");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                writeToThisInitialvelBatchFile.close();
            } catch (IOException ex) {
                Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        FileOutputStream writeToThisInitialrho0BatchFile = null;
        try {
            writeToThisInitialrho0BatchFile = new FileOutputStream("initialRho0.plt");

            PrintStream gnuplotScriptInitialrho0 = new PrintStream(writeToThisInitialrho0BatchFile);
            gnuplotScriptInitialrho0.printf("set title \"The Initial Rho0.\"\n"
                    + "set ylabel \"rho0\"\n"
                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                    + "plot \"" + initialName + "\" using 2:8 title \"Rest Mass\"\n"
                    + "set output \"" + initialRho0 + "\"\n"
                    + "set terminal eps\n"
                    + "replot");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                writeToThisInitialrho0BatchFile.close();
            } catch (IOException ex) {
                Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        FileOutputStream writeToThisInitialepsBatchFile = null;
        try {
            writeToThisInitialepsBatchFile = new FileOutputStream("initialEps.plt");

            PrintStream gnuplotScriptInitialeps = new PrintStream(writeToThisInitialepsBatchFile);
            gnuplotScriptInitialeps.printf("set title \"The Initial Energy.\"\n"
                    + "set ylabel \"energy\"\n"
                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                    + "plot \"" + initialName + "\" using 2:9 title \"Energy\"\n"
                    + "set output \"" + initialEps + "\"\n"
                    + "set terminal eps\n"
                    + "replot");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                writeToThisInitialepsBatchFile.close();
            } catch (IOException ex) {
                Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        FileOutputStream writeToThisInitialpBatchFile = null;
        try {
            writeToThisInitialpBatchFile = new FileOutputStream("initialP.plt");

            PrintStream gnuplotScriptInitialp = new PrintStream(writeToThisInitialpBatchFile);
            gnuplotScriptInitialp.printf("set title \"The Initial Pressure.\"\n"
                    + "set ylabel \"pressure\"\n"
                    + "set xlabel \"Length of Shock Tube (km)\"\n"
                    + "plot \"" + initialName + "\" using 2:10 title \"Pressure\"\n"
                    + "set output \"" + initialP + "\"\n"
                    + "set terminal eps\n"
                    + "replot");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                writeToThisInitialpBatchFile.close();
            } catch (IOException ex) {
                Logger.getLogger(gnuplotBatchFileSetup.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
}
