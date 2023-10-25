/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Gregory Mohammed
 */
public class exactRiemannSolution {

    riemannExact rexact1;
    riemannExact rexact2;
    riemannExact rexact3;
    exactNode enode;
    public double sdrho0;
    public double sdp;
    public double sdv;
    public double sdu;
    double diffrho0;
    double diffp;
    double diffv;
    double diffu;
    double sumrho0;
    double sump;
    double sumu;
    double sumv;
    double[] rho0i;
    double[] pi;
    double[] ui;
    double[] vi;

    public exactRiemannSolution() {
        rexact1 = new riemannExact();
        rexact2 = new riemannExact();
        rexact3 = new riemannExact();
        enode = new exactNode();

        sdrho0 = 0.0;
        sdp = 0.0;
        sdv = 0.0;
        sdu = 0.0;

        diffrho0 = 0.0;
        diffp = 0.0;
        diffv = 0.0;
        diffu = 0.0;

        sumrho0 = 0.0;
        sump = 0.0;
        sumu = 0.0;
        sumv = 0.0;
    }

    public void MakeExactRieandSDSod(IO io, eqs eqs, double totalTime, grid lastGrid,
            stateVariables lstate, stateVariables rstate) throws Exception {
        try {
            String cellsfilename = "Riemann_Exact_Solution.csv";
            FileOutputStream writeToThisFileForCells = new FileOutputStream(cellsfilename);
            PrintStream dataToBeWrittenForCells = new PrintStream(writeToThisFileForCells);
            String header0 = "#x(1) v(2) rho0(3) u(4) p(5) diffrho0(6) diffp(7) diffv(8) diffu(9) \n";
            dataToBeWrittenForCells.printf(header0);

            double x_discon = 0.0;

            rexact1.Solve(totalTime, x_discon, eqs, io.forRiemannLeft, io.forRiemannRight);

            ArrayList<exactNode> nlist = rexact1.MakeMeshValues(eqs, lastGrid.clist.size(), io.lx, io.rx, lstate, rstate);
            
            rho0i = new double[lastGrid.clist.size()];
            pi = new double[lastGrid.clist.size()];
            ui = new double[lastGrid.clist.size()];
            vi = new double[lastGrid.clist.size()];

            for (int i = 0; i < lastGrid.clist.size(); i++) {

                /// Values at nodes for graphing purposes
                //double x = lastGrid.nlist.get(i).x;

                //rexact1.GetExactValues(x, eqs, enode, io.forRiemannLeft, io.forRiemannRight);

                /// Calculate positive differences
                diffrho0 = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._rho0 - nlist.get(i).rho0), 2.0));
                diffp = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._p - nlist.get(i).p), 2.0));
                diffv = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._v - nlist.get(i).vel), 2.0));
                diffu = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._eps - nlist.get(i).u), 2.0));

                sumrho0 = sumrho0 + diffrho0;
                sump = sump + diffp;
                sumv = sumv + diffv;
                sumu = sumu + diffu;

                rho0i[i] = diffrho0;
                pi[i] = diffp;
                vi[i] = diffv;
                ui[i] = diffu;

                dataToBeWrittenForCells.printf(lastGrid.clist.get(i).xmid(i, lastGrid.nlist) + " " + nlist.get(i).vel
                        + " " + nlist.get(i).rho0 + " " + nlist.get(i).u + " " + nlist.get(i).p
                        + " " + diffrho0 + " " + diffp + " " + diffv + " " + diffu
                        + "\n");
                dataToBeWrittenForCells.flush();
            }

            for (int i = 0; i < lastGrid.clist.size(); i++) {
                sdrho0 = Math.sqrt((Math.pow(((sumrho0 / io.resolution) - rho0i[i]), 2.0)) / io.resolution);
                sdp = Math.sqrt((Math.pow(((sump / io.resolution) - pi[i]), 2.0)) / io.resolution);
                sdu = Math.sqrt((Math.pow(((sumu / io.resolution) - ui[i]), 2.0)) / io.resolution);
                sdv = Math.sqrt((Math.pow(((sumv / io.resolution) - vi[i]), 2.0)) / io.resolution);
            }

            String sdfilename = "standard_deviations.csv";
            FileOutputStream writeToSDFileForCells = new FileOutputStream(sdfilename);
            PrintStream dataToBeWrittenForSD = new PrintStream(writeToSDFileForCells);
            String header = "#Resolution(1) sdv(2) sdrho0(3) sdu(4) sdp(5) \n";
            dataToBeWrittenForSD.printf(header);

            dataToBeWrittenForSD.printf(io.resolution + " " + sdv + " " + sdrho0 + " " + sdu + " " + sdp + "\n");

            dataToBeWrittenForSD.close();
            dataToBeWrittenForCells.close();
        } catch (IOException ex) {
            System.err.println(ex);
            Logger.getLogger(grid.class.getName()).log(Level.SEVERE, null, ex);
        }
        /// =================================================================
    }

    public void MakeExactRieandSDExperimental(IO io, eqs eqs, double totalTime, grid lastGrid,
            stateVariables lstate, stateVariables rstate) throws Exception {
        try {
            String cellsfilename = "Riemann_Exact_Solution.csv";
            FileOutputStream writeToThisFileForCells = new FileOutputStream(cellsfilename);
            PrintStream dataToBeWrittenForCells = new PrintStream(writeToThisFileForCells);
            String header0 = "#x(1) v(2) rho0(3) u(4) p(5) squarediffrho0(6) squarediffp(7) squarediffv(8) squarediffu(9) \n";
            dataToBeWrittenForCells.printf(header0);

            double x_discon = 0.0;

            rexact2.Solve(totalTime, x_discon, eqs, lstate, rstate);

            rho0i = new double[lastGrid.clist.size()];
            pi = new double[lastGrid.clist.size()];
            ui = new double[lastGrid.clist.size()];
            vi = new double[lastGrid.clist.size()];

            for (int i = 0; i < lastGrid.clist.size(); i++) {

                /// Values at nodes for graphing purposes
                double x = lastGrid.clist.get(i)._x;

                rexact2.GetExactValues(x, eqs, enode, lstate, rstate);

                /// Calculate positive differences
                diffrho0 = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._rho0 - enode.rho0), 2.0));
                diffp = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._p - enode.p), 2.0));
                diffv = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._v - enode.vel), 2.0));
                diffu = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._eps - enode.u), 2.0));

                sumrho0 = sumrho0 + diffrho0;
                sump = sump + diffp;
                sumv = sumv + diffv;
                sumu = sumu + diffu;

                rho0i[i] = diffrho0;
                pi[i] = diffp;
                vi[i] = diffv;
                ui[i] = diffu;

                dataToBeWrittenForCells.printf(x + " " + enode.vel
                        + " " + enode.rho0 + " " + enode.u + " " + enode.p
                        + " " + diffrho0 + " " + diffp + " " + diffv + " " + diffu
                        + "\n");
                dataToBeWrittenForCells.flush();
            }

            for (int i = 0; i < lastGrid.clist.size(); i++) {
                sdrho0 = Math.sqrt((Math.pow(((sumrho0 / io.resolution) - rho0i[i]), 2.0)) / io.resolution);
                sdp = Math.sqrt((Math.pow(((sump / io.resolution) - pi[i]), 2.0)) / io.resolution);
                sdu = Math.sqrt((Math.pow(((sumu / io.resolution) - ui[i]), 2.0)) / io.resolution);
                sdv = Math.sqrt((Math.pow(((sumv / io.resolution) - vi[i]), 2.0)) / io.resolution);
            }

            String sdfilename = "standard_deviations.csv";
            FileOutputStream writeToSDFileForCells = new FileOutputStream(sdfilename);
            PrintStream dataToBeWrittenForSD = new PrintStream(writeToSDFileForCells);
            String header = "#Resolution(1) sdv(2) sdrho0(3) sdu(4) sdp(5) \n";
            dataToBeWrittenForSD.printf(header);

            dataToBeWrittenForSD.printf(io.resolution + " " + sdv + " " + sdrho0 + " " + sdu + " " + sdp + "\n");

            dataToBeWrittenForSD.close();
            dataToBeWrittenForCells.close();
        } catch (IOException ex) {
            System.err.println(ex);
            Logger.getLogger(grid.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void MakeExactRieandSDKKT(IO io, eqs eqs, double totalTime, geometrize geo, grid lastGrid, timeData ptime,
            stateVariables lstate, stateVariables rstate) throws Exception {
        try {
            String cellsfilename = "Riemann_Exact_Solution.csv";
            FileOutputStream writeToThisFileForCells = new FileOutputStream(cellsfilename);
            PrintStream dataToBeWrittenForCells = new PrintStream(writeToThisFileForCells);
            String header0 = "#x(1) v(2) rho0(3) u(4) p(5) squarediffrho0(6) squarediffp(7) squarediffv(8) squarediffu(9) \n";
            dataToBeWrittenForCells.printf(header0);

            //double T = 0.25;
            double x_discon = 0.0;

            rexact3.Solve(totalTime, x_discon, eqs, lstate, rstate);

            rho0i = new double[lastGrid.clist.size()];
            pi = new double[lastGrid.clist.size()];
            ui = new double[lastGrid.clist.size()];
            vi = new double[lastGrid.clist.size()];

            if (io.showExact == 0) {
                /// Plot un-geometrized ...
                for (int i = 0; i < lastGrid.clist.size(); i++) {

                    /// Values at nodes for graphing purposes
                    double x = lastGrid.clist.get(i)._x;

                    if (i == (lastGrid.clist.size() / io.resolution)) {
                        lstate.rho0 = lastGrid.clist.get(i - 1).cur._rho0;
                        lstate.p = lastGrid.clist.get(i - 1).cur._p;
                        lstate.u = lastGrid.clist.get(i - 1).cur._eps;
                        lstate.vel = lastGrid.clist.get(i - 1).cur._v;
                        rstate.rho0 = lastGrid.clist.get(i + 1).cur._rho0;
                        rstate.p = lastGrid.clist.get(i + 1).cur._p;
                        rstate.u = lastGrid.clist.get(i + 1).cur._eps;
                        rstate.vel = lastGrid.clist.get(i + 1).cur._v;
                        rexact3.GetExactValues(x, eqs, enode, lstate, rstate);
                    }

                    /// Calculate positive differences
                    diffrho0 = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._rho0 - enode.rho0), 2.0));
                    diffp = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._p - enode.p), 2.0));
                    diffv = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._v - enode.vel), 2.0));
                    diffu = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._eps - enode.u), 2.0));

                    sumrho0 = sumrho0 + diffrho0;
                    sump = sump + diffp;
                    sumv = sumv + diffv;
                    sumu = sumu + diffu;

                    rho0i[i] = diffrho0;
                    pi[i] = diffp;
                    vi[i] = diffv;
                    ui[i] = diffu;

                    geo.makeUnGeometric(0, enode.rho0, enode.p, enode.u, enode.vel, x);

                    dataToBeWrittenForCells.printf(geo.xUn + " " + geo.vUn
                            + " " + geo.rho0Un + " " + geo.epsUn + " " + geo.pUn
                            + " " + diffrho0 + " " + diffp + " " + diffv + " " + diffu
                            + "\n");
                    dataToBeWrittenForCells.flush();
                }

                for (int i = 0; i < lastGrid.clist.size(); i++) {
                    sdrho0 = Math.sqrt((Math.pow(((sumrho0 / io.resolution) - rho0i[i]), 2.0)) / io.resolution);
                    sdp = Math.sqrt((Math.pow(((sump / io.resolution) - pi[i]), 2.0)) / io.resolution);
                    sdu = Math.sqrt((Math.pow(((sumu / io.resolution) - ui[i]), 2.0)) / io.resolution);
                    sdv = Math.sqrt((Math.pow(((sumv / io.resolution) - vi[i]), 2.0)) / io.resolution);
                }
            } else if (io.showExact == 1) {
                /// Plot geometrized ...
                for (int i = 0; i < lastGrid.clist.size(); i++) {

                    /// Values at nodes for graphing purposes
                    double x = lastGrid.clist.get(i)._x;

                    if (i == (lastGrid.clist.size() / io.resolution)) {
                        lstate.rho0 = lastGrid.clist.get(i - 1).cur._rho0;
                        lstate.p = lastGrid.clist.get(i - 1).cur._p;
                        lstate.u = lastGrid.clist.get(i - 1).cur._eps;
                        lstate.vel = lastGrid.clist.get(i - 1).cur._v;
                        rstate.rho0 = lastGrid.clist.get(i + 1).cur._rho0;
                        rstate.p = lastGrid.clist.get(i + 1).cur._p;
                        rstate.u = lastGrid.clist.get(i + 1).cur._eps;
                        rstate.vel = lastGrid.clist.get(i + 1).cur._v;
                        rexact3.GetExactValues(x, eqs, enode, lstate, rstate);
                    }

                    /// Calculate positive differences
                    diffrho0 = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._rho0 - enode.rho0), 2.0));
                    diffp = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._p - enode.p), 2.0));
                    diffv = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._v - enode.vel), 2.0));
                    diffu = Math.sqrt(Math.pow((lastGrid.clist.get(i).cur._eps - enode.u), 2.0));

                    sumrho0 = sumrho0 + diffrho0;
                    sump = sump + diffp;
                    sumv = sumv + diffv;
                    sumu = sumu + diffu;

                    rho0i[i] = diffrho0;
                    pi[i] = diffp;
                    vi[i] = diffv;
                    ui[i] = diffu;

                    geo.makeUnGeometric(0, enode.rho0, enode.p, enode.u, enode.vel, x);

                    dataToBeWrittenForCells.printf(geo.xUn + " " + geo.vUn
                            + " " + geo.rho0Un + " " + geo.epsUn + " " + geo.pUn
                            + " " + diffrho0 + " " + diffp + " " + diffv + " " + diffu
                            + "\n");
                    dataToBeWrittenForCells.flush();
                }

                for (int i = 0; i < lastGrid.clist.size(); i++) {
                    sdrho0 = Math.sqrt((Math.pow(((sumrho0 / io.resolution) - rho0i[i]), 2.0)) / io.resolution);
                    sdp = Math.sqrt((Math.pow(((sump / io.resolution) - pi[i]), 2.0)) / io.resolution);
                    sdu = Math.sqrt((Math.pow(((sumu / io.resolution) - ui[i]), 2.0)) / io.resolution);
                    sdv = Math.sqrt((Math.pow(((sumv / io.resolution) - vi[i]), 2.0)) / io.resolution);
                }
            } else {
                System.err.println("ERROR :: exactRiemannSolution : MakeExactRieandSDKKT :\n"
                        + "Unknown value for plotGeo. Check your param.dat file under plotGeo.");
                System.exit(1);
            }
            String sdfilename = "standard_deviations.csv";
            FileOutputStream writeToSDFileForCells = new FileOutputStream(sdfilename);
            PrintStream dataToBeWrittenForSD = new PrintStream(writeToSDFileForCells);
            String header = "#Resolution(1) sdv(2) sdrho0(3) sdu(4) sdp(5) \n";
            dataToBeWrittenForSD.printf(header);

            dataToBeWrittenForSD.printf(io.resolution + " " + sdv + " " + sdrho0 + " " + sdu + " " + sdp + "\n");

            dataToBeWrittenForSD.close();
            dataToBeWrittenForCells.close();
        } catch (IOException ex) {
            System.err.println(ex);
            Logger.getLogger(grid.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
