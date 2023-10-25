/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

/**
 *
 * @author gregory
 */
public class Main {

    public static void main(String[] args) throws Exception {
        greg_const gc = new greg_const();
        gnuplotBatchFileSetup gnuplot;

        eqs eqs = new eqs();

        geometrize geo = new geometrize();
        exactRiemannSolution exactRieAndSDs = new exactRiemannSolution();

        stateVariables lstate = new stateVariables();
        stateVariables rstate = new stateVariables();

        grid grid = new grid();
        evolve doThis = new evolve();

        timeData ptime = new timeData();
        lastDump lastdump = new lastDump();

        double holdTime = 0.0;

        IO io = new IO(eqs, lstate, rstate);
        artVisControl artvis = new artVisControl(io);
        neutrinoMethods neutrinoM = new neutrinoMethods(gc);

        ptime.Make(io);
        /**
         * Initialize NOTE: initial data is read by the Grid class. Use gnuplot
         * to graph the initial data to check for correctness.
         */
        if (io.showInitialSlice == 1) {
            grid.MakeInitial(eqs, io, gc, geo, neutrinoM);
            String initialForms = doThis.PutSlice(io, ptime, geo, grid);
            gnuplot = new gnuplotBatchFileSetup();
            gnuplot.gnuplotBatchScript(initialForms);
            /// Dump the initial slice
            doThis.PutSlice(io, ptime, geo, grid);
            /// Print to console the IO values.
            //io.printIO(lstate, rstate);
        } else if (io.showInitialSlice == 0) {
            grid.MakeInitial(eqs, io, gc, geo, neutrinoM);
            /// Print to console the IO values.
            //io.printIO(lstate, rstate);
        } else {
            System.err.println("ERROR :: Main :: \n"
                    + "showInitialSlice in param.dat is an unrecognized value. \n"
                    + "Please emake 0 => Do not show or 1 => Show.");
            System.exit(1);
        }

        /// "LastDump" class decides on frequency of complete time slice dumps
        lastdump.LastDump(io, ptime);

        /*
         * StopWatch times the execution in milliseconds.
         */
        StopWatch timer = new StopWatch();
        timer.start();

        ///=============================================================================
        /// Evolve
        ///=============================================================================
        for (int n_time = 1; n_time <= ptime.n_max; n_time++) {
            System.out.print(".");

            ptime.makeDelta(eqs, grid);
            
            if (io.method == 0) {
                /// Evolve using upwind.
                doThis.Evolve(0, artvis, eqs, io, neutrinoM, grid, ptime);
            } else if (io.method == 1) {
                /// Evolve using Godunov.
                doThis.Evolve(1, artvis, eqs, io, neutrinoM, grid, ptime);
            } else {
                System.err.println("ERROR :: In Main :: Invalid switch for the method."
                        + " Choose either 1 or 0 according to comments in the param.dat file.");
            }

            holdTime = timer.getElapsedTime();

            /// update all times
            ptime.Update();

            /// Check for end of run (physical or cpu time) reached.

            switch (ptime.testForMax()) {
                case timeData.OUT_OF_TIME:
                    
                    if (io.dataSwitch == 0) {
                        if (io.showExact == 1) {
                            exactRieAndSDs.MakeExactRieandSDSod(io, eqs, ptime.coord, grid, lstate, rstate);
                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate1 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate1);
                            System.exit(0);
                        } else if (io.showExact == 0) {
                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate1 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate1);
                            System.exit(0);
                        } else {
                            System.err.println("ERROR :: Main ::\n"
                                    + "io.showExact is an unrecognized value. Make it 1 to show the exact solution\n"
                                    + "or 0 to suppress that graphic.");
                            System.exit(1);
                        }
                    } else if (io.dataSwitch == 1) {
                        if (io.showExact == 1) {
                            exactRieAndSDs.MakeExactRieandSDExperimental(io, eqs, ptime.coord, grid, lstate, rstate);

                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate1 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate1);
                            System.exit(0);
                        } else if (io.showExact == 0) {
                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate1 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate1);
                            System.exit(0);
                        } else {
                            System.err.println("ERROR :: Main ::\n"
                                    + "io.showExact is an unrecognized value. Make it 1 to show the exact solution\n"
                                    + "or 0 to suppress that graphic.");
                            System.exit(1);
                        }
                    } else if (io.dataSwitch == 2) {
                        if (io.showExact == 1) {
                            exactRieAndSDs.MakeExactRieandSDKKT(io, eqs, ptime.coord, geo, grid, ptime, lstate, rstate);

                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate1 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate1);
                            System.exit(0);
                        } else if (io.showExact == 0) {
                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate1 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate1);
                            System.exit(0);
                        } else {
                            System.err.println("ERROR :: Main ::\n"
                                    + "io.showExact is an unrecognized value. Make it 1 to show the exact solution\n"
                                    + "or 0 to suppress that graphic.");
                            System.exit(1);
                        }
                    } else {
                        System.err.println("ERROR :: Main :: \n"
                                + "Unrecognized dataSwitch. Check your param.dat file.");
                        System.exit(1);
                    }
                    break;
                case timeData.OUT_OF_TIME_STEPS:
                    if (io.dataSwitch == 0) {
                        if (io.showExact == 1) {
                            exactRieAndSDs.MakeExactRieandSDSod(io, eqs, ptime.coord, grid, lstate, rstate);

                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate2 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate2);
                            System.exit(0);
                        } else if (io.showExact == 0) {
                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate2 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate2);
                            System.exit(0);
                        } else {
                            System.err.println("ERROR :: Main ::\n"
                                    + "io.showExact is an unrecognized value. Make it 1 to show the exact solution\n"
                                    + "or 0 to suppress that graphic.");
                            System.exit(1);
                        }
                    } else if (io.dataSwitch == 1) {
                        if (io.showExact == 1) {
                            exactRieAndSDs.MakeExactRieandSDExperimental(io, eqs, ptime.coord, grid, lstate, rstate);

                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate2 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate2);
                            System.exit(0);
                        } else if (io.showExact == 0) {
                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate2 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate2);
                            System.exit(0);
                        } else {
                            System.err.println("ERROR :: Main ::\n"
                                    + "io.showExact is an unrecognized value. Make it 1 to show the exact solution\n"
                                    + "or 0 to suppress that graphic.");
                            System.exit(1);
                        }
                    } else if (io.dataSwitch == 2) {
                        if (io.showExact == 1) {
                            exactRieAndSDs.MakeExactRieandSDKKT(io, eqs, ptime.coord, geo, grid, ptime, lstate, rstate);

                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate2 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate2);
                            System.exit(0);
                        } else if (io.showExact == 0) {
                            System.out.println("main: Normal Exit.");
                            System.out.println("The elapsed time for the program is " + (holdTime / 1000) + " seconds.");
                            String approximate2 = doThis.PutSlice(io, ptime, geo, grid);
                            timer.stop();
                            gnuplot = new gnuplotBatchFileSetup(io, approximate2);
                            System.exit(0);
                        } else {
                            System.err.println("ERROR :: Main ::\n"
                                    + "io.showExact is an unrecognized value. Make it 1 to show the exact solution\n"
                                    + "or 0 to suppress that graphic.");
                            System.exit(1);
                        }
                    } else {
                        System.err.println("ERROR :: Main :: \n"
                                + "Unrecognized dataSwitch. Check your param.dat file.");
                        System.exit(1);
                    }
                    break;
            }

            /// time step too small?

            if (ptime.delta < ptime.delta_min) {
                System.err.println("main: ERROR: time step too small\n" + grid + "\n");
                doThis.PutSlice(io, ptime, geo, grid);
                System.exit(1);
            }

            /// Check to see if a time slice dump is due.

            lastdump.Update();
            if (lastdump.Test(ptime) == lastdump.DUMP) {
                doThis.PutSlice(io, ptime, geo, grid);
                lastdump.Reset(ptime);
            }

            if (ptime.coord >= ptime.coord_max) {
                System.out.println("Successful end of run.");
                /// dump the last slice to a standard final configuration file
                doThis.PutSlice(io, ptime, geo, grid);
                System.exit(0);
            }
        }   /// end of evolution loop
        System.err.println("ERROR :: Max number of time steps reached.");
        doThis.PutSlice(io, ptime, geo, grid);
    }
}
