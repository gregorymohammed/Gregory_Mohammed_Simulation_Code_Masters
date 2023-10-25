package Gregory_Mohammed_Simulation_Code_Masters;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Gregory Mohammed
 */
public class evolve {

    public int corrector_level;
    public int n_corrector;
    public int verify;
    public int UpwindorGodunov;

    /**
     * ==================================================================
     */
    public evolve() {
        corrector_level = 0;
        n_corrector = 0;
        verify = 0;
        UpwindorGodunov = 0;
    }

    /**
     * ==================================================================
     */
    /**
     * ===================================================================
     * Switch to determine if an Upwind or Godunov evolution.
     * ===================================================================
     */
    public void Evolve(int choice, artVisControl artvis, eqs eqs, IO io, neutrinoMethods neutrinoM, grid existingGrid, timeData ptime) {
        switch (choice) {
            case 0:
                evolveUpwind(artvis, eqs, io, neutrinoM, existingGrid, ptime);
                break;
            case 1:
                evolveGod(artvis, eqs, io, neutrinoM, existingGrid, ptime);
                break;
            default:
                System.err.println("grid::Evolve: ERROR: unknown evolution.");
                System.exit(1);
        }
    }

    /**
     * ================================================================== THE
     * EVOLUTION METHODS FOR THE GODUNOV.
     * ==================================================================
     */
    public void evolveGod(artVisControl artvis, eqs eqs, IO io, neutrinoMethods neutrinoM, grid existingGrid, timeData ptime) {
        /**
         * =================================================== Exchange
         * old/current data (only the elements need this).
         * ===================================================
         */
        existingGrid.gridExchange();
        /**
         * =================================================== Evolve the
         * elements and nodes.
         * ===================================================
         */
        existingGrid.makeRiemannFluxPredictor(ptime, artvis, eqs, io);
        existingGrid.predict(neutrinoM, ptime, io);

        existingGrid.makePrimitiveCur(eqs);
        existingGrid.makeSigma(io);

        existingGrid.makeRiemannFluxCorrector(ptime, artvis, eqs, io);
        existingGrid.correct(neutrinoM, ptime, io);

        existingGrid.makePrimitiveCur(eqs);
        existingGrid.makeSigma(io);

        /**
         * =================================================== There can be
         * oscillations and difficulties so do some very simple fixups here.
         * ===================================================
         */
        for (int ie = 0; ie < existingGrid.clist.size(); ++ie) {
            if (existingGrid.clist.get(ie).cur._E < 0.0) {
                existingGrid.clist.get(ie).cur._E = 0.0;
            }
            if (existingGrid.clist.get(ie).cur._eps < 0.0) {
                existingGrid.clist.get(ie).cur._eps = 0.0;
            }
        }
    }

    /**
     * ================================================================== THE
     * EVOLUTION METHODS FOR THE UPWIND.
     * ==================================================================
     */
    public void evolveUpwind(artVisControl artvis, eqs eqs, IO io, neutrinoMethods neutrinoM,
            grid existingGrid, timeData ptime) {
        /**
         * exchange old/current data (only the elements need this)
         */
        existingGrid.gridExchange();
        /**
         * Predict
         */
        existingGrid.makeFluxPredictor(artvis, eqs, ptime);
        existingGrid.predict(neutrinoM, ptime, io);
        existingGrid.makePrimitiveCur(eqs);
        existingGrid.makeSigma(io);
        /**
         * Correct
         */
        for (corrector_level = 1; corrector_level <= io.corrector; ++corrector_level) {
            existingGrid.makeFluxCorrector(artvis, eqs, ptime);

            existingGrid.correct(neutrinoM, ptime, io);
            existingGrid.makePrimitiveCur(eqs);

            existingGrid.makeSigma(io);
        }
        /**
         * There can be oscillations and difficulties so do some very simple
         * fixups here.
         */
        for (int ie = 0; ie < existingGrid.clist.size(); ++ie) {
            if (existingGrid.clist.get(ie).cur._E < 0.0) {
                existingGrid.clist.get(ie).cur._E = 0.0;
            }
            if (existingGrid.clist.get(ie).cur._eps < 0.0) {
                existingGrid.clist.get(ie).cur._eps = 0.0;
            }
        }
    }

    /**
     * =================================================================== Put
     * Write the evolution values into a file.
     * ===================================================================
     */
    public String PutSlice(IO io, timeData ptime, geometrize geo, grid existingGrid) {
        if (io.plotGeo == 0) {
            String cellsfilename = String.format("ShockType" + io.initShock + "_Method" + io.method + "_Neutrinos" + io.neutrinos + "_ReconType" + io.reconstructionType + "_Res" + io.resolution + "_t%08d.csv", ptime.n);
            try {
                FileOutputStream writeToThisFileForCells = new FileOutputStream(cellsfilename);
                PrintStream dataToBeWrittenForCells = new PrintStream(writeToThisFileForCells);
                String header0 = "#CellList\n" + "#Coordinate Time=" + ptime.coord + "\n" + "#Id(1) x(2) Volume(3) S(4) D(5) E(6) 3-Vel(7) rho0(8) eps(9) p(10) left(12) right(13) \n"; //a(11)
                dataToBeWrittenForCells.printf(header0);
                
                for (int i = 0; i < existingGrid.clist.size(); i++) {
                    geo.makeUnGeometric(1,
                            existingGrid.clist.get(i).cur._rho0,
                            existingGrid.clist.get(i).cur._p,
                            existingGrid.clist.get(i).cur._eps,
                            existingGrid.clist.get(i).cur._v,
                            existingGrid.clist.get(i)._x);

                    String header1 = existingGrid.clist.get(i)._id
                            + " " + existingGrid.clist.get(i).xmid(i, existingGrid.nlist)
                            + " " + existingGrid.clist.get(i).Volume(i, existingGrid.nlist)
                            + " " + existingGrid.clist.get(i).cur._S
                            + " " + existingGrid.clist.get(i).cur._D
                            + " " + existingGrid.clist.get(i).cur._E
                            + " " + geo.vUn
                            + " " + geo.rho0Un
                            + " " + geo.epsUn
                            + " " + geo.pUn
                            + " " + existingGrid.nlist.get(i).id
                            + " " + existingGrid.nlist.get(i + 1).id
                            + "\n"; //+ " " + geo.aUn 
                    dataToBeWrittenForCells.printf(header1);
                    dataToBeWrittenForCells.flush();
                }
                dataToBeWrittenForCells.close();

            } catch (IOException ex) {
                System.err.println(ex);
                Logger.getLogger(grid.class.getName()).log(Level.SEVERE, null, ex);
            }
            return cellsfilename;
        } else if (io.plotGeo == 1) {
            String cellsfilename = String.format("ShockType" + io.initShock + "_Method" + io.method + "_Neutrinos" + io.neutrinos + "_ReconType" + io.reconstructionType + "_Res" + io.resolution + "_t%08d.csv", ptime.n);
            try {
                FileOutputStream writeToThisFileForCells = new FileOutputStream(cellsfilename);
                PrintStream dataToBeWrittenForCells = new PrintStream(writeToThisFileForCells);
                String header0 = "#CellList\n" + "#Coordinate Time=" + ptime.coord + "\n" + "#Id(1) x(2) Volume(3) S(4) D(5) E(6) 3-Vel(7) rho0(8) eps(9) p(10) left(12) right(13) \n"; //a(11)
                dataToBeWrittenForCells.printf(header0);
                
                for (int i = 0; i < existingGrid.clist.size(); i++) {
                    String header1 = existingGrid.clist.get(i)._id
                            + " " + existingGrid.clist.get(i).xmid(i, existingGrid.nlist)
                            + " " + existingGrid.clist.get(i).Volume(i, existingGrid.nlist)
                            + " " + existingGrid.clist.get(i).cur._S
                            + " " + existingGrid.clist.get(i).cur._D
                            + " " + existingGrid.clist.get(i).cur._E
                            + " " + existingGrid.clist.get(i).cur._v
                            + " " + existingGrid.clist.get(i).cur._rho0
                            + " " + existingGrid.clist.get(i).cur._eps
                            + " " + existingGrid.clist.get(i).cur._p
                            + " " + existingGrid.nlist.get(i).id
                            + " " + existingGrid.nlist.get(i + 1).id
                            + "\n"; //+ " " + geo.aUn 
                    dataToBeWrittenForCells.printf(header1);
                    dataToBeWrittenForCells.flush();
                }
                dataToBeWrittenForCells.close();

            } catch (IOException ex) {
                System.err.println(ex);
                Logger.getLogger(grid.class.getName()).log(Level.SEVERE, null, ex);
            }
            return cellsfilename;
        } else {
            System.err.println("ERROR :: evolve : PutSlice : Unknown unit selection.\n"
                    + "Check your param.dat file under plotGeo.");
            System.exit(1);
            return null;
        }
    }
}
