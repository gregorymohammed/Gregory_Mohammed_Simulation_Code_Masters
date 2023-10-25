/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Gregory_Mohammed_Simulation_Code_Masters;

/**
 *
 * @author gregory
 */
public class lastDump {
    //=============================================================================

    /** Local class to handle dump times during the evolution.

    Tests various quantities (interval since last dump,
    proper time elapsed, ...) and sets a flag if a dump is indicated.
     */
    int DUMP;
    int NODUMP;
    int interval;        // number of time steps since last dump
    double time;         // proper time at centre of last dump
    int interval_max;    // max number of times steps between dumps
    double time_max;     // max proper time at centre between dumps
    double epsdmp;

    public lastDump() {
        DUMP = 0;
        NODUMP = 1;
        interval = 0;
        time = 0.0;
        interval_max = 0;
        time_max = 0.0;
        epsdmp = 0.0;
    }

    public void Reset(timeData ptime) {
        interval = 0;
        time = ptime.coord;
    }

    public void Update() {
        interval++;
    }

    public void LastDump(IO io, timeData ptime) {
        epsdmp = io.epsdmp;
        interval_max = io.dmpinterval;

        time_max = epsdmp * ptime.coord_max;   /// time between dumps

        interval = 0;
        time = 0.0;
    }

    public int Test(timeData ptime) {
        if (interval >= interval_max) {
            return DUMP;
        }
        if (Math.abs(ptime.coord - time) > time_max) {
            return DUMP;
        }
        return NODUMP;
    }
}
