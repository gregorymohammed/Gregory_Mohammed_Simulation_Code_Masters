package Gregory_Mohammed_Simulation_Code_Masters;

public class artVisControl {

    /**
     * ==============================================
     * Just a struct.
     * ==============================================
     */
    public double _k1;
    public double _k2;
    /**
     * ==============================================
     * Constructor.
     * ==============================================
     */
    public artVisControl(IO io) {
        _k1 = io.artvisk1;
        _k2 = io.artvisk2;
    }
}
