package Gregory_Mohammed_Simulation_Code_Masters;

public class basic {

    /**
     * ==============================================
     * Constructor.
     * ==============================================
     */
    public basic() {
    }

    /**
     * ==============================================
     * Methods.
     * ==============================================
     */
    public double signum(double x) {
        if (x < 0.0) {
            return -1.0;
        } else if (x == 0.0) {
            return 0.0;
        } else if (x > 0.0) {
            return 1.0;
        } else {
            /**
             * Must be that x is Not-a-Number.
             */
            return x;
        }
    }

    public double min(double x1, double x2, double x3) {
        return MIN(MIN(x1, x2), x3);
    }

    public double max(double x1, double x2, double x3) {
        return MAX(MAX(x1, x2), x3);
    }

    public double MAX(final double a, final double b) {
        return (a > b ? a : b);
    }

    public double MIN(final double a, final double b) {
        return (a <= b ? a : b);
    }

    public double minmod(double a, double b, double c) {
        if ((a > 0.0)
                && (b > 0.0)
                && (c > 0.0)) {
            return +min(a, b, c);
        } else if ((a < 0.0)
                && (b < 0.0)
                && (c < 0.0)) {
            return max(a, b, c);
        } else {
            return 0.0;
        }
    }

//    public double signum(double a) {
//        return ((a) >= 0.0 ? (+1) : (-1));
//    }

    public double DSIGN(double a, double b) {
        return Math.abs(a) * signum(b);
    }
}
