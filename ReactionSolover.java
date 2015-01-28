import static java.lang.Math.*;

public class ReactionSolover {

    /*
    Equations:
    dX / dt = D * nabla^2 X - W(X, T)
    dT / dt = lambda / (rho * c) * nabla^2 T + rho * Q * W(X, T)
     */
    public static PhysicalValue[] getVariables() {
        return new PhysicalValue[]{
                // Definitely constants
                new PhysicalValue("E",      8 * 1e4,            "J/mol"),
                new PhysicalValue("R",      8.314,              "J/(mol * K)"),
                new PhysicalValue("alpha",  1),
                new PhysicalValue("Q",      7 * 1e5,            "J/kg"),
                new PhysicalValue("rho",    830,                "kg/m^3"),
                new PhysicalValue("T0",     293,                "K"),
                new PhysicalValue("c",      1980,               "J/(kg * K)"),
                new PhysicalValue("lambda", 0.13,               "J/(m * sec * K)"),

                // May be changed
                new PhysicalValue("K",      1.6 * 1e-6,         "1/sec"),
                new PhysicalValue("D",      0.13 / 830 / 1980,  "m^2 / sec"),    // initially: lambda / rho / c
                new PhysicalValue("Tw",     300,                "K"),
                new PhysicalValue("dt",     1e-4,               "sec"),
                new PhysicalValue("dz",     1e-4,               "m"),
                new PhysicalValue("l",      1,                  "m"),
                new PhysicalValue("steps",  1000),
        };
    }

    private static PhysicalValue findValue(PhysicalValue[] values, String name) {
        for (PhysicalValue v : values) {
            if (v.name.equals(name)) {
                return v;
            }
        }
        throw new RuntimeException("No such value found: " + name);
    }

    private static class WCalculator
    {
        WCalculator(double K, double alpha, double E, double R) {
            this.K = K;
            this.alpha = alpha;
            this.E = E;
            this.R = R;
        }

        public double getW(double X, double T) {
            return K * Math.pow(X, alpha) * Math.exp(-E / R / T);
        }

        private double K, alpha, E, R;
    }
    /* The simplest scheme
    X_i^{n+1} - X_i^{n} = dt * (D * (X_{i+1}^{n+1} - 2 * X_i^{n+1} + X_{i-1}^{n+1}) / dz^2 - K * X_{i}^{n+1} * (X_i_n)^(alpha - 1) * exp(-E  / (R * T_i^n))
    T_i^{n+1} - T_i^{n} = dt * (lambda / (rho * c) * (T^{n+1}_{i+1} - 2 * T^n_{i} + T^n_{i-1}) / dz^2 - K * Q / c * (x_i^{n+1})^alpha * exp(-E / (R * T_i^n))
     */
    public static Graph[] solove(PhysicalValue[] physicalValues) {
        double D = findValue(physicalValues, "D").value;
        double E = findValue(physicalValues, "E").value;
        double R = findValue(physicalValues, "R").value;
        double alpha = findValue(physicalValues, "alpha").value;
        double K = findValue(physicalValues, "K").value;
        double T0 = findValue(physicalValues, "T0").value;
        double Tw = findValue(physicalValues, "Tw").value;

        double dt = findValue(physicalValues, "dt").value;
        double dz = findValue(physicalValues, "dz").value;
        double l = findValue(physicalValues, "l").value;
        int z_steps = (int) Math.round(l / dz);
        int t_steps = findValue(physicalValues, "steps").intValue();

        double[][][] d = new double[3][t_steps][z_steps];
        double[][] X = d[0];
        double[][] T = d[1];
        double[][] W = d[2];

        WCalculator wCalculator = new WCalculator(K, alpha, E, R);

        for (int i = 0; i < z_steps; i++) {
            X[0][i] = i == 0 ? 1 : 0;
            T[0][i] = i == 0 ? Tw : T0;
            W[0][i] = wCalculator.getW(X[0][i], T[0][i]);
        }

        for (int n = 0; n < t_steps; n++) {
            X[n][0] = 0;
            T[n][0] = Tw;
            W[n][0] = wCalculator.getW(X[n][0], T[n][0]);

            X[n][z_steps - 1] = 1;
            T[n][z_steps - 1] = T0;
            W[n][z_steps - 1] = wCalculator.getW(X[n][z_steps - 1], T[n][z_steps - 1]);
        }

//        for (int t = 0; t < n; t++) {
//            for (int x = 0; x < m; x++) {
//                d[0][t][x] = sin(a * x + b * t);
//                d[1][t][x] = cos(a * x + b * t + c * x * t);
//            }
//        }

        return new Graph[]{
                new Graph(X, "X"),
                new Graph(T, "T"),
                new Graph(W, "W"),
        };
    }

}
