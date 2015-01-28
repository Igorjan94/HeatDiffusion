import static java.lang.Math.*;

import java.awt.Color;

public class ReactionSolover {

    public static PhysicalValue[] getVariables() {
        return new PhysicalValue[]{
                // These are left for sample
                new PhysicalValue("a", 0.007),
                new PhysicalValue("b", 2),
                new PhysicalValue("c", 3),
                new PhysicalValue("time", 1000, "sec"),
                new PhysicalValue("length", 1000, "px"),

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
                new PhysicalValue("dz",     1e-4,               "m"),
                new PhysicalValue("dt",     1e-4,               "sec"),
                new PhysicalValue("l",      1,                  "m"),
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

    public static Graph[] solove(PhysicalValue[] physicalValues) {
        double a = findValue(physicalValues, "a").value;
        double b = findValue(physicalValues, "b").value;
        double c = findValue(physicalValues, "c").value;
        int n = findValue(physicalValues, "time").intValue();
        int m = findValue(physicalValues, "length").intValue();

        double[][][] d = new double[3][n][m];

        for (int t = 0; t < n; t++) {
            for (int x = 0; x < m; x++) {
                d[0][t][x] = sin(a * x + b * t);
                d[1][t][x] = cos(a * x + b * t + c * x * t);
                d[2][t][x] = d[0][t][x] + d[1][t][x];
            }
        }

        Graph[] graph = new Graph[]{
                new Graph(d[0], "sin"),
                new Graph(d[1], "cos"),
                new Graph(d[2], "sin + cos"),
        };

        return graph;
    }

}
