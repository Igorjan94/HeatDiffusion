import java.awt.Color;

public class Graph {
    public final int foreground, background;
    public final double[][] data;
    public final String ordinate;
    public final double min, max;

    public Graph(double[][] data, String ordinate) {
        this(data, ordinate, Color.BLACK, Color.WHITE);
    }

    public Graph(double[][] data, double min, double max, String ordinate) {
        this(data, min, max, ordinate, Color.BLACK, Color.WHITE);
    }

    public Graph(double[][] data, String ordinate, Color foreground, Color background) {
        double min = Double.POSITIVE_INFINITY, max = -min;

        for (double[] layer : data) {
            for (double val : layer) {
                min = Math.min(min, val);
                max = Math.max(max, val);
            }
        }
        this.data = data;
        this.min = min;
        this.max = max;
        this.ordinate = ordinate;
        this.foreground = foreground.getRGB();
        this.background = background.getRGB();
    }

    public Graph(double[][] data, double min, double max, String ordinate, Color foreground, Color background) {
        this.data = data;
        this.min = min;
        this.max = max;
        this.ordinate = ordinate;
        this.foreground = foreground.getRGB();
        this.background = background.getRGB();
    }

}
