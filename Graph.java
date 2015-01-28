import java.awt.Color;

public class Graph {
	public final int foreground, background;
	public final double[][] data;
	public final String ordinate;

	public Graph(double[][] data, String ordinate) {
		this(data, ordinate, Color.BLACK, Color.WHITE);
	}

	public Graph(double[][] data, String ordinate, Color foreground, Color background) {
		this.data = data;
		this.ordinate = ordinate;
		this.foreground = foreground.getRGB();
		this.background = background.getRGB();
	}

}
