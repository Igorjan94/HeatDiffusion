import static java.lang.Math.max;
import static java.lang.Math.min;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.KeyEventDispatcher;
import java.awt.KeyboardFocusManager;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

public class Visualizer extends JFrame {
	static final long serialVersionUID = 0xE1A;

	JLabel label = new JLabel();
	BufferedImage canvas;
	Graph[] graph;

	int t = 0, n;

	int speed = 1, logSpeed = 0;
	int border = 20;

	void draw() {

		int cw = canvas.getWidth(), ch = canvas.getHeight();
		Graphics g = canvas.getGraphics();
		int gw = cw, gh = ch / n;

		for (int i = 0; i < n; i++) {

			for (int x = 0; x < gw; x++) {
				for (int y = 0; y < gh; y++) {
					canvas.setRGB(x, y + gh * i, graph[i].background);
				}
			}

			int m = graph[i].data.length;
			int cur = ((t % m) + m) % m;
			double[] d = graph[i].data[cur];

			double min = Double.POSITIVE_INFINITY, max = -min;
			for (double v : d) {
				min = Math.min(min, v);
				max = Math.max(max, v);
			}

			int h = Math.max(border * 3, gh - border * 2), w = Math.max(border * 3, gw - border / 2);

			double dw = d.length + 1, dh = (max - min) * 1.0001;

			for (int dx = 0; dx < d.length; dx++) {
				double dy = d[dx] - min;
				if (Double.isNaN(dy) || Double.isInfinite(dy)) {
					continue;
				}

				int x = (int) ((dx * w) / dw);
				if (0 <= x && x < w) {
					int y = (int) ((dy * h) / dh);
					if (0 <= y && y < h) {
						y = h - y - 1 + gh * i + border;
						x = x + border / 2;
						canvas.setRGB(x, y, graph[i].foreground);
					}
				}
			}
			g.setColor(new Color(graph[i].foreground));
			g.drawString(Double.toString(max), border / 2, gh * i + border);
			g.drawString(graph[i].ordinate, cw / 2, gh * i + border);
			g.drawString(Double.toString(min), border / 2, gh * (i + 1));

		}

		g.dispose();
		repaint();
	}

	public void onResize() {
		int cw = Math.max(200, getWidth()) - 14, ch = Math.max(300, getHeight()) - 37;
		canvas = new BufferedImage(cw, ch, BufferedImage.TYPE_INT_RGB);
		label.setIcon(new ImageIcon(canvas));

		draw();
	}

	public Visualizer(Graph[] graph) {
		this.graph = graph;
		n = graph.length;
		addComponentListener(new ComponentAdapter() {
			@Override
			public void componentResized(ComponentEvent e) {
				onResize();
			}
		});

		KeyboardFocusManager manager = KeyboardFocusManager.getCurrentKeyboardFocusManager();
		manager.addKeyEventDispatcher(new KeyEventDispatcher() {
			@Override
			public boolean dispatchKeyEvent(KeyEvent e) {
				if (e.getID() == KeyEvent.KEY_PRESSED) {
					int keyCode = e.getKeyCode();
					if ((keyCode == 37 || keyCode == 39)) {
						t += (keyCode - 38) * (speed);
						setTitle("time = " + t);
						draw();
					}

					if (keyCode == 'R') {
						t = 0;
						setTitle("time = " + t);
						draw();
					}

					if (keyCode == 38 || keyCode == 40) {
						logSpeed = max(0, min(10, logSpeed - keyCode + 39));
						speed = 1 << logSpeed;
						setTitle("speed = " + speed);
					}

				}
				return false;
			}
		});

		add(label);
		setBounds(320, 240, 640, 480);
	}
}
