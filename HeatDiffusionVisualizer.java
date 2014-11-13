import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Graphics;
import java.awt.KeyEventDispatcher;
import java.awt.KeyboardFocusManager;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.WindowConstants;

import static java.lang.Math.*;

public class HeatDiffusionVisualizer extends JFrame {


	BufferedImage canvas = new BufferedImage(42, 23, BufferedImage.TYPE_INT_RGB);
	JLabel graph = new JLabel();

    static int n;
    static double D, K, E, R, Q, T_0, T_m, C, rho, lambda, alpha, beta, gamma, length, sigma_h, sigma_r, u;
	int s = 0, textH = 23;

	String[] varName = new String[] { "n(c)", "K(1/c)", "E(Дж/моль)", "R(Дж/моль/К)", "Q(Дж/кг)", "ρ(кг/м3)", "T_0(К)", "C(Дж/кг/К)", "λ(Дж/м/с/К)" };
	JTextField[] var = new JTextField[varName.length];

	Component[][] component = new Component[2][];

	double[][][] f = new double[0][0][0];

	Random rnd = new Random();
	int t = 0, maxTime = 0;

	int speed = 1, logSpeed = 0;

	double[] curve;

	public void onResize() {
		int ch = max(getHeight() - textH * (component.length + 2), textH);
		int cw = max(getWidth(), textH);
		curve = new double[cw];
		canvas = new BufferedImage(cw, ch, BufferedImage.TYPE_INT_RGB);
		graph.setBounds(0, 0, cw, ch);
		graph.setIcon(new ImageIcon(canvas));

		for (int i = 0; i < component.length; i++) {
			for (int j = 0; j < component[i].length; j++) {
				int w = max(textH, cw / component[i].length);
				component[i][j].setBounds(j * w, ch + i * textH, w, textH);
			}
		}
		draw();
	}

    public static double R(double x, double T_0, double T_m)
    {
        double c = (x - T_0) / (T_m - T_0) - 0.06;
        if (c < 0)
            c += 1.0;
        return c;
    }

	public void draw() {
		int w = canvas.getWidth(), h = canvas.getHeight();
		for (int x = 0; x < w; x++) {
			for (int y = 0; y < h; y++) {
				canvas.setRGB(x, y, Color.HSBtoRGB((float) R(h - y, 0, h), 1.0f, 1.0f));
			}
		}

		repaint();
	}

	void reCalc() {
		try {

			n       = Integer.parseInt  (var[0].getText());
			K       = Double.parseDouble(var[1].getText());
			E       = Double.parseDouble(var[2].getText());
			R       = Double.parseDouble(var[3].getText());
			Q       = Double.parseDouble(var[4].getText());
			rho     = Double.parseDouble(var[5].getText());
			T_0     = Double.parseDouble(var[6].getText());
			C       = Double.parseDouble(var[7].getText());
			lambda  = Double.parseDouble(var[8].getText());
            T_m     = T_0 + Q / C;
            alpha   = 1;
            beta    = R * T_m / E;
            gamma   = R * T_m * T_m * C / E / Q;
            u       = sqrt(2 * K * lambda * R * T_m * T_m * Math.exp(-E / R / T_m) / Q / rho / E);
            sigma_h = lambda / rho / C / u;
            sigma_r = sigma_h * beta;
            length  = sigma_h * 1000;

			if (2 > n || n < 16) {
				throw new Exception("2 > n || n < 16");
			}
			maxTime = n;

			int curW = canvas.getWidth();
			double curH = canvas.getHeight();
			double[] l = new double[curW];

			Entry<Integer, Integer> last = null;

			t = 0;
			draw();
			setTitle("Recalc " + n);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public HeatDiffusionVisualizer() {
		this.setLayout(null);

		Container content = getContentPane();
		content.add(graph);
        component[0] = new Component[varName.length];

        for (int i = 0; i < varName.length; i++) {
            JTextField name = new JTextField();
            name.setText(varName[i]);
            name.setEditable(false);
            content.add(component[0][i] = name);
        }
        component[1] = new Component[var.length];
        var[0] = new JTextField("1000");
        var[1] = new JTextField("1600000");
        var[2] = new JTextField("80000");
        var[3] = new JTextField("8.314");
        var[4] = new JTextField("700000");
        var[5] = new JTextField("830");
        var[6] = new JTextField("293");
        var[7] = new JTextField("1980");
        var[8] = new JTextField("0.13");
        for (int i = 0; i < varName.length; i++) {
            content.add(component[1][i] = var[i]);
        }

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
					if ((keyCode == 37 || keyCode == 39) && 0 < maxTime) {
						t = (t + maxTime * 4 + (keyCode - 38) * (speed % maxTime)) % maxTime;
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
		this.setSize(640, 480);
	}

	static final long serialVersionUID = 0xE1A;

	public static void main(String[] args) {
		HeatDiffusionVisualizer hdv = new HeatDiffusionVisualizer();
		hdv.setVisible(true);
		hdv.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
	}
}
