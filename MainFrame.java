import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JTextField;
import javax.swing.WindowConstants;

public class MainFrame extends JFrame {
    static final long serialVersionUID = 0xE1A;

    int border = 10;

    public MainFrame() {
        setLayout(null);

        final PhysicalValue[] pv = ReactionSolover.getVariables();
        final int n = pv.length;

        final JTextField[][] var = new JTextField[n][3];
        final JButton button = new JButton("GO");
        add(button);

        button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                PhysicalValue[] v = new PhysicalValue[n];
                for (int i = 0; i < n; i++) {
                    v[i] = new PhysicalValue(pv[i].name, Double.parseDouble(var[i][1].getText()), pv[i].unit);
                }
                Visualizer visualizer = new Visualizer(ReactionSolover.solove(v));
                visualizer.setVisible(true);
            }
        });

        for (int i = 0; i < n; i++) {
            var[i][0] = new JTextField(pv[i].name);
            var[i][0].setEditable(false);
            add(var[i][0]);

            var[i][1] = new JTextField(Double.toString(pv[i].value));
            add(var[i][1]);

            var[i][2] = new JTextField(pv[i].unit);
            var[i][2].setEditable(false);
            add(var[i][2]);
        }

        addComponentListener(new ComponentAdapter() {
            @Override
            public void componentResized(ComponentEvent e) {
                int w = Math.max(80, getWidth()) - 14, h = Math.max(60, getHeight()) - 37;
                int dw = Math.max(border * 2, w / 3 - border);
                int dh = Math.max(border * 2, h / (n + 1) - border);
                {
                    int curY = border / 2;
                    for (int i = 0; i < n; i++) {
                        int curX = border / 2;
                        for (int j = 0; j < 3; j++) {
                            var[i][j].setBounds(curX, curY, dw, dh);
                            curX += border + dw;
                        }
                        curY += border + dh;
                    }
                    button.setBounds(border / 2, curY, w - border, dh);
                }
            }
        });

        setBounds(320, 240, 700, 240);
    }

    public static void main(String[] args) {
        MainFrame frame = new MainFrame();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    }

}
