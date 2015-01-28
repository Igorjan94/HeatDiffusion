import static java.lang.Math.*;

import java.awt.Color;
public class ReactionSolover {
	
	public static PhysicalValue[] getVariables() {
		PhysicalValue[] physicalValue = new PhysicalValue[] {
				new PhysicalValue("a", 0.007),
				new PhysicalValue("b", 2),
				new PhysicalValue("c", 3),			
				new PhysicalValue("time", 1000, "sec"),
				new PhysicalValue("length", 1000, "px"),
		};
		return physicalValue;
	}

	public static Graph[] solove(PhysicalValue[] physicalValue) {
		double a = physicalValue[0].value;
		double b = physicalValue[1].value;
		double c = physicalValue[2].value;
		int    n = physicalValue[3].intValue();
		int    m = physicalValue[4].intValue();		
		
		double[][][] d = new double[3][n][m];		
				
		for (int t = 0; t < n; t++) {
			for (int x = 0; x < m; x++) {
				d[0][t][x] = sin(a * x + b * t);
				d[1][t][x] = cos(a * x + b * t + c * x * t);
				d[2][t][x] = d[0][t][x] + d[1][t][x];
			}
		}
		
		Graph[] graph = new Graph[] {  
				new Graph(d[0], "sin"),
				new Graph(d[1], "cos"),
				new Graph(d[2], "sin + cos"),			
		};		
		
		return graph;
	}

}
