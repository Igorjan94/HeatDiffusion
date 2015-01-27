import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

public class PhysicalValue implements Serializable {
	static final long serialVersionUID = 0xE1A;
	public final String name, unit;
	public final double value;

	public PhysicalValue(String name) {
		this(name, 0.0, "");
	}

	public PhysicalValue(String name, String unit) {
		this(name, 0.0, unit);
	}

	public PhysicalValue(String name, double value) {
		this(name, value, "");
	}

	public PhysicalValue(String name, double value, String unit) {
		if (name == null) {
			throw new NullPointerException("Name can't be null.");
		}
		if (unit == null) {
			throw new NullPointerException("Unit can't be null.");
		}
		this.name = name;
		this.unit = unit;
		this.value = value;
	}

	public int intValue() {
		return (int) Math.round(value);
	}

	public String toString() {
		if (unit.equals("")) {
			return name + " = " + value;
		} else {
			return name + " = " + value + " " + unit;
		}
	}

}
