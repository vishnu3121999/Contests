import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

public class BigDecimal_ {
    private static final BigDecimal TWO = BigDecimal.valueOf(2);
    private static final BigDecimal ONE = BigDecimal.ONE;
    private static final MathContext mc = new MathContext(20, RoundingMode.HALF_UP); // Precision


    public static void main(String[] args) {
        // Creating BigDecimal from a string (preferred)
        BigDecimal num1 = new BigDecimal("975461057789971042");
        BigDecimal num2 = new BigDecimal("987654321");

        // Addition
        BigDecimal sum = num1.add(num2);
        System.out.println("Sum: " + sum);

        // Subtraction
        BigDecimal difference = num1.subtract(num2);
        System.out.println("Difference: " + difference);

        // Multiplication
        BigDecimal product = num1.multiply(num2);
        System.out.println("Product: " + product);

        // Division with scale and rounding mode
        BigDecimal division = num1.divide(num2, 10, RoundingMode.HALF_UP); // 10 decimal places
        System.out.println("Division: " + division);

        // Square root approximation using Math.sqrt on BigDecimal (not natively supported in BigDecimal)
        BigDecimal sqrt = sqrt(num1, 20); // Precision up to 10 decimal places
        System.out.println("Square Root: " + sqrt);
        System.out.println("Square Root2: " + num1.sqrt(mc).setScale(0,RoundingMode.CEILING));

        BigDecimal number = new BigDecimal(Math.pow(2,47)+"");
        BigDecimal logValue = log(number);
        System.out.println("Natural log (ln): " + logValue);

        // To calculate log base 2
        BigDecimal LN2 = log(new BigDecimal("2"));
        BigDecimal log10Value = logValue.divide(LN2, mc); // log10(x) = ln(x) / ln(10)
        System.out.println("Log base 10: " + log10Value);
    }

    // Helper function to calculate square root of a BigDecimal
    public static BigDecimal sqrt(BigDecimal value, int scale) {
        BigDecimal two = BigDecimal.valueOf(2);
        BigDecimal x = new BigDecimal(Math.sqrt(value.doubleValue()));
        BigDecimal lastX = BigDecimal.ZERO;

        // Newton's method
        while (!x.equals(lastX)) {
            lastX = x;
            x = value.divide(x, scale, RoundingMode.HALF_UP);
            x = x.add(lastX);
            x = x.divide(two, scale, RoundingMode.HALF_UP);
        }
        return x;
    }

    // Natural logarithm (ln) using an approximation
    public static BigDecimal log(BigDecimal value) {
        if (value.compareTo(BigDecimal.ZERO) <= 0) {
            throw new IllegalArgumentException("Logarithm is undefined for zero or negative numbers.");
        }

        int scale = value.scale();
        int precision = 20 + scale;  // More precision for more accurate results
        MathContext mc = new MathContext(precision);

        BigDecimal x = value;
        int k = 0;

        // Normalize x to range [0.5, 1) by factoring out powers of 2
        while (x.compareTo(BigDecimal.ONE) > 0) {
            x = x.divide(TWO, mc);
            k++;
        }
        while (x.compareTo(BigDecimal.ONE) < 0) {
            x = x.multiply(TWO, mc);
            k--;
        }

        // Use log(x) = log(1 + y) with y = x - 1 and iterate with Taylor series
        BigDecimal y = x.subtract(ONE);
        BigDecimal result = BigDecimal.ZERO;
        BigDecimal term = y;

        // Apply the Taylor series: ln(1 + y) = y - y^2/2 + y^3/3 - y^4/4 + ...
        for (int n = 1; n < precision; n++) {
            if (n % 2 == 1) {  // Odd term
                result = result.add(term.divide(BigDecimal.valueOf(n), mc));
            } else {  // Even term (negative)
                result = result.subtract(term.divide(BigDecimal.valueOf(n), mc));
            }
            term = term.multiply(y, mc);  // Update term to y^n
        }

        // Adjust for the factor of 2^k that we factored out initially
        return result.add(BigDecimal.valueOf(k).multiply(BigDecimal.valueOf(Math.log(2)), mc));
    }

}
