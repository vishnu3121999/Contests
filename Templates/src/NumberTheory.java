import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class NumberTheory {

    public static long gcd(long a, long b)
    {
        if(a == 0L)
            return b;
        return gcd(b%a, a);
    }

    public static long lcm(long a, long b) {
        return (a / gcd(a, b)) * b; // To prevent overflow
    }
    // Returns factorial of n
    static long fact(int n)
    {
        if(n==0 || n==1)
            return 1;
        long res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i;
        return res;
    }

    static long factMod(int n, int mod)
    {
        if(n==0 || n==1)
            return 1;
        long res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i%mod;
        return res;
    }
    // Function to calculate the power of a number (a) raised to the power of b modulo mod
    public static long power(long a, long b, long mod) {
        long result = 1;
        a = a %mod;
        while (b > 0) {
            // If the current bit of b is set, multiply the result by a
            if ((b & 1) == 1)
                result = (result * a) % mod;

            // Square the value of a and reduce it modulo mod
            a = (a * a) % mod;

            // Right shift b to move to the next bit
            b >>= 1;
        }
        return result;
    }

    public static ArrayList<Integer> getPrimeFactors(int n)
    {
        ArrayList<Integer> pf = new ArrayList<>();
        // Print the number of 2s that divide n
        while (n%2==0)
        {
            pf.add(2);
            n /= 2;
        }
        // n must be odd at this point.  So we can skip one element (Note i = i +2)
        for (int i = 3; i <= Math.sqrt(n); i+= 2)
        {
            // While i divides n, print i and divide n
            while (n%i == 0)
            {
                pf.add(i);
                n /= i;
            }
        }
        // This condition is to handle the case when n is a prime number greater than 2
        if (n > 2)
            pf.add(n);
        return pf;
    }

    static int phi(int n)
    {
        float result = n;
        for (int p = 2; p * p <= n; ++p) {
            if (n % p == 0) {
                while (n % p == 0)
                    n /= p;
                result *= (1.0 - (1.0 / (float)p));
            }
        }
        if (n > 1)
            result -= result / n;
        return (int)result;
    }

    static ArrayList<Integer> getFactors(int n){
        ArrayList<Integer> l = new ArrayList<>();
        for (int i = 1; i*i <=n ; i++) {
            if(n%i==0){
                l.add(i);
                if(i!=n/i){
                    l.add(n/i);
                }
            }
        }
        return l;
    }

    public static boolean isPrime(long n) {
        if (n <= 1) return false;
        if (n <= 3) return true;
        if (n % 2 == 0 || n % 3 == 0) return false;
        for (long i = 5; i * i <= n; i += 6) {
            if (n % i == 0 || n % (i + 2) == 0) return false;
        }
        return true;
    }

    // Modular inverse using Fermat's Little Theorem (when mod is prime)
    public static long modularInverse(long a, long mod) {
        return power(a, mod - 2, mod);
    }

    // nCr % mod (precomputed factorials and modular inverses required)
    public static long nCr(long n, long r, long mod) {
        if (r > n) return 0;
        long[] fact = new long[(int) (n + 1)];
        fact[0] = 1;
        for (int i = 1; i <= n; i++) fact[i] = (fact[i - 1] * i) % mod;

        return (fact[(int) n] * modularInverse(fact[(int) r], mod) % mod
                * modularInverse(fact[(int) (n - r)], mod) % mod) % mod;
    }




    // Method to calculate the GCD of all elements in a long array
    public static long gcd(long[] arr) {
        long result = arr[0];
        for (int i = 1; i < arr.length; i++) {
            result = gcd(result, arr[i]);
            if (result == 1) break; // Early stop if GCD is 1
        }
        return result;
    }

    // Method to calculate the GCD of all elements in an ArrayList<Long>
    public static long gcd(ArrayList<Long> list) {
        long result = list.get(0);
        for (int i = 1; i < list.size(); i++) {
            result = gcd(result, list.get(i));
            if (result == 1) break; // Early stop if GCD is 1
        }
        return result;
    }

    public static long lcm(long[] arr) {
        long result = arr[0];
        for (int i = 1; i < arr.length; i++) {
            result = lcm(result, arr[i]);
            if (result == 0) break; // Early exit if LCM is zero
        }
        return result;
    }

    // Sieve of Eratosthenes for generating primes up to N
    public static ArrayList<Integer> sieve(int n) {
        boolean[] isPrime = new boolean[n + 1];
        Arrays.fill(isPrime, true);
        isPrime[0] = isPrime[1] = false;
        for (int i = 2; i * i <= n; i++) {
            if (isPrime[i]) {
                for (int j = i * i; j <= n; j += i) isPrime[j] = false;
            }
        }
        ArrayList<Integer> primes = new ArrayList<>();
        for (int i = 2; i <= n; i++) if (isPrime[i]) primes.add(i);
        return primes;
    }

    public static ArrayList<Integer> getDigits(long n) {
        ArrayList<Integer> digits = new ArrayList<>();
        while (n > 0) {
            digits.add((int) (n % 10));
            n /= 10;
        }
        Collections.reverse(digits); // If you need most significant to least significant
        return digits;
    }
}
