import java.io.*;
import java.util.*;
import static java.lang.Math.*;
import java.util.function.*;
import java.util.stream.Collectors;



public class Main implements Runnable {
    int mod = 1000000007;
    int mod2 = 998244353;



    FastReader fr = new FastReader(System.in);
    public static void main(String[] args) throws FileNotFoundException {
        new Thread(null, new Main(), "Main", 1 << 28).start();
    }

    @Override
    public void run() {


        int T = fr.nextInt();
        while(T-->0) {

            int n = fr.nextInt();
            print(2*n-1);
            for (int i = 1; i <= n-1; i++) {
                print(i,1,i);
                print(i,i+1,n);
            }
            print(n,1,n);

        }
    }

    private boolean isPossibel(int[] a, int[] b, int x, int mid, int n) {

        int idx = n;
        int k = mid;
        for (int i = n - 1; i >= 0; i--) {
            if(k==0)break;
            if(k>0 && x-k+1>=b[i]){
                k--;
            }
            if(k==0){
                idx=i;break;
            }
        }
//        print(mid,idx);
        for (int i = 0; i < idx; i++) {
            if(x<a[i])return false;

        }
        for (int i = idx; i <n ; i++) {
            if(x<a[i])return false;
            if(x>=b[i])
                x--;
        }

        return true;

    }

    private int cal(int prev, int i) {

        String p = Integer.toBinaryString(prev);
        String s = new StringBuilder(Integer.toBinaryString(i)).toString();
        s= "0".repeat(p.length()-s.length())+s;
//        print(prev,i,s);
        boolean x= false;
        for (int j = 0; j <=s.length()-1; j++) {
            if(s.charAt(j)==p.charAt(j) && !x)continue;
            else x=true;
            if(s.charAt(j)=='0'){
                int or = (1<<(s.length()-1-j))|i;
//                print(s,j,or);
                if(or>=prev){
                    return or;
                }
                else i=or;
            }
        }
        return 0;

    }


    //b  --> searchInsert or Ceil
    // if target exists return its index else return the index of element greater than target
    public int ceil(long[] nums, long target) {
        int l=0,h=nums.length;
        while(l<h){
            int mid = l+(h-l)/2;
            if(nums[mid]>=target)h=mid;
            else l=mid+1;
        }
        return l;
    }
    // c
    // if target exists return its index else return the index of element smaller than target
    public int floor(long[] nums, long target) {
        int l=0,h=nums.length;
        while(l<h){
            int mid = l+(h-l)/2;
            if(nums[mid]>target)h=mid;
            else l=mid+1;
        }
        return l-1;
    }

    public static String sortString(String input) {
        char[] chars = input.toCharArray();
        Arrays.sort(chars);
        return new String(chars);
    }

    public LinkedHashMap<Integer, Integer> sortByValue(HashMap<Integer, Integer> hm, int comparator_0ForASC_1forDESC) {
        Comparator<Map.Entry<Integer, Integer>> comparator =
                comparator_0ForASC_1forDESC == 0
                        ? Comparator.comparing((Map.Entry<Integer, Integer> e) -> e.getValue())
                        : Comparator.comparing((Map.Entry<Integer, Integer> e) -> e.getValue()).reversed();

        return hm.entrySet()
                .stream()
                .sorted(comparator)
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        Map.Entry::getValue,
                        (e1, e2) -> e1,
                        LinkedHashMap::new));
    }

    public static void push(Map<Integer, Integer> map, Integer k, Integer v)
    {
        //map[k] += v;
        if(!map.containsKey(k))
            map.put(k, v);
        else
            map.put(k, map.get(k)+v);
    }
    public static void pull(Map<Integer, Integer> map, int k, int v)
    {
        //assumes map[k] >= v
        //map[k] -= v
        int lol = map.get(k);
        if(lol == v)
            map.remove(k);
        else
            map.put(k, lol-v);
    }
    public static void push(Map<Long, Long> map, long k, long v)
    {
        //map[k] += v;
        if(!map.containsKey(k))
            map.put(k, v);
        else
            map.put(k, map.get(k)+v);
    }
    public static void pull(Map<Long, Long> map, long k, long v)
    {
        //assumes map[k] >= v
        //map[k] -= v
        long lol = map.get(k);
        if(lol == v)
            map.remove(k);
        else
            map.put(k, lol-v);
    }
    public static void push(Map<String, Integer> map, String k, Integer v)
    {
        //map[k] += v;
        if(!map.containsKey(k))
            map.put(k, v);
        else
            map.put(k, map.get(k)+v);
    }
    public static void push(Map<Character, Integer> map, char k, int v) {
        if (!map.containsKey(k))
            map.put(k, v);
        else
            map.put(k, map.get(k) + v);
    }

    public static void pull(Map<Character, Integer> map, char k, int v) {
        // Assumes map.get(k) >= v
        int current = map.get(k);
        if (current == v)
            map.remove(k);
        else
            map.put(k, current - v);
    }

    public static int ceil(int n,int m){
        return (int) Math.ceil(n/(double)m);
    }
    public static long ceil(long n,long m){
        return (long) Math.ceil(n/(double)m);
    }

    static void print(Object... args) {
        for (int i = 0; i < args.length; i++) {
            Object o = args[i];
            if (o.getClass().isArray()) {
                if (o instanceof int[])
                    System.out.println(Arrays.toString((int[]) o));
                else if (o instanceof long[])
                    System.out.println(Arrays.toString((long[]) o));
                else if (o instanceof double[])
                    System.out.println(Arrays.toString((double[]) o));
                else if (o instanceof float[])
                    System.out.println(Arrays.toString((float[]) o));
                else if (o instanceof boolean[])
                    System.out.println(Arrays.toString((boolean[]) o));
                else if (o instanceof char[])
                    System.out.println(Arrays.toString((char[]) o));
            }
            else {
                System.out.print(o);
                if (i != args.length - 1) {
                    System.out.print(" ");
                }
            }
        }
        System.out.println();
    }

    public static void print(int[] arr) {
        StringBuilder sb = new StringBuilder();
        for (int val : arr) sb.append(val).append(" ");
        System.out.println(sb);
    }
    public static void print(long[] arr) {
        StringBuilder sb = new StringBuilder();
        for (long val : arr) sb.append(val).append(" ");
        System.out.println(sb);
    }
    public static <T> void print(ArrayList<T> list) {
        StringBuilder sb = new StringBuilder();
        for (T val : list) sb.append(val).append(" ");
        System.out.println(sb);
    }

    public int[] readArr(int N) {
        int[] arr = new int[N];
        for (int i = 0; i < N; i++)
            arr[i] = fr.nextInt();
        return arr;
    }
    public int[][] read2DArr(int N,int M) {
        int[][] arr = new int[N][M];
        for (int i = 0; i < N; i++){
            for (int j = 0; j < M; j++) {
                arr[i][j]=fr.nextInt();
            }
        }
        return arr;
    }
    public long[] readlongArr(int N) {
        long[] arr = new long[N];
        for (int i = 0; i < N; i++)
            arr[i] = fr.nextLong();
        return arr;
    }


    static class FastReader {
        BufferedReader br;
        StringTokenizer st;

        public FastReader(InputStream inputStream)
        {
            br = new BufferedReader(
                    new InputStreamReader(inputStream));
        }

        String next()
        {
            while (st == null || !st.hasMoreElements()) {
                try {
                    st = new StringTokenizer(br.readLine());
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            }
            return st.nextToken();
        }

        int nextInt() { return Integer.parseInt(next()); }

        long nextLong() { return Long.parseLong(next()); }

        double nextDouble()
        {
            return Double.parseDouble(next());
        }

        String nextLine()
        {
            String str = "";
            try {
                if(st.hasMoreTokens()){
                    str = st.nextToken("\n");
                }
                else{
                    str = br.readLine();
                }
            }
            catch (IOException e) {
                e.printStackTrace();
            }
            return str;
        }
    }
}


class Pair<T, U> {
    T key;
    U value;

    Pair(T key, U value) {
        this.key = key;
        this.value = value;
    }

    @Override
    public String toString() {
        return "Pair{" +
                "key=" + key +
                ", value=" + value +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Pair<?, ?> pair = (Pair<?, ?>) o;
        return Objects.equals(key, pair.key) && Objects.equals(value, pair.value);
    }

    @Override
    public int hashCode() {
        return Objects.hash(key, value);
    }
}

// ArrayUtils
class AU {

    public static List<Integer> toList(int[] arr) {
        return Arrays.stream(arr)
                .boxed()    // Convert each int to an Integer
                .collect(Collectors.toList());
    }

    public static List<Long> toList(long[] arr) {
        return Arrays.stream(arr)
                .boxed()    // Convert each long to a Long
                .collect(Collectors.toList());
    }

    public static void reverse(int[] arr) {
        int left = 0, right = arr.length - 1;
        while (left < right) {
            int temp = arr[left];
            arr[left] = arr[right];
            arr[right] = temp;
            left++;
            right--;
        }
    }

    public static void swap(int[] a, int i, int j) {
        int t = a[i]; a[i] = a[j]; a[j] = t;
    }

    // Method to find the minimum value in a long array
    public static long min(long[] arr) {
        long minVal = Long.MAX_VALUE;
        for (long val : arr) {
            if (val < minVal) {
                minVal = val;
            }
        }
        return minVal;
    }

    // Method to find the minimum value in an ArrayList<Long>
    public static long min(ArrayList<Long> list) {
        return Collections.min(list);
    }


    // Method to find the maximum value in a long array
    public static long max(long[] arr) {
        long maxVal = Long.MIN_VALUE;
        for (long val : arr) {
            if (val > maxVal) {
                maxVal = val;
            }
        }
        return maxVal;
    }

    // Method to find the maximum value in an ArrayList<Long>
    public static long max(ArrayList<Long> list) {
        return Collections.max(list);
    }


    // Method to find the index of the minimum value in a long array
    public static int minI(long[] arr) {
        int index = 0;
        long minVal = Long.MAX_VALUE;
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] < minVal) {
                minVal = arr[i];
                index = i;
            }
        }
        return index;
    }

    // Method to find the index of the maximum value in a long array
    public static int maxI(long[] arr) {
        int index = 0;
        long maxVal = Long.MIN_VALUE;
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] > maxVal) {
                maxVal = arr[i];
                index = i;
            }
        }
        return index;
    }


    // Method to find the minimum value in an int array
    public static int min(int[] arr) {
        int minVal = Integer.MAX_VALUE;
        for (int val : arr) {
            if (val < minVal) {
                minVal = val;
            }
        }
        return minVal;
    }

    // Method to find the maximum value in an int array
    public static int max(int[] arr) {
        int maxVal = Integer.MIN_VALUE;
        for (int val : arr) {
            if (val > maxVal) {
                maxVal = val;
            }
        }
        return maxVal;
    }

    // Method to find the index of the minimum value in an int array
    public static int minI(int[] arr) {
        int index = 0;
        int minVal = Integer.MAX_VALUE;
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] < minVal) {
                minVal = arr[i];
                index = i;
            }
        }
        return index;
    }

    // Method to find the index of the maximum value in an int array
    public static int maxI(int[] arr) {
        int index = 0;
        int maxVal = Integer.MIN_VALUE;
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] > maxVal) {
                maxVal = arr[i];
                index = i;
            }
        }
        return index;
    }

    public static long sum(long[] arr) {
        long sum = 0;
        for (long val : arr) {
            sum += val;
        }
        return sum;
    }
    public static long sum(ArrayList<Long> list) {
        long sum = 0;
        for (long val : list) {
            sum += val;
        }
        return sum;
    }
    public static long sum(int[] arr) {
        long sum = 0;
        for (int val : arr) {
            sum += val;
        }
        return sum;
    }

    public static long[] prefixSum(int[] arr) {
        int n = arr.length;
        // Create an array with one extra element.
        long[] ps = new long[n + 1];
        ps[0] = 0; // This dummy element represents the sum before any elements.
        for (int i = 0; i < n; i++) {
            ps[i + 1] = ps[i] + arr[i];
        }
        return ps;
    }
    public static long[] suffixSum(int[] arr) {
        int n = arr.length;
        // Create an array with one extra element.
        long[] ss = new long[n + 1];
        ss[n] = 0; // This dummy element represents the sum after all elements.
        for (int i = n - 1; i >= 0; i--) {
            ss[i] = ss[i + 1] + arr[i];
        }
        return ss;
    }


    public static void sort(int[] arr)
    {
        //because Arrays.sort() uses quicksort which is O(n^2) worstcase
        //Collections.sort() uses merge sort
        ArrayList<Integer> ls = new ArrayList<Integer>();
        for(int x: arr)
            ls.add(x);
        Collections.sort(ls);
        for(int i=0; i < arr.length; i++)
            arr[i] = ls.get(i);
    }


}

//NumberTheory
class NT {

    public static long gcd(long a, long b)
    {
        if(a == 0L)
            return b;
        return gcd(b%a, a);
    }

    public static long lcm(long a, long b) {
        return (a / gcd(a, b)) * b; // To prevent overflow
    }

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

    // O(Sqrt(N)) , for arrays O(N*Sqrt(N))
    public static ArrayList<Long> getPrimeFactors(long n)
    {
        ArrayList<Long> pf = new ArrayList<>();
        // Print the number of 2s that divide n
        while (n%2==0)
        {
            pf.add(2L);
            n /= 2;
        }
        // n must be odd at this point.  So we can skip one element (Note i = i +2)
        for (long i = 3; i <= Math.sqrt(n); i+= 2)
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

    // for arrays O(N*Sqrt(N)/logN)
    public static List<Long> getPrimeFactors(long n, List<Integer> primesUntilSqRootN) {
        List<Long> pf = new ArrayList<>();
        long tmp = n;
        for (int p : primesUntilSqRootN) {
            long pp = (long)p * p;
            if (pp > tmp) break;
            while (tmp % p == 0) {
                pf.add((long)p);
                tmp /= p;
            }
        }
        if (tmp > 1) {
            // tmp is now either 1 or a prime > sqrt(n)
            pf.add(tmp);
        }
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

    static ArrayList<Long> getFactors(long n){
        ArrayList<Long> l = new ArrayList<>();
        for (long i = 1; i*i <=n ; i++) {
            if(n%i==0){
                l.add(i);
                if(i!=n/i){
                    l.add(n/i);
                }
            }
        }
        Collections.sort(l);
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

    public static long modDiv(long a, long b, long mod) {
        // bring a into [0,mod)
        long am = ((a % mod) + mod) % mod;
        // compute b^{-1} mod
        long invB = modularInverse(b, mod);
        return (am * invB) % mod;
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
    public static int gcd(ArrayList<Integer> list) {
        int result = list.get(0);
        for (int i = 1; i < list.size(); i++) {
            result = (int) gcd(result, list.get(i));
            if (result == 1) break; // Early stop if GCD is 1
        }
        return result;
    }
    public static long gcd2(ArrayList<Long> list) {
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

    public class ExtendedEuclid {
        // Helper class to store results
        static class Result {
            int gcd;  // GCD of a and b
            int x;    // Coefficient for a
            int y;    // Coefficient for b

            Result(int gcd, int x, int y) {
                this.gcd = gcd;
                this.x = x;
                this.y = y;
            }
        }

        // Extended Euclidean Algorithm
        public static Result extendedEuclid(int a, int b) {
            if (b == 0) {
                // Base case: gcd(a, 0) = a, x = 1, y = 0
                return new Result(a, 1, 0);
            }

            // Recursive case
            Result next = extendedEuclid(b, a % b);

            // Update x and y using the results of recursion
            int x = next.y;
            int y = next.x - (a / b) * next.y;

            return new Result(next.gcd, x, y);
        }

        public static boolean findAnySolution(int a, int b, int c, int[] solution) {
            Result result = extendedEuclid(Math.abs(a), Math.abs(b));
            int g = result.gcd;

            // No solution if c is not divisible by gcd(a, b)
            if (c % g != 0) {
                return false;
            }

            // Scale solution
            int x0 = result.x * (c / g);
            int y0 = result.y * (c / g);

            // Adjust signs if a or b is negative
            if (a < 0) x0 = -x0;
            if (b < 0) y0 = -y0;

            // Store the solution
            solution[0] = x0;
            solution[1] = y0;

            return true;
        }

        public static List<int[]> findAllSolutions(int a, int b, int c,
                                                   int xmin, int xmax,
                                                   int ymin, int ymax) {
            List<int[]> solutions = new ArrayList<>();

            // Find any solution
            Result result = extendedEuclid(Math.abs(a), Math.abs(b));
            int g = result.gcd;

            // No solution if c is not divisible by gcd(a, b)
            if (c % g != 0) {
                return solutions; // Empty list
            }

            // Scale solution
            int x0 = result.x * (c / g);
            int y0 = result.y * (c / g);

            // Adjust signs if a or b is negative
            if (a < 0) x0 = -x0;
            if (b < 0) y0 = -y0;

            // Step sizes
            int stepX = b / g;
            int stepY = a / g;

            // Calculate range for k based on x bounds
            int kxMin = (int) Math.ceil((xmin - x0) / (double) stepX);
            int kxMax = (int) Math.floor((xmax - x0) / (double) stepX);

            // Calculate range for k based on y bounds
            int kyMin = (int) Math.ceil((y0 - ymax) / (double) stepY);
            int kyMax = (int) Math.floor((y0 - ymin) / (double) stepY);

            // Find overlapping range for k
            int kMin = Math.max(kxMin, kyMin);
            int kMax = Math.min(kxMax, kyMax);

            // Generate all solutions for valid k
            for (int k = kMin; k <= kMax; k++) {
                int x = x0 + k * stepX;
                int y = y0 - k * stepY;
                solutions.add(new int[]{x, y});
            }

            return solutions;
        }

    }

    public static String getBitMask(int n) {
        String s = Integer.toBinaryString(n);
        return "0".repeat(32-s.length())+s;
    }

    public static int findMex(int[] arr) {
        int n = arr.length;
        boolean[] present = new boolean[n + 1];
        for (int x : arr) {
            if (x >= 0 && x <= n) {
                present[x] = true;
            }
        }
        for (int i = 0; i <= n; i++) {
            if (!present[i]) {
                return i;
            }
        }
        return n + 1;
    }

    public static boolean isPerfectSquare(long num) {
        if (num < 0) return false;
        long sqrt = (long) Math.sqrt(num);
        return sqrt * sqrt == num;
    }

    // Math.sqrt is O(1) but it has floating point issues for big numbers
    // Math.sqrt(big_non-perfect_square number) maybe be rounded to long with 0 as decimal part.
    // ex : sqrt(975461057789971042) = 987654321.0 (but its actually 987654321.00000000050625)
    // Below method returns
    //          -1, if integer sqrt doesnt exist
    //          int sqrt, if it exists
    public static long integerSqrt(long x) {
        if (x < 0) {
            return -1;
        }
        if (x < 2) {
            return x;
        }

        long left = 1;
        long right = x / 2 + 1; // An upper bound for sqrt(x) when x > 1

        while (left < right) {
            long mid = left + (right - left + 1) / 2;

            // Instead of mid*mid which might overflow,
            // use mid <= x / mid to compare.
            if (mid <= x / mid) {
                left = mid;
            } else {
                right = mid - 1;
            }
        }
        if (left * left == x) {
            return left;
        } else {
            return -1;
        }
    }

}

// Zero-based Indexing of Vertices
class Graph_ {

    private int V; // Number of vertices
    public List<List<int[]>> adj; // Adjacency list: each int[] contains {neighbor, weight}

    public Graph_(int V) {
        this.V = V;
        adj = new ArrayList<>(V);
        for (int i = 0; i < V; i++) {
            adj.add(new ArrayList<>());
        }
    }

    public void addEdge(int u, int v, int weight) {
        adj.get(u).add(new int[]{v, weight});
    }

    // Print the adjacency list representation of the graph
    public void printGraph() {
        for (int u = 0; u <V; u++) {
            System.out.print("Vertex " + u + ":");
            for (int[] neighbor : adj.get(u)) {
                int v = neighbor[0];
                int weight = neighbor[1];
                System.out.print(" -> " + v + " (weight " + weight + ")");
            }
            System.out.println();
        }
    }

    // In some problesm edges are not always given in parent->child format, they could be child->parent also
    // To directify the graph from root, use below, so that you dont need to include parent or visited arr as
    //  the method parameter for every methods like getSize, dp etc
    public void convertToTree(int root) {
        boolean[] visited = new boolean[V];
        convertToTreeHelper(root, -1, visited);
    }

    // Helper method that performs DFS from 'node', with 'parent' being the parent vertex.
    // It updates the adj list for each node so that it only includes children (excluding the parent).
    private void convertToTreeHelper(int node, int parent, boolean[] visited) {
        visited[node] = true;
        List<int[]> children = new ArrayList<>();
        // Iterate over the current neighbors of node.
        for (int[] edge : adj.get(node)) {
            int neighbor = edge[0];
            // Skip the edge back to the parent.
            if (neighbor == parent) {
                continue;
            }
            children.add(edge);
            if (!visited[neighbor]) {
                convertToTreeHelper(neighbor, node, visited);
            }
        }
        // Update the node's adjacency list to include only its children.
        adj.set(node, children);
    }

    public Graph_ reverseGraph() {
        Graph_ rev = new Graph_(V);
        for (int u = 0; u < V; u++) {
            for (int[] edge : adj.get(u)) {
                int v = edge[0];
                int w = edge[1];
                rev.addEdge(v, u, w);
            }
        }
        return rev;
    }

    // 1. Dijkstra's Algorithm
    // Time Complexity: O((V + E) log V) using a priority queue
    // Space Complexity: O(V) for the distance array and O(V) for the priority queue
    public int[] dijkstra(int src) {
        int[] dist = new int[V];
        Arrays.fill(dist, Integer.MAX_VALUE);
        dist[src] = 0;

        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[1]));
        pq.add(new int[]{src, 0});

        while (!pq.isEmpty()) {
            int[] curr = pq.poll();
            int u = curr[0];
            int d = curr[1];

            if (d > dist[u]) continue;

            for (int[] neighbor : adj.get(u)) {
                int v = neighbor[0];
                int weight = neighbor[1];

                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    pq.add(new int[]{v, dist[v]});
                }
            }
        }

        return dist;
    }

    // 2. Bellman-Ford Algorithm
    // Time Complexity: O(V * E)
    // Space Complexity: O(V) for the distance array
    public int[] bellmanFord(int src) {
        int[] dist = new int[V];
        Arrays.fill(dist, Integer.MAX_VALUE);
        dist[src] = 0;

        for (int i = 1; i < V; i++) { // Relax all edges V-1 times
            for (int u = 0; u < V; u++) {
                for (int[] edge : adj.get(u)) {
                    int v = edge[0];
                    int weight = edge[1];
                    if (dist[u] != Integer.MAX_VALUE && dist[u] + weight < dist[v]) {
                        dist[v] = dist[u] + weight;
                    }
                }
            }
        }

        // Check for negative-weight cycles
        for (int u = 0; u < V; u++) {
            for (int[] edge : adj.get(u)) {
                int v = edge[0];
                int weight = edge[1];
                if (dist[u] != Integer.MAX_VALUE && dist[u] + weight < dist[v]) {
                    System.out.println("Graph contains negative-weight cycle");
                    return null;
                }
            }
        }

        return dist;
    }

    // 3. Floyd-Warshall Algorithm
    // Time Complexity: O(V^3)
    // Space Complexity: O(V^2) for the distance matrix
    public int[][] floydWarshall() {
        int[][] dist = new int[V][V];

        for (int i = 0; i < V; i++) {
            Arrays.fill(dist[i], Integer.MAX_VALUE);
            dist[i][i] = 0;
        }

        for (int u = 0; u < V; u++) {
            for (int[] edge : adj.get(u)) {
                int v = edge[0];
                dist[u][v] = edge[1];
            }
        }

        for (int k = 0; k < V; k++) {
            for (int i = 0; i < V; i++) {
                for (int j = 0; j < V; j++) {
                    if (dist[i][k] != Integer.MAX_VALUE && dist[k][j] != Integer.MAX_VALUE && dist[i][k] + dist[k][j] < dist[i][j]) {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }

        return dist;
    }


    // 4. Prim's Algorithm (for Minimum Spanning Tree)
    // Time Complexity: O((V + E) log V) using a priority queue
    // Space Complexity: O(V + E) for the priority queue and adjacency list
    public int[] prim() {
        int[] dist = new int[V]; // To track the minimum edge weight to reach each vertex
        int[] parent = new int[V]; // To track the MST
        boolean[] inMST = new boolean[V]; // To track included vertices in MST

        Arrays.fill(dist, Integer.MAX_VALUE);
        dist[0] = 0;
        parent[0] = -1; // Start from any vertex (0 in this case)

        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[1]));
        pq.add(new int[]{0, 0});

        while (!pq.isEmpty()) {
            int u = pq.poll()[0];
            inMST[u] = true;

            for (int[] neighbor : adj.get(u)) {
                int v = neighbor[0];
                int weight = neighbor[1];

                if (!inMST[v] && weight < dist[v]) {
                    dist[v] = weight;
                    pq.add(new int[]{v, dist[v]});
                    parent[v] = u;
                }
            }
        }

        return parent; // Returns the MST as an array where parent[i] gives the parent of vertex i in the MST
    }

//    public static void main(String[] args) {
//        Graph_ graph = new Graph_(5);
//
//    }
}

// Zero-based indexing
class DSU_ {
    private int[] parent;
    public int[] size;

    // Time Complexity: O(n)
    // Space Complexity: O(n)
    public DSU_(int n) {
        parent = new int[n];
        size = new int[n];
        for (int i = 0; i < n; i++) {
            parent[i] = i;
            size[i] = 1;  // Initially, each set has size 1
        }
    }

    // Find the root of the set containing x with path compression
    // Time Complexity: O(log n) (amortized)
    public int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);  // Path compression
        }
        return parent[x];
    }

    // Union two sets by size
    // Time Complexity: O(1) (amortized)
    public void union(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX != rootY) {
            // Attach the smaller tree under the larger tree
            if (size[rootX] < size[rootY]) {
                parent[rootX] = rootY;
                size[rootY] += size[rootX];
            } else {
                parent[rootY] = rootX;
                size[rootX] += size[rootY];
            }
        }
    }

    // Check if two elements are in the same set
    public boolean isConnected(int x, int y) {
        return find(x) == find(y);
    }

    public void printDSU() {
        Map<Integer, List<Integer>> groups = new HashMap<>();
        for (int i = 0; i < parent.length; i++) {
            int root = find(i);
            groups.computeIfAbsent(root, k -> new ArrayList<>()).add(i);
        }

        // Separate into multi-element and singletons
        List<Map.Entry<Integer, List<Integer>>> multi = new ArrayList<>();
        List<Integer> singletons = new ArrayList<>();
        for (Map.Entry<Integer, List<Integer>> e : groups.entrySet()) {
            if (e.getValue().size() > 1) {
                multi.add(e);
            } else {
                singletons.add(e.getValue().get(0));
            }
        }

        // Print multi-element sets, if any
        if (!multi.isEmpty()) {
            System.out.println("Multi-element sets:");
            for (Map.Entry<Integer, List<Integer>> e : multi) {
                System.out.println("  Root " + e.getKey() + ": " + e.getValue());
            }
        }

        // Print single-element sets together
        if (!singletons.isEmpty()) {
            System.out.println("Single-element sets: " + singletons);
        }
    }

}

class Recursion {
    private static final ThreadLocal<Integer> depth = ThreadLocal.withInitial(() -> 0);

    private static void increaseDepth() {
        depth.set(depth.get() + 1);
    }

    private static void decreaseDepth() {
        depth.set(depth.get() - 1);
    }

    public static void logEnter(Object... args) {
        System.out.println(indent() + "→ Entering " + formatArgs(args));
        increaseDepth();
    }

    public static void logExit(Object... args) {
        decreaseDepth();
        System.out.println(indent() + "← Exiting " + formatArgs(args));
    }

    public static void logBaseCase(Object result) {
        System.out.println(indent() + "★ Base case reached, returning " + result);
    }
    private static String formatArgs(Object... args) {
        return "(" + String.join(", ", java.util.Arrays.stream(args)
                .map(String::valueOf).toArray(String[]::new)) + ")";
    }

    private static String indent() {
        return "│ ".repeat(depth.get());  // Adds vertical bars for better depth visualization
    }

    // O(2^n · n)
    public static List<List<Integer>> powerSet(int[] a) {
        int n = a.length;
        List<List<Integer>> all = new ArrayList<>(1 << n);
        for (int mask = 0; mask < (1 << n); mask++) {
            List<Integer> subset = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    subset.add(a[i]);
                }
            }
            all.add(subset);
        }
        return all;
    }

    // (in-place)
    public static boolean nextPermutation(int[] a) {
        int n = a.length, i = n - 2;
        while (i >= 0 && a[i] >= a[i + 1]) i--;
        if (i < 0) return false;
        int j = n - 1;
        while (a[j] <= a[i]) j--;
        AU.swap(a, i, j);
        // reverse suffix
        for (int l = i + 1, r = n - 1; l < r; l++, r--) AU.swap(a, l, r);
        return true;
    }

    // ALL PERMUTATIONS using nextPermutation, O(n! · n)
    public static List<List<Integer>> allPermutations(int[] a) {
        List<List<Integer>> all = new ArrayList<>();
        Arrays.sort(a);
        do {
            List<Integer> p = new ArrayList<>();
            for (int x : a) p.add(x);
            all.add(p);
        } while (nextPermutation(a));
        return all;
    }
}

// Zero-based Indexing with L,R inclusive in [L,R]
class SparseTable {
    private final int N;
    private final int K;
    private final int[][] st;       // now [N][K]
    private final int[] log2;
    private final IntBinaryOperator op;

    /**
     * @param arr the input array (length N)
     * @param op  an idempotent operator, e.g. Integer::min, Integer::max, gcd, xor…
     */
    public SparseTable(int[] arr, IntBinaryOperator op) {
        this.N = arr.length;
        this.K = 32 - Integer.numberOfLeadingZeros(N);  // max levels needed
        this.op = op;
        this.st = new int[N][K];
        this.log2 = new int[N + 1];
        buildLog();
        buildST(arr);
    }

    // Precompute floor(log2(i)) for i = 1..N
    private void buildLog() {
        log2[1] = 0;
        for (int i = 2; i <= N; i++) {
            log2[i] = log2[i / 2] + 1;
        }
    }

    // Build in O(N log N)
    private void buildST(int[] arr) {
        // level k = 0: intervals of length 1
        for (int i = 0; i < N; i++) {
            st[i][0] = arr[i];
        }

        // levels k = 1..K-1
        for (int k = 1; k < K; k++) {
            int len = 1 << k;       // interval length = 2^k
            int half = len >> 1;    // half = 2^(k-1)
            for (int i = 0; i + len <= N; i++) {
                // combine [i, i+2^(k-1)-1] and [i+2^(k-1), i+2^k-1]
                st[i][k] = op.applyAsInt(st[i][k - 1], st[i + half][k - 1]);
            }
        }
    }

    /**
     * O(1) query over [L..R] inclusive.
     * @param L left index (0-based)
     * @param R right index (0-based)
     */
    public int query(int L, int R) {
        int len = R - L + 1;
        int k = log2[len];
        // overlap two blocks of length 2^k
        return op.applyAsInt(
                st[L][k],
                st[R - (1 << k) + 1][k]
        );
    }

//    public static void main(String[] args) {
//        int[] a = {5, 2, 4, 7, 1, 3, 6};
//
//        SparseTable stMin = new SparseTable(a, Integer::min);
//        System.out.println(stMin.query(1, 4)); // 1
//
//        SparseTable stMax = new SparseTable(a, Integer::max);
//        System.out.println(stMax.query(2, 5)); // 7
//
//        SparseTable stGcd = new SparseTable(a, (x, y) -> gcd(x, y));
//        System.out.println(stGcd.query(0, 3)); // 1
//    }

}

/**
 * Zero index based
 * point-update and range-query
 * Allowed associative combine function and their identities:
 * sum - 0
 * product - 1
 * min - Integer.MAX_VALUE
 * max - Integer.MIN_VALUE
 * gcd - 0
 * lcm - 1
 * xor - 0
 * & - ?
 * | - ?
 * Time Complexity:
 *  - Build: O(n)
 *  - Update: O(log n)
 *  - Query: O(log n)
 * Space: O(4*n)
 */
class SegmentTree<T> {
    private final int n;
    private final Object[] tree;
    private final BinaryOperator<T> combiner;
    private final T identity;


    public SegmentTree(T[] arr, BinaryOperator<T> combiner, T identity) {
        this.n = arr.length;
        this.combiner = combiner;
        this.identity = identity;
        this.tree = new Object[4 * n];
        build(1, 0, n - 1, arr);
    }

    @SuppressWarnings("unchecked")
    public SegmentTree(int[] arr, BinaryOperator<Integer> combiner, Integer identity) {
        this((T[]) Arrays.stream(arr).boxed().toArray(Integer[]::new), (BinaryOperator<T>) combiner, (T) identity);
    }
    @SuppressWarnings("unchecked")
    public SegmentTree(long[] arr, BinaryOperator<Long> combiner, Long identity) {
        this((T[]) Arrays.stream(arr).boxed().toArray(Long[]::new), (BinaryOperator<T>) combiner, (T) identity);
    }

    @SuppressWarnings("unchecked")
    private void build(int node, int start, int end, T[] arr) {
        if (start == end) {
            tree[node] = arr[start];
        } else {
            int mid = (start + end) >>> 1;
            build(node << 1, start, mid, arr);
            build(node << 1 | 1, mid + 1, end, arr);
            T left = (T) tree[node << 1];
            T right = (T) tree[node << 1 | 1];
            tree[node] = combiner.apply(left, right);
        }
    }

    public void update(int idx, T value) {
        update(1, 0, n - 1, idx, value);
    }

    @SuppressWarnings("unchecked")
    private void update(int node, int start, int end, int idx, T value) {
        if (start == end) {
            tree[node] = value;
        } else {
            int mid = (start + end) >>> 1;
            if (idx <= mid) {
                update(node << 1, start, mid, idx, value);
            } else {
                update(node << 1 | 1, mid + 1, end, idx, value);
            }
            T left = (T) tree[node << 1];
            T right = (T) tree[node << 1 | 1];
            tree[node] = combiner.apply(left, right);
        }
    }

    @SuppressWarnings("unchecked")
    public T query(int l, int r) {
        return query(1, 0, n - 1, l, r);
    }

    @SuppressWarnings("unchecked")
    private T query(int node, int start, int end, int l, int r) {
        if (r < start || end < l) {
            return identity;
        }
        if (l <= start && end <= r) {
            return (T) tree[node];
        }
        int mid = (start + end) >>> 1;
        T left = query(node << 1, start, mid, l, r);
        T right = query(node << 1 | 1, mid + 1, end, l, r);
        return combiner.apply(left, right);
    }

    @SuppressWarnings("unchecked")
    public void printSegmentTree() {
        System.out.println("Segment Tree Detailed (node [start..end] = value):");
        printNode(1, 0, n - 1, "");
    }

    @SuppressWarnings("unchecked")
    private void printNode(int node, int start, int end, String indent) {
        T val = (T) tree[node];
        System.out.printf("%s[%d..%d] = %s%n", indent, start, end, val);
        if (start == end) return;
        int mid = (start + end) >>> 1;
        printNode(node << 1, start, mid, indent + "  ");
        printNode(node << 1 | 1, mid + 1, end, indent + "  ");
    }


    // Example usage:
    // Integer[] data = {2,1,5,3,4};
    // SegmentTree<Integer> st = new SegmentTree<>(data, Integer::sum, 0);
    // System.out.println(st.query(1,3));   // sum from index 1 to 3
}




