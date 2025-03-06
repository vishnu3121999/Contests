import java.io.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.*;
import java.util.logging.LogManager;
import java.util.stream.Collectors;
import static java.lang.Math.*;


public class Main implements Runnable {
    int mod = 1000000007;
    int mod2 = 998244353;
    FastReader fr = new FastReader(System.in);
    public static void main(String[] args) throws FileNotFoundException {
        new Thread(null, new Main(), "Main", 1 << 28).start();
    }

    @Override
    public void run() {
//        print(Generators.generateArray(new ArrayProps(5,1,10,true)));
        print(pow(2,25));

//        int T = fr.nextInt();
//        while(T-->0) {
//            int n = fr.nextInt();
//            StringBuilder sb = new StringBuilder();
//            HashMap<Integer, Integer> map = new HashMap<>();
//            sb.append("! ");
//            for (int i = 2; i <= n; i++) {
//
//                int a = 1;
//                int b = i;
//
//                int x = ask(a,b);
//                while(x!=a){
//                    a=x;
//                    x=ask(a,b);
//                }
//                if(map.containsKey(b))continue;
//                map.put(b,x);
//                sb.append(x).append(" ").append(b).append(" ");
////                    print(x,b);
//            }
//            print(sb);
//            System.out.flush();
//        }

    }

    private int ask(int a, int b) {
        System.out.println("? "+a+" "+b);
        System.out.flush();
        return fr.nextInt();
    }


    public int upper(ArrayList<Integer> nums, int target) {
        int l=0,h=nums.size();
        while(l<h){
            int mid = l+(h-l)/2;
            if(nums.get(mid)>target)h=mid;
            else l=mid+1;
        }
        return l-1;
    }

    public int lower(ArrayList<Integer> nums, int target) {
        int l=0,h=nums.size();
        while(l<h){
            int mid = l+(h-l)/2;
            if(nums.get(mid)>=target)h=mid;
            else l=mid+1;
        }
        return l;
    }

    static int binarSearch(ArrayList<Integer> l,int target){
        int n = l.size();
        int i=0,j=n;
        while(i<j){
            int mid = i+(j-i)/2;
            if(l.get(mid)>target)j=mid;
            else i=mid+1;
        }
        return i-1;
    }

    static int bs1(int[] arr ,long target){
        int n = arr.length;
        int l =0,h=n;
        while(l<h){
            int mid = l+(h-l)/2;
            if(arr[mid]>=target)h=mid;
            else l=mid+1;
        }
        return l-1;
    }

    static int bs2(int[] arr ,long target){
        int n = arr.length;
        int l =0,h=n;
        while(l<h){
            int mid = l+(h-l)/2;
            if(arr[mid]>=target)h=mid;
            else l=mid+1;
        }
        return l;
    }


    public static HashMap<Integer, Integer> sortByValue(HashMap<Integer, Integer> hm)
    {
        HashMap<Integer, Integer> temp
                = hm.entrySet()
                .stream()
                .sorted((i1, i2)
                        -> i1.getValue().compareTo(
                        i2.getValue()))
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        Map.Entry::getValue,
                        (e1, e2) -> e1,
                        LinkedHashMap::new));
        return temp;
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

    static void print(Object... args){
        for(var o:args){
            if(o.getClass().isArray()){
                if (o instanceof int[])     System.out.println(Arrays.toString((int[]) o));
                if (o instanceof long[])    System.out.println(Arrays.toString((long[]) o));
                if (o instanceof double[])  System.out.println(Arrays.toString((double[]) o));
                if (o instanceof float[])   System.out.println(Arrays.toString((float[]) o));
                if (o instanceof boolean[]) System.out.println(Arrays.toString((boolean[]) o));
                if (o instanceof char[])    System.out.println(Arrays.toString((char[]) o));
            }
            else System.out.print(o+" ");
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

    // Method to compute the sum of elements in a long array
    public static long sum(long[] arr) {
        long sum = 0;
        for (long val : arr) {
            sum += val;
        }
        return sum;
    }

    // Method to compute the sum of elements in an ArrayList<Long>
    public static long sum(ArrayList<Long> list) {
        long sum = 0;
        for (long val : list) {
            sum += val;
        }
        return sum;
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

    // Method to compute the sum of elements in an int array
    public static int sum(int[] arr) {
        int sum = 0;
        for (int val : arr) {
            sum += val;
        }
        return sum;
    }

    public static long[] prefixSum(int[] ai) {
        long sum = 0;
        long[] ps = new long[ai.length];
        for (int i = 0; i < ai.length; i++) {
            sum += ai[i];
            ps[i] = sum;
        }
        return ps;
    }

    public static long[] suffixSum(long[] ar) {
        long sum = 0;
        long[] ss = new long[ar.length];
        for (int i = ar.length - 1; i >= 0; i--) {
            sum += ar[i];
            ss[i] = sum;
        }
        return ss;
    }
    public static long[] suffixSum(int[] ar) {
        long sum = 0;
        long[] ss = new long[ar.length];
        for (int i = ar.length - 1; i >= 0; i--) {
            sum += ar[i];
            ss[i] = sum;
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

    public static boolean isPerfectSquare(long num) {
        if (num < 0) return false;
        long sqrt = (long) Math.sqrt(num);
        return sqrt * sqrt == num;
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

class DSU_ {
    public int[] parent;
    public int[] size;

    // Time Complexity: O(n)
    // Space Complexity: O(n)
    public DSU_(int n) {
        parent = new int[n];  // Zero-based indexing
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

}