import java.util.Arrays;

public class DSU_ {
    private int[] parent;
    private int[] rank;


    // Time Complexity: O(n)
    // Space Complexity: O(n)
    public DSU_(int n) {
        parent = new int[n];  // Zero-based indexing
        rank = new int[n];
        for (int i = 0; i < n; i++) {
            parent[i] = i;
            rank[i] = 0;
        }
    }

    // Find the root of the set containing x with path compression
    // Time Complexity: O(log n) (amortized)
    // Space Complexity: O(log n) (due to recursive stack)
    public int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);  // Path compression
        }
        return parent[x];
    }

    // Union two sets by rank
    // Time Complexity: O(1) (amortized)
    // Space Complexity: O(1)
    public void union(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX != rootY) {
            // Union by rank
            if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;
            } else if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;
            } else {
                parent[rootY] = rootX;
                rank[rootX]++;
            }
        }
    }

    // Check if two elements are in the same set
    // Time Complexity: O(log n) (amortized)
    // Space Complexity: O(1)
    public boolean isConnected(int x, int y) {
        return find(x) == find(y);
    }

    // Main method for testing
    public static void main(String[] args) {
        int n = 5;  // Number of elements
        DSU_ ds = new DSU_(n);

        ds.union(0, 4);
        ds.union(1, 3);
        ds.union(2, 0);
        ds.union(3, 2);
        ds.union(4, 1);
        for (int i = 0; i <n; i++) {
            ds.find(i);
        }
        System.out.println(Arrays.toString(ds.parent));
//        System.out.println("0 and 3 are connected: " + ds.isConnected(0, 3)); // true
//        System.out.println("0 and 4 are connected: " + ds.isConnected(0, 4)); // false
    }
}


