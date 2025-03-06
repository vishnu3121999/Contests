import java.util.HashMap;
import java.util.Map;

public class DSU_HashMap_<T> {
    private Map<T, T> parent;
    private Map<T, Integer> rank;

    // Constructor to create and initialize a disjoint set
    // Time Complexity: O(n)
    // Space Complexity: O(n)
    public DSU_HashMap_() {
        parent = new HashMap<>();
        rank = new HashMap<>();
    }

    // Create a set with a single element
    // Time Complexity: O(1)
    // Space Complexity: O(1)
    public void makeSet(T x) {
        parent.put(x, x);
        rank.put(x, 0);
    }

    // Find the root of the set containing x with path compression
    // Time Complexity: O(log n) (amortized)
    // Space Complexity: O(log n) (due to recursive stack)
    public T find(T x) {
        if (!parent.get(x).equals(x)) {
            parent.put(x, find(parent.get(x)));  // Path compression
        }
        return parent.get(x);
    }

    // Union two sets by rank
    // Time Complexity: O(1) (amortized)
    // Space Complexity: O(1)
    public void union(T x, T y) {
        T rootX = find(x);
        T rootY = find(y);

        if (!rootX.equals(rootY)) {
            // Union by rank
            int rankX = rank.get(rootX);
            int rankY = rank.get(rootY);

            if (rankX > rankY) {
                parent.put(rootY, rootX);
            } else if (rankX < rankY) {
                parent.put(rootX, rootY);
            } else {
                parent.put(rootY, rootX);
                rank.put(rootX, rankX + 1);
            }
        }
    }

    // Check if two elements are in the same set
    // Time Complexity: O(log n) (amortized)
    // Space Complexity: O(1)
    public boolean isConnected(T x, T y) {
        return find(x).equals(find(y));
    }

    // Main method for testing
    public static void main(String[] args) {
        DSU_HashMap_<String> ds = new DSU_HashMap_<>();

        // Create sets
        ds.makeSet("A");
        ds.makeSet("B");
        ds.makeSet("C");
        ds.makeSet("D");
        ds.makeSet("E");

        ds.union("A", "B");
        ds.union("C", "D");
        ds.union("B", "C");

        System.out.println("A and D are connected: " + ds.isConnected("A", "D")); // true
        System.out.println("A and E are connected: " + ds.isConnected("A", "E")); // false
    }
}
