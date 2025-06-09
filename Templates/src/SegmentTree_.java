public class SegmentTree_ {
    // Array-based iterative segment tree for range sum queries and point updates
    private final int size;
    private final int[] tree;

    // Initialize with the length of the input array
    public SegmentTree_(int n) {
        int s = 1;
        while (s < n) s <<= 1;
        size = s;
        tree = new int[size * 2];
    }

    // Build the segment tree from an array
    public void build(int[] arr) {
        int n = arr.length;
        for (int i = 0; i < n; i++) {
            tree[size + i] = arr[i];
        }
        for (int i = size - 1; i > 0; i--) {
            tree[i] = tree[i * 2] + tree[i * 2 + 1];
        }
    }

    // Point update: set arr[index] = value
    public void update(int index, int value) {
        int pos = index + size;
        tree[pos] = value;
        for (pos >>= 1; pos > 0; pos >>= 1) {
            tree[pos] = tree[pos * 2] + tree[pos * 2 + 1];
        }
    }

    // Range sum query on [l, r] (inclusive)
    public int query(int l, int r) {
        int res = 0;
        l += size;
        r += size;
        while (l <= r) {
            if ((l & 1) == 1) res += tree[l++];
            if ((r & 1) == 0) res += tree[r--];
            l >>= 1;
            r >>= 1;
        }
        return res;
    }
}
