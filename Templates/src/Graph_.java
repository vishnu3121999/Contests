import java.util.*;


// Zero-based Indexing of Vertices
public class Graph_ {

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
