#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <cmath>
#include <limits>

using namespace std;

class Graph {
public:
    int V;  // Number of vertices
    vector<vector<int>> adj;  // Adjacency list
    const double dampingFactor = 0.85;
    const double epsilon = 1e-6;  // Convergence threshold for iterative methods

    Graph(int V) : V(V), adj(V) {}

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);  // For undirected graph
    }

    // Degree Centrality
    void degreeCentrality() {
        for (int i = 0; i < V; ++i) {
            cout << "Degree Centrality of node " << i << " : " << adj[i].size() << endl;
        }
    }

    // Closeness Centrality
    double closenessCentrality(int node) {
        vector<int> dist(V, -1);
        queue<int> q;
        q.push(node);
        dist[node] = 0;
        int totalDist = 0;
        int reachableNodes = 0;

        while (!q.empty()) {
            int curr = q.front();
            q.pop();

            for (int neighbor : adj[curr]) {
                if (dist[neighbor] == -1) {
                    q.push(neighbor);
                    dist[neighbor] = dist[curr] + 1;
                    totalDist += dist[neighbor];
                    reachableNodes++;
                }
            }
        }

        if (reachableNodes == 0) return 0.0;
        return reachableNodes / static_cast<double>(totalDist);
    }

    // Betweenness Centrality using Brandes' Algorithm
    vector<double> betweennessCentrality() {
        vector<double> betweenness(V, 0.0);

        for (int s = 0; s < V; ++s) {
            stack<int> S;
            vector<vector<int>> pred(V);  // Predecessors of each node
            vector<int> sigma(V, 0);      // Number of shortest paths
            vector<int> dist(V, -1);      // Distance from source
            vector<double> delta(V, 0.0); // Dependency score

            sigma[s] = 1;
            dist[s] = 0;
            queue<int> Q;
            Q.push(s);

            while (!Q.empty()) {
                int v = Q.front();
                Q.pop();
                S.push(v);
                for (int w : adj[v]) {
                    if (dist[w] < 0) {
                        Q.push(w);
                        dist[w] = dist[v] + 1;
                    }
                    if (dist[w] == dist[v] + 1) {
                        sigma[w] += sigma[v];
                        pred[w].push_back(v);
                    }
                }
            }

            while (!S.empty()) {
                int w = S.top();
                S.pop();
                for (int v : pred[w]) {
                    delta[v] += (sigma[v] / static_cast<double>(sigma[w])) * (1.0 + delta[w]);
                }
                if (w != s) {
                    betweenness[w] += delta[w];
                }
            }
        }

        for (double &val : betweenness) {
            val /= 2.0;  // Since the graph is undirected
        }

        return betweenness;
    }

    // Eigenvector Centrality
    vector<double> eigenvectorCentrality(int maxIterations = 100, double tolerance = 1e-6) {
        vector<double> centrality(V, 1.0);  // Initialize all centralities to 1
        vector<double> centralityOld(V, 0.0);

        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            centralityOld = centrality;

            // Update centrality based on neighbors
            for (int i = 0; i < V; ++i) {
                centrality[i] = 0.0;
                for (int neighbor : adj[i]) {
                    centrality[i] += centralityOld[neighbor];
                }
            }

            // Normalize
            double norm = 0.0;
            for (double c : centrality) norm += c * c;
            norm = sqrt(norm);
            for (double &c : centrality) c /= norm;

            // Check for convergence
            double diff = 0.0;
            for (int i = 0; i < V; ++i) {
                diff += abs(centrality[i] - centralityOld[i]);
            }
            if (diff < tolerance) break;
        }
        return centrality;
    }

    // PageRank
    vector<double> pageRank(int maxIterations = 100) {
        vector<double> rank(V, 1.0 / V);  // Initialize rank to 1/N
        vector<double> rankOld(V, 0.0);

        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            rankOld = rank;

            for (int i = 0; i < V; ++i) {
                rank[i] = (1.0 - dampingFactor) / V;
                for (int neighbor : adj[i]) {
                    rank[i] += dampingFactor * rankOld[neighbor] / adj[neighbor].size();
                }
            }

            // Check for convergence
            double diff = 0.0;
            for (int i = 0; i < V; ++i) {
                diff += abs(rank[i] - rankOld[i]);
            }
            if (diff < epsilon) break;
        }
        return rank;
    }

    // Katz Centrality
    vector<double> katzCentrality(double alpha = 0.01, double beta = 1.0, int maxIterations = 100) {
        vector<double> centrality(V, beta);  // Initialize all centralities to beta
        vector<double> centralityOld(V, 0.0);

        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            centralityOld = centrality;

            for (int i = 0; i < V; ++i) {
                centrality[i] = beta;
                for (int neighbor : adj[i]) {
                    centrality[i] += alpha * centralityOld[neighbor];
                }
            }

            // Check for convergence
            double diff = 0.0;
            for (int i = 0; i < V; ++i) {
                diff += abs(centrality[i] - centralityOld[i]);
            }
            if (diff < epsilon) break;
        }
        return centrality;
    }

    void computeCentralities() {
        cout << "\n--- Degree Centrality ---" << endl;
        degreeCentrality();

        cout << "\n--- Closeness Centrality ---" << endl;
        for (int i = 0; i < V; ++i) {
            cout << "Closeness Centrality of node " << i << " : " << closenessCentrality(i) << endl;
        }

        cout << "\n--- Betweenness Centrality ---" << endl;
        vector<double> betweenness = betweennessCentrality();
        for (int i = 0; i < V; ++i) {
            cout << "Betweenness Centrality of node " << i << " : " << betweenness[i] << endl;
        }

        cout << "\n--- Eigenvector Centrality ---" << endl;
        vector<double> eigen = eigenvectorCentrality();
        for (int i = 0; i < V; ++i) {
            cout << "Eigenvector Centrality of node " << i << " : " << eigen[i] << endl;
        }

        cout << "\n--- PageRank ---" << endl;
        vector<double> pagerank = pageRank();
        for (int i = 0; i < V; ++i) {
            cout << "PageRank of node " << i << " : " << pagerank[i] << endl;
        }

        cout << "\n--- Katz Centrality ---" << endl;
        vector<double> katz = katzCentrality();
        for (int i = 0; i < V; ++i) {
            cout << "Katz Centrality of node " << i << " : " << katz[i] << endl;
        }
    }
};

int main() {
    Graph g(6);

    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(3, 4);
    g.addEdge(4, 5);

    g.computeCentralities();

    return 0;
}
