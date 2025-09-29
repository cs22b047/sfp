#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <array>
#include <map>
#include <climits>
using namespace std;

array<vector<int>, 2> make_esdg(int V , int E , vector<vector<array<int,5>>> edges){
    vector<int> start(E + 1, 0);
    vector<int> adj;
    int edge_count = 0;
    
    // For each edge in the temporal graph
    for(int u = 0; u < V; u++){
        for(auto& edge : edges[u]){
            int edge_id = edge[0];    // Current edge ID
            int vu = edge[1];         // Source vertex
            int vv = edge[2];         // Destination vertex  
            int alpha = edge[3];      // Starting time
            int omega = edge[4] + edge[3];      // Arrival time
            
            // Find dependencies for this edge
            // Look at all outgoing edges from destination vertex vv
            if(vv < V && !edges[vv].empty()){
                int low = 0;
                int high = static_cast<int>(edges[vv].size()) - 1;  // Fix 1: Cast to int
                int earliest = -1; // Use -1 to indicate infinity/not found
                
                // Binary search for earliest valid connection
                while(low <= high){
                    int mid = (low + high) / 2;
                    auto& next_edge = edges[vv][mid];
                    int alpha_prime = next_edge[3]; // Starting time of next edge
                    
                    // Check if next edge starts after current edge ends
                    if(alpha_prime >= omega){ 
                        earliest = mid;
                        high = mid - 1; // Look for earlier valid edge
                    } else {
                        low = mid + 1;  // Look for later valid edge
                    }
                }
                
                // If we found a valid earliest connection
                if(earliest != -1){
                    map<int,vector<array<int, 5>>> adj_edges;
                    for(int i = earliest; i < static_cast<int>(edges[vv].size()); i++){  // Fix 2: Cast to int
                        auto& next_edge = edges[vv][i];
                        if(next_edge[3] >= omega){ // Only consider edges starting after current edge ends
                            adj_edges[next_edge[2]].push_back(next_edge); // Add edge ID to adjacency list
                        }
                    }
                    for(auto& entry : adj_edges){
                        vector<array<int,5>> v_next = entry.second;
                        sort(v_next.begin(), v_next.end(), [](const array<int,5>& a, const array<int,5>& b) {
                            return a[3]+a[4] < b[3]+b[4]; // Sort by arrival time
                        });
                        int t = INT_MAX;
                        for(auto& next_edge : v_next){
                            if(next_edge[4] < t){
                                adj.push_back(next_edge[0]);
                                t = next_edge[4];
                                edge_count++;
                            }
                        }
                    }
                }
            }
            
            start[edge_id + 1] = edge_count;
        }
    }
    
    return {start, adj};
}

int main(int argc ,char* argv[]) {
    int V, E;
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <graph_file> <esdg_file> <edges_file>" << endl;
        return 1;
    }
    vector<array<int, 4>> graph; // Graph structure to hold edges
    vector<vector<array<int, 5>>> edges;
    
    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Failed to open file: " << argv[1] << endl;
        return 1;
    }

    infile >> V >> E; // Read V, E
    graph.resize(E);
    for (int i = 0; i < E; ++i) {
        int x , u, v, d, t;
        infile >> x >> u >> v >> d >> t; // Read each edge
        graph[i] = {u, v, d, t}; // Store edge in the graph
    }
    infile.close();

    sort(graph.begin(), graph.end(), [](const array<int, 4>& a, const array<int, 4>& b) {
        return a[0] < b[0] || (a[0] == b[0] && a[2] < b[2]); // Sort by start vertex and starting time
    });

    // Initialize edges structure with proper sorting
    edges.resize(V);
    int id = 0;
    for (const auto& edge : graph) {
        int u = edge[0], v = edge[1], d = edge[2], t = edge[3];
        edges[u].push_back({id, u, v, d, t});
        id++;
    }

    auto esdg = make_esdg(V, static_cast<int>(graph.size()), edges); // Fix 5: Cast to int
    vector<int> start = esdg[0];
    vector<int> adj = esdg[1];

    ofstream outfile(argv[2]);
    if (!outfile) {
        cerr << "Failed to open file for writing: " << argv[2] << endl;
        return 1;
    }

    outfile << V << " " << graph.size() << " " << adj.size()<< endl; 
    for (int i = 0; i < static_cast<int>(graph.size()) + 1; ++i) {  // Fix 6: Cast to int
        outfile << start[i] << " "; // Write start array
    }
    outfile << endl;
    for(int i = 0; i < static_cast<int>(adj.size()); ++i) {  // Fix 7: Cast to int
        outfile << adj[i] << " "; // Write adjacency list
    }
    outfile << endl;
    outfile.close();

    ofstream outfile2(argv[3]);
    if (!outfile2) {    
        cerr << "Failed to open file for writing: " << argv[3] << endl;
        return 1;
    }
    outfile2 << V << " " << graph.size() << " " << adj.size() << endl; // Write V and actual edge count
    int edge_id = 0;
    int i=1;
    for (const auto& edge : graph) {
        outfile2 << edge_id << " " << edge[0] << " " << edge[1] << " " << edge[2] << " " << edge[3] << " " << start[i] << endl; // Write each edge
        edge_id++;  // Fix 8: Increment edge_id
        i++;
    }
    for(int i = 0; i < static_cast<int>(adj.size()); ++i) {  // Fix 7: Cast to int
        outfile2 << adj[i] << " "; // Write adjacency list
    }
    outfile2 << endl;
    outfile2.close();
    
    return 0;  // Fix 9: Add return statement
}