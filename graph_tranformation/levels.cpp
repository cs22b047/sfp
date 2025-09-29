#include <iostream>
#include <vector>
#include <fstream>
#include <array>
#include <algorithm>
using namespace std;

vector<vector<int>> label_edges_from_esdg(int& V, int& E, int& L, vector<vector<int>>& outgoing_edges, vector<vector<int>>& incoming_edges, vector<int>& roots) {
    vector<vector<int>> levels_by_group; // Each index represents a level, contains vertices at that level
    
    // Create three sets as specified in algorithm
    vector<int> current;
    vector<int> next;
    vector<bool> visited(V, false);
    vector<bool> in_next(V, false); // Track which vertices are already in next to avoid duplicates
    
    // Add all root vertices to current set (step 2-4)
    for (int root : roots) {
        current.push_back(root);
    }
    
    int lvl = 1; // Start from level 1 (step 5)
    
    // Main loop while current is not empty (step 6)
    while (!current.empty()) {
        // Create a new level group for current level
        vector<int> current_level_vertices;
        
        // Process each vertex in current set (step 7) 
        for (int v : current) {
            visited[v] = true; // Mark as visited (step 8)
        }
        for (int v : current) {
            current_level_vertices.push_back(v); // Add to current level group (step 9)
            
            // Check each outgoing vertex (step 10)
            for (int v_prime : outgoing_edges[v]) {
                // Skip if already visited or already in next
                if (visited[v_prime] || in_next[v_prime]) {
                    continue;
                }
                
                // Check if all incoming vertices of v_prime are visited (step 11)
                bool all_incoming_visited = true;
                for (int incoming_v : incoming_edges[v_prime]) {
                    if (!visited[incoming_v]) {
                        all_incoming_visited = false;
                        break;
                    }
                }
                
                // If all incoming vertices are visited, add to next (step 12)
                if (all_incoming_visited) {
                    next.push_back(v_prime);
                    in_next[v_prime] = true; // Mark as added to next
                }
            }
        }
        
        sort(current_level_vertices.begin(), current_level_vertices.end());
        // Add current level vertices to the result
        levels_by_group.push_back(current_level_vertices);
        
        // Move to next level (step 16-18)
        current = next;
        next.clear();
        
        // Reset in_next flags for vertices that are now in current
        for (int v : current) {
            in_next[v] = false;
        }
        
        lvl++;
    }
    
    return levels_by_group; // Return vertices grouped by level (step 20)
}

int main(int argc, char* argv[]) {
    // Check if correct number of arguments are provided
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <esdg_file> <edges_file> <output_file>" << endl;
        return 1;
    }
    
    string esdg_filename = argv[1];
    string edges_filename = argv[2];
    string output_filename = argv[3];
    
    int V, E, Vi, R;
    vector<int> start, adj;
    vector<vector<int>> outgoing_edges, incoming_edges;
    
    // Read ESDG file
    ifstream infile(esdg_filename);
    if (!infile) {
        cerr << "Failed to open file: " << esdg_filename << endl;
        return 1;
    }
    
    infile >> Vi >> V >> E; // Read V, E, Vi
    start.resize(V + 1);
    adj.resize(E);
    
    for (int i = 0; i < V + 1; ++i) {
        infile >> start[i]; // Read start array
    }
    for (int i = 0; i < E; ++i) {
        infile >> adj[i]; // Read adj array
    }
    infile.close();

    // Build outgoing and incoming edges
    for(int i=0; i<V; i++) {
        outgoing_edges.push_back(vector<int>());
        incoming_edges.push_back(vector<int>());
    }
    for (int i = 0; i < V; ++i) {
        for (int j = start[i]; j < start[i + 1]; j++) {
            outgoing_edges[i].push_back(adj[j]);
            incoming_edges[adj[j]].push_back(i);
        }
    }
    
    vector<int> roots;
    for (int i = 0; i < V; ++i) {
        if (incoming_edges[i].empty()) { // If no incoming edges, it's a root
            roots.push_back(i);
        }
    }

    // Get levels from the algorithm
    vector<vector<int>> levels = label_edges_from_esdg(V, E, Vi, outgoing_edges, incoming_edges, roots);
    
    // Read edges file
    ifstream edges_file(edges_filename);
    if (!edges_file) {
        cerr << "Failed to open file: " << edges_filename << endl;
        return 1;
    }
    
    vector<array<int, 4>> edges;
    int temp_Vi, temp_V, temp_E;
    edges_file >> temp_Vi >> temp_V >> temp_E; // Read the first line from edges file
    for (int i = 0; i < E; ++i) {
        int id, u, v, d, t, offset;
        edges_file >> id >> u >> v >> d >> t >> offset; // Read each edge
        edges.push_back({u, v, d, t});
    }
    edges_file.close();
    
    // Write output file
    ofstream outfile(output_filename);
    if (!outfile) {
        cerr << "Failed to open file for writing: " << output_filename << endl;
        return 1;
    }

    outfile << Vi << " " << V << " " << E << " " << levels.size() << endl; // Write Vi, V, E, and number of levels
    int j = 0;
    for (int i = 0; i < levels.size(); i++) {
        for (int v : levels[i]) {
            outfile << j << " " << edges[v][0] << " " << edges[v][1] << " " << edges[v][2] << " " << edges[v][3] << " " << i << endl;
            j++;
        }
    }
    outfile.close();
    
    return 0;
}