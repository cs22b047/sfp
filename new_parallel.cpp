#include <vector>
#include <array>
#include <utility>
#include <fstream>
#include <limits.h>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <omp.h>
#include <mutex>
#include <map>
using namespace std;

vector<pair<int, int>> short_fastest_path(int V,int Vi, vector<vector<int>>& outgoing_edges, vector<vector<int>>& incoming_edges, vector<array<int,4>>& vertex_data, int source, int threads) {
    // Initialize data structures
    vector<int> current, next;
    vector<bool> visited(V, false);  // Visited vertices
    vector<bool> in_next(V, false); // Track vertices in next
    vector<pair<int,int>> LE(V,{-1,-1});  // List of pairs for each vertex
    vector<pair<int, int>> journey(Vi, {INT_MAX, INT_MAX});  // (journey_time, total_travel_time)
    
    // Thread-local next vectors for each thread
    vector<vector<int>> thread_next(threads);
    
    // Mutexes for thread safety
    mutex journey_mutex;
    mutex le_mutex;
    mutex next_mutex;
    
    journey[source] = {0, 0}; // Journey to self is (0,0)
    // Add all root vertices to current (vertices with no incoming edges)
    for (int v = 0; v < V; v++) {
        if (incoming_edges[v].empty()) {
            current.push_back(v);
        }
    }
    
    while (!current.empty()) {
        // Mark all vertices in current as visited (parallelized)
        int s = current.size();
        for (int i = 0; i < s; i++) {
            visited[current[i]] = true;
        }
        
        // Clear thread-local next vectors and reset in_next flags
        for (int t = 0; t < threads; t++) {
            thread_next[t].clear();
        }
        fill(in_next.begin(), in_next.end(), false);
        
        // Process each vertex in current (parallelized)
        #pragma omp parallel num_threads(threads)
        {
            int thread_id = omp_get_thread_num();
            int chunk_size = (s + threads - 1) / threads; // Ceiling division
            int start_idx = thread_id * chunk_size;
            int end_idx = min(start_idx + chunk_size, s);
            
            for (int i = start_idx; i < end_idx; i++) {
                int v = current[i];
                if(vertex_data[v]==array<int,4>{140, 139, 66480, 240}) {
                    //cout<<LE[v].first<<", "<<LE[v].second<<endl;
                }
                //cout<<vertex_data[v][0]<<", "<<source<<endl;
                // If vertex starts from source

                if (vertex_data[v][0] == source) {
                    {
                        lock_guard<mutex> lock(le_mutex);
                        LE[v] = pair<int,int>(vertex_data[v][2], vertex_data[v][3]);
                    }
                    if (vertex_data[v]==array<int,4>{140, 139, 66480, 240}) {
                        //cout<<LE[v].first<<", "<<LE[v].second<<endl;
                    }
                    {
                        lock_guard<mutex> lock(journey_mutex);
                        if (journey[vertex_data[v][1]].first > vertex_data[v][3]) {
                            journey[vertex_data[v][1]] = {vertex_data[v][3], vertex_data[v][3]};
                        }
                    }
                }

                for (int u :outgoing_edges[v]) {
                    // Check if LE[v] is set
                    // Check if all incoming vertices of u are visited
                    bool all_incoming_visited = true;
                    for (int incoming : incoming_edges[u]) {
                        if (!visited[incoming]) {
                            all_incoming_visited = false;
                            break;
                        }
                    }
                    
                    if (all_incoming_visited) {
                        // Use atomic operation to check and set in_next flag
                        lock_guard<mutex> lock(next_mutex);
                        if (!in_next[u]) {
                            in_next[u] = true;
                            thread_next[thread_id].push_back(u);
                        }
                    }
                }
                if (LE[v] == pair<int,int>(-1,-1)) continue; // Skip if LE[u] is already set
                // Process outgoing edges
                for (int u : outgoing_edges[v]) {
                    // Update LE[u] if conditions are met
                    bool should_update = false;
                    pair<int,int> current_le_v, current_le_u;
                    
                    {
                        lock_guard<mutex> lock(le_mutex);
                        current_le_v = LE[v];
                        current_le_u = LE[u];
                        
                        if (current_le_u == pair<int,int>(-1,-1)) {
                            should_update = true;
                        } else if (current_le_u.first < current_le_v.first) {
                            should_update = true;
                        } else if (current_le_u.first == current_le_v.first && 
                                  current_le_u.second > current_le_v.second + vertex_data[u][3]) {
                            should_update = true;
                        }
                        if (should_update) {
                            LE[u] = pair<int,int>(current_le_v.first, current_le_v.second + vertex_data[u][3]);
                        }
                    }
                    
                    if (should_update) {
                        pair<int,int> new_le_u;
                        {
                            lock_guard<mutex> lock(le_mutex);
                            new_le_u = LE[u];
                        }
                        
                        int new_journey_time = vertex_data[u][3] + vertex_data[u][2] - new_le_u.first;
                        
                        {
                            lock_guard<mutex> lock(journey_mutex);
                            if (journey[vertex_data[u][1]].first > new_journey_time || (journey[vertex_data[u][1]].first == new_journey_time && journey[vertex_data[u][1]].second > new_le_u.second)) {
                                journey[vertex_data[u][1]] = {new_journey_time, new_le_u.second};
                            }
                        }
                    }
                }
            }
        }
        
        // Combine all thread-local next vectors into the main next vector
        next.clear();
        for (int t = 0; t < threads; t++) {
            for (int vertex : thread_next[t]) {
                next.push_back(vertex);
            }
        }
        
        // Move to next level
        current = next;
        next.clear();
    }
    
    return journey;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <graph_file> <query_file> <output_file> [num_threads]" << endl;
        return 1;
    }

    string graph_filename = argv[1];
    string query_filename = argv[2];
    string output_filename = argv[3];
    int no_threads = (argc > 4) ? atoi(argv[4]) : 1;
    
    // Set OpenMP thread count
    omp_set_num_threads(no_threads);

    vector<int> start, adj;
    vector<array<int, 4>> edges;

    start.push_back(0); // Initialize start vector with 0

    ifstream infile(graph_filename);
    if (!infile) {
        cerr << "Failed to open file: " << graph_filename << endl;
        return 1;
    }

    int V, E, Ex;
    infile >> V >> E >> Ex ; // Read V, E, L, X

    for (int i = 0; i < E; ++i) {
        int id, u, v, d, t, offset;

        infile >> id >> u >> v >> d >> t >> offset; // Read id, u, v, d, t, offset
        edges.push_back({u, v, d, t});
        start.push_back(offset); // Store offset in start vector
    }
    for (int i = 0; i < Ex; ++i) {
        int s;
        infile >> s; // Read start array
        adj.push_back(s);
    }

    infile.close();
    // cout<<V<<", "<<E<<", "<<Ex<<endl;
    // cout<<edges[0][0]<<" " <<edges[0][1]<<", "<<edges[0][2]<<", "<<edges[0][3]<<endl;

    vector<int> queryVertices;
    ifstream infile2(query_filename);

    if (!infile2) {
        cerr << "Failed to open file: " << query_filename << endl;
        return 1;
    }
    int numQueries;
    infile2 >> numQueries; // Read the number of queries
    for (int i = 0; i < numQueries; ++i) {
        int s;
        infile2 >> s;
        queryVertices.push_back(s);
    }

    infile2.close();
    vector<vector<int>> outgoing_edges, incoming_edges;
    for(int i=0; i<E; i++) {
        outgoing_edges.push_back(vector<int>());
        incoming_edges.push_back(vector<int>());
    }
    for (int i = 0; i < E; ++i) {
        for (int j = start[i]; j < start[i + 1]; j++) {
            outgoing_edges[i].push_back(adj[j]);
            incoming_edges[adj[j]].push_back(i);
        }
    }

    float totalTime, queryTime;
    totalTime = 0;
    double start_time, end_time;
    ofstream outFile(output_filename, ios::out | ios::trunc);
    if (!outFile) {
        cerr << "Error opening file for writing.\n";
        return 1;
    }
    
    for (int i = 0; i < queryVertices.size(); i++) {
        start_time = omp_get_wtime();
        vector<pair<int,int>> result = short_fastest_path(E,V, outgoing_edges, incoming_edges, edges, queryVertices[i], no_threads);
        end_time = omp_get_wtime();

        for (int j = 0; j < result.size(); j++) {
            outFile << queryVertices[i] << "--" << j << ": (" << result[j].first << ", " << result[j].second << ")\n";
        }
        queryTime = (end_time - start_time);
        totalTime += queryTime;
    }
    outFile.close();
    float averageTime = (totalTime / queryVertices.size()) * 1000;
    cout << "Average query time(" << output_filename << "):" << averageTime << " ms" << endl;
    return 0;
}