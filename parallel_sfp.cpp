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
#include <cstdint> 

using namespace std;

vector<unsigned long long> short_fastest_path(int V,int Vi, vector<vector<int>>& outgoing_edges, vector<vector<int>>& incoming_edges, vector<array<int,4>>& vertex_data, int source, int threads){
    // Initialize data structures
    vector<int> current;
    vector<pair<int,int>> LE(V,{-1,-1});  // List of pairs for each vertex
    vector<unsigned long long> journey(Vi, ((unsigned long long)INT_MAX << 32) | (uint32_t)INT_MAX );  // (journey_time, total_travel_time)
    vector<int> in_cnt(V, 0);
    journey[source] = {0}; // Journey to self is (0,0)
    // Thread-local next vectors for each thread
    vector<vector<int>> thread_next(threads);
    
    // Add all root vertices to current (vertices with no incoming edges)
    for (int v = 0; v < V; v++) {
        if (incoming_edges[v].empty()) {
            if(vertex_data[v][0] == source) {
                LE[v] = pair<int,int>(vertex_data[v][2], vertex_data[v][3]);
                journey[vertex_data[v][1]] = min(journey[vertex_data[v][1]], ((unsigned long long)vertex_data[v][3] << 32) | (uint32_t)vertex_data[v][3]);
            }
            for(auto& it : outgoing_edges[v]) {
                in_cnt[it]++;
                if(in_cnt[it]==incoming_edges[it].size()) {
                    current.push_back(it);
                }
            }
        }
    }
    
    while (!current.empty()) {
        int size = current.size();
        // Clear thread-local next vectors and reset in_next flags
        for (int t = 0; t < threads; t++) {
            thread_next[t].clear();
        }

        #pragma omp parallel num_threads(threads)
        {
            int thread_id = omp_get_thread_num();
            int chunk_size = (size + threads - 1) / threads; // Ceiling division
            int start_idx = thread_id * chunk_size;
            int end_idx = min(start_idx + chunk_size, size);
            
            for (int i = start_idx; i < end_idx; i++) {
                int v = current[i];

                // Process outgoing edges
                for (auto& it : outgoing_edges[v]) {
                    int prev_val;
                    #pragma omp atomic capture
                    prev_val = in_cnt[it]++;

                    if (prev_val + 1 == (int)incoming_edges[it].size()) {
                        thread_next[thread_id].push_back(it);
                    }
                }
            
                // Source initialization
                if (vertex_data[v][0] == source) {
                    LE[v] = {vertex_data[v][2], vertex_data[v][3]};
                    unsigned long long val = ((unsigned long long)vertex_data[v][3] << 32) | (uint32_t)vertex_data[v][3];

                    // #pragma omp critical
                    // {
                    //    journey[vertex_data[v][1]] = min(journey[vertex_data[v][1]], val);
                    // }
                
                    int idx = vertex_data[v][1];
                    unsigned long long old_val;
                                   
                    do {
                       #pragma omp atomic read
                       old_val = journey[idx];
                       if (val >= old_val) break;
                       #pragma omp atomic write
                       journey[idx] = val;
                    } while (val < old_val);

                }
            
                for (auto& u : incoming_edges[v]) {
                    if (LE[u] == pair<int, int>(-1, -1)) continue;
                    if (LE[v] == pair<int, int>(-1, -1) || LE[u].first > LE[v].first ||
                        (LE[u].first == LE[v].first && LE[v].second > LE[u].second + vertex_data[v][3])) {
                        
                        LE[v] = {LE[u].first, LE[u].second + vertex_data[v][3]};
                        int new_journey_time = vertex_data[v][2] + vertex_data[v][3] - LE[v].first;
                        unsigned long long val = ((unsigned long long)new_journey_time << 32) | (uint32_t)LE[v].second;
                        
                        // #pragma omp critical
                        // {
                        //     journey[vertex_data[v][1]] = min(journey[vertex_data[v][1]], val);
                        // }

                        int idx = vertex_data[v][1];
                        unsigned long long old_val;
                                       
                        do {
                           #pragma omp atomic read
                           old_val = journey[idx];
                           if (val >= old_val) break;
                           #pragma omp atomic write
                           journey[idx] = val;
                        } while (val < old_val);
                    }
                }
            }
        }

        current.clear();
        for (int t = 0; t < threads; t++) {
            for (int vertex : thread_next[t]) {
                    current.push_back(vertex);
            }
        }
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
        vector<unsigned long long> result = short_fastest_path(E,V, outgoing_edges, incoming_edges, edges, queryVertices[i], no_threads);
        end_time = omp_get_wtime();

        for (int j = 0; j < result.size(); j++) {
            outFile << queryVertices[i] << "--" << j << ": (" << ((result[j] >> 32) & 0xFFFFFFFF) << ", " << (result[j] & 0xFFFFFFFF) << ")\n";
        }
        queryTime = (end_time - start_time);
        totalTime += queryTime;
    }
    outFile.close();
    float averageTime = (totalTime / queryVertices.size()) * 1000;
    cout << "Average query time(" << output_filename << "):" << averageTime << " ms" << endl;
    return 0;
}