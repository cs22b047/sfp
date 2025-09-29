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

array<int, 3> find_predecessor(const array<int, 4>& edge, const vector<array<int, 3>>& Lvi) {
    int ti = edge[2]; 
    int r = Lvi.size();
    if (r == 0) {
        return {-1, -1, -1};
    }
    int j = 1;
    while (j <= r && Lvi[r - j][1] > ti) {
        j *= 2;
    }
    int low = max(r - j, 0);
    int high = r - j / 2 - 1;
    int ans = -1;

    while (low <= high) {
        int mid = (low + high) / 2;
        if (Lvi[mid][1] <= ti) {
            ans = mid;
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    if (ans != -1) {
        return Lvi[ans];
    }
    return {-1, -1, -1};
}

array<int,3> dominates(array<int, 3> a, array<int, 3> b) {
    if (a[0] > b[0] && a[1] <= b[1]) {
        return a; // a dominates b
    } else if (a[0] == b[0] && a[1] < b[1] && a[2] <= b[2]) {
        return a; // a dominates b
    } else if (a[0] == b[0] && a[1] == b[1] && a[2] <= b[2]) {
        return a; // a dominates b
    }
    if (a[0] < b[0] && a[1] >= b[1]) {
        return b; // b dominates a
    } else if (a[0] == b[0] && a[1] > b[1] && a[2] >= b[2]) {
        return b; // b dominates a
    } else if (a[0] == b[0] && a[1] == b[1] && a[2] > b[2]) {
        return b; // b dominates a
    }
    else {
        return {-1, -1, -1}; // Neither dominates the other
    }
}

vector<pair<int,int>> short_fastest_path(vector<vector<vector<array<int, 4>>>>& graph, int S, int V, int E, int no_threads) {
    vector<pair<int,int>> journey;
    vector<vector<array<int, 3>>> Lists;
    vector<mutex> list_mutexes(V); // Mutex for each vertex's list
    
    for (int i = 0; i < V; i++) {
        vector<array<int, 3>> item;
        Lists.push_back(item);
    }
    for (int i = 0; i < V; i++) {
        journey.push_back(make_pair(INT_MAX, INT_MAX));
    }
    journey[S] = make_pair(0, 0);

    // Get the maximum number of levels across all threads
    int max_levels = 0;
    for (int thread = 0; thread < no_threads; thread++) {
        max_levels = max(max_levels, (int)graph[thread].size());
    }

    // Process each level (for thread synchronization)
    for (int level = 0; level < max_levels; level++) {
        // Parallel processing within each level
        #pragma omp parallel num_threads(no_threads)
        {
            int thread_id = omp_get_thread_num();
            
            // Process edges for this thread at this level
            if (level < graph[thread_id].size()) {
                for (int edge_idx = 0; edge_idx < graph[thread_id][level].size(); edge_idx++) {
                    array<int, 4> current_edge = graph[thread_id][level][edge_idx];
                    
                    if (current_edge[1] == S) {
                        continue;
                    }
                    if (current_edge[0] == S) {
                        lock_guard<mutex> lock(list_mutexes[current_edge[0]]);
                        Lists[current_edge[0]].push_back({current_edge[2], current_edge[2], 0});
                    }
                    
                    array<int,3> y;
                    {
                        lock_guard<mutex> lock(list_mutexes[current_edge[0]]);
                        y = find_predecessor(current_edge, Lists[current_edge[0]]);

                    }
                    
                    if (y != array<int,3>{-1, -1, -1}) {
                        array<int, 3> T = {y[0], current_edge[2] + current_edge[3], y[2] + current_edge[3]};
  
                        // Lock the target vertex's list for modification
                        {
                            lock_guard<mutex> lock(list_mutexes[current_edge[1]]);
                            auto& targetList = Lists[current_edge[1]];
                            
                            auto it = targetList.begin();
                            while (it != targetList.end() && (it->at(0) < T[0] || (it->at(0) == T[0] && it->at(1) < T[1]))) {
                                ++it;
                            }

                            // Step 2: Check if T is dominated by any immediate neighbor (prev or next)
                            bool isDominated = false;

                            if (it != targetList.begin()) {
                                auto prev = it - 1;
                                if (dominates(*prev, T) == *prev) {
                                    isDominated = true;
                                }
                            }

                            if (it != targetList.end()) {
                                if (dominates(*it, T) == *it) {
                                    isDominated = true;
                                }
                            }

                            // Step 3: Insert T only if it is not dominated
                            if (!isDominated) {
                                it = targetList.insert(it, T);

                                // Step 4: Remove dominated previous element
                                while (it != targetList.begin()) {
                                    auto prev = it - 1;
                                    if (dominates(T, *prev) == T) {
                                        it = targetList.erase(prev);
                                    } else {
                                        break;
                                    }
                                }

                                // Step 5: Remove all dominated next elements
                                auto next = it + 1;
                                while (next != targetList.end() && dominates(T, *next) == T) {
                                    next = targetList.erase(next);
                                }
                            }

                            if (T[1]-T[0] < journey[current_edge[1]].first || (T[1]-T[0] == journey[current_edge[1]].first && T[2] < journey[current_edge[1]].second)) {
                                journey[current_edge[1]] = make_pair(T[1]-T[0], T[2]);
                            }
                        }
                    }
                }
            }
        }
        // Implicit barrier: all threads synchronize here before moving to next level
    }

    return journey;
}    

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <graph_file> <query_file> <output_file>" << endl;
        return 1;
    }

    string graph_filename = argv[1];
    string query_filename = argv[2];
    int no_threads = 1;
    
    // Set OpenMP thread count
    omp_set_num_threads(no_threads);
    
    ifstream infile(graph_filename);
    if (!infile) {
        cerr << "Failed to open file: " << graph_filename << endl;
        return 1;
    }

    int V, E, L, X;
    infile >> V >> E >> L >> X; // Read V, E, L, X

    vector<vector<vector<array<int, 4>>>> edges;
    for (int i = 0; i < no_threads; ++i) {
        vector<vector<array<int, 4>>> item;
        edges.push_back(item);
    }  
    int level = -1;
    int thread_id = 0;
    map<int,int> alloc;
    for (int i = 0; i < E; ++i) {
        int id, u, v, d, t, lvl;

        infile >> id >> u >> v >> d >> t >> lvl; // Read id, u, v, d, t, lvl
        if (lvl > level) {
            for (int j = 0; j < no_threads; ++j) {
                edges[j].push_back(vector<array<int, 4>>());
            } 
            alloc.clear(); // Clear allocation map for new level
            level = lvl;
            thread_id = 0; // Reset thread_id for the new level
        }
        if (lvl == level) {
            if (alloc.find(v) == alloc.end()) {
                alloc[v] = thread_id;
                edges[thread_id][lvl].push_back({u, v, d, t});
                thread_id = (thread_id + 1) % no_threads;
            } else {
                edges[alloc[v]][lvl].push_back({u, v, d, t});
            }
             // Cycle through threads
        }
    }

    infile.close();

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
    
    
    float totalTime, queryTime;
    totalTime = 0;
    double start, end;
    string output_filename = argv[3];
    ofstream outFile(output_filename, ios::out | ios::trunc);
    if (!outFile) {
        cerr << "Error opening file for writing.\n";
        return 1;
    }
    
    for (int i = 0; i < queryVertices.size(); i++) {
        start = omp_get_wtime();
        vector<pair<int,int>> result = short_fastest_path(edges, queryVertices[i], V, E, no_threads);
        end = omp_get_wtime();

        for (int j = 0; j < result.size(); j++) {
            outFile << queryVertices[i] << "--" << j << ": (" << result[j].first << ", " << result[j].second << ")\n";
        }
        queryTime = (end - start);
        totalTime += queryTime;
    }
    outFile.close();
    float averageTime = (totalTime / queryVertices.size()) * 1000;
    cout << "Average query time(" << output_filename << "):" << averageTime << " ms" << endl;
    return 0;
}