#include <vector>
#include <array>
#include <utility>
#include <fstream>
#include <limits.h>
#include <iostream>
#include <algorithm>
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
    } else if (a[0] == b[0] && a[1] == b[1] && a[2] >= b[2]) {
        return b; // b dominates a
    }
    else {
        return {-1, -1, -1}; // Neither dominates the other
    }
    
}

vector<pair<int,int>> short_fastest_path(vector<array<int, 4>>& graph,int S, int V, int E) {
    vector<pair<int,int>> journey;
    vector<vector<array<int, 3>>> Lists;
    for (int i = 0; i < V; i++) {
        vector<array<int, 3>> item;
        Lists.push_back(item);
    }
    for (int i = 0; i < V; i++) {
        journey.push_back(make_pair(INT_MAX,INT_MAX));
    }
    journey[S] = make_pair(0, 0);

    for(int i = 0; i < E; i++) {
        if (graph[i][1] == S) {
            continue;
        }
        if (graph[i][0] == S) {
            Lists[graph[i][0]].push_back({graph[i][2], graph[i][2], 0});
        }
        
        array<int,3> y =find_predecessor(graph[i],Lists[graph[i][0]]);
        if (y!= array<int,3>{-1, -1, -1}) {
            array<int, 3> T = {y[0], graph[i][2] + graph[i][3], y[2] + graph[i][3]};
            auto& targetList = Lists[graph[i][1]];
            
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


            if (T[1]-T[0] < journey[graph[i][1]].first || (T[1]-T[0] == journey[graph[i][1]].first && T[2] < journey[graph[i][1]].second)) {
                journey[graph[i][1]] = make_pair(T[1]-T[0], T[2]);
            }
        }
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

    ifstream infile(graph_filename);
    if (!infile) {
        cerr << "Failed to open file: " << graph_filename << endl;
        return 1;
    }

    int V, E;
    infile >> V >> E;

    vector<array<int, 4>> edges;

    for (int i = 0; i < E; ++i) {
        int u, v, d, t;
        infile >> u >> v >> d >> t;
        edges.push_back({u, v, d, t});
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
    
    sort(edges.begin(), edges.end(), [](const array<int, 4>& a, const array<int, 4>& b) {
        return (a[2]) < (b[2]);
    });
    float totalTime ,queryTime;
    totalTime = 0;
    clock_t start , end;
    string output_filename = argv[3]; // Change this to your output filename if needed
    ofstream outFile(output_filename, ios::out | ios::trunc);
    if (!outFile) {
        cerr << "Error opening file for writing.\n";
        return 1;
    }
    for (int i = 0; i < queryVertices.size(); i++) {
        start = clock();
        vector<pair<int,int>> result = short_fastest_path(edges, queryVertices[i], V, E);
        end =clock();

        for (int j = 0;  j< result.size(); j++) {
            outFile << queryVertices[i]<< "--" << j << ": (" << result[j].first << ", " << result[j].second << ")\n";
        }
        queryTime = double(end - start) / CLOCKS_PER_SEC;
        totalTime += queryTime;
    }
    outFile.close();
    float averageTime = (totalTime / queryVertices.size())*1000;
    cout << "Average query time("<<output_filename<<"):" << averageTime << " ms" << endl;
    return 0;
}
