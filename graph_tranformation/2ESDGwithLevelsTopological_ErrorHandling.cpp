/*
g++ -o addLevel 2ESDGwithLevelsTopological.cpp
./addLevel ../../../Data/dependencies/sample1_dep.txt ../../../Data/sample1_esdg_level1.txt

*/

#include<iostream>
#include <list>
#include <queue>
#include<fstream>
#include<vector>
#include<algorithm>
#include<stack>
using namespace std;

typedef struct conn
{
	int u,v,departure,duration;
}conn;

class Graph
{
	int vCount; // number of vertices in Temporal Graph
    int nCount; // number of nodes in ESDG Graph or connections in Temporal Graph
    int totalDependencies; // number of edges in ESDG Graph;
    vector<conn> connectionsVector;
    vector<int> dependency;
    vector<int> offset;
    vector<int> *adj;
    vector<int> level;
    vector<bool> processed; // Track whether a vertex has been processed

public:
    Graph(int n, int m,int td)
    {
        vCount=n;
        nCount=m;
        totalDependencies=td;
        adj = new vector<int>[nCount];
    }
    void addEdge(int u,int v);
    void readESDGFromFile(ifstream &connections);
    void topologicalSort();
    void printLevelNumbers(ofstream &connectionsWithLevels);
    bool isCyclicUtil(int v, vector<bool>& visited, vector<bool>& recStack);
    bool isCyclic();
	bool detectCycle();
    bool detectCycleUtil(int v, vector<bool>& visited, vector<int>& parent, vector<bool>& recStack);
    void printCycle(int start, int end, const vector<int>& parent);
}; // class declaration

void Graph::addEdge(int u,int v)
{
	if (u < 0 || u >= nCount || v < 0 || v >= nCount) 
	{
        throw invalid_argument("Vertex index out of bounds");
    }
    adj[u].push_back(v);
}

void Graph::readESDGFromFile(ifstream &connections)
{
	int id;

	cout<<"In readESDGFromFile "<<vCount<<" "<<nCount<<" "<<totalDependencies<<endl;

    for(int i=0;i<nCount;i++)
	{
		conn currConnection;
		connections >> id >> currConnection.u >> currConnection.v >> currConnection.departure >> currConnection.duration;
		connectionsVector.push_back(currConnection);

		int currOffset;
		connections >> currOffset;

		offset.push_back(currOffset);
	}

	cout<<nCount<<endl;
	offset.push_back(totalDependencies);

	for(int i=0;i<totalDependencies;i++)
	{
		int currDependency;
		connections >> currDependency;
		dependency.push_back(currDependency);
	}

	for(int i=0;i<nCount;i++)
	{
		for(int j=offset[i];j<offset[i+1];j++)
			addEdge(i,dependency[j]);
	}
}


bool Graph::isCyclicUtil(int v, vector<bool>& visited, vector<bool>& recStack) {
    if (!visited[v]) {
        visited[v] = true;
        recStack[v] = true;

        for (int neighbor : adj[v]) {
            if (!visited[neighbor] && isCyclicUtil(neighbor, visited, recStack))
                return true;
            else if (recStack[neighbor])
                return true;
        }
    }
    recStack[v] = false;
    return false;
}

bool Graph::isCyclic() {
    vector<bool> visited(nCount, false);
    vector<bool> recStack(nCount, false);

    for (int i = 0; i < nCount; i++)
        if (isCyclicUtil(i, visited, recStack))
            return true;

    return false;
}


bool Graph::detectCycle() {
    vector<bool> visited(nCount, false);
    vector<int> parent(nCount, -1);   // Stores the parent of each node for cycle reconstruction
    vector<bool> recStack(nCount, false);

    for (int i = 0; i < nCount; i++) {
        if (!visited[i]) {
            if (detectCycleUtil(i, visited, parent, recStack)) {
                return true;
            }
        }
    }
    return false;
}

bool Graph::detectCycleUtil(int v, vector<bool>& visited, vector<int>& parent, vector<bool>& recStack) {
    visited[v] = true;
    recStack[v] = true;

    for (int neighbor : adj[v]) {
        if (!visited[neighbor]) {
            parent[neighbor] = v;
            if (detectCycleUtil(neighbor, visited, parent, recStack)) {
                return true;
            }
        } else if (recStack[neighbor]) {
            // Cycle detected, print the cycle
            printCycle(neighbor, v, parent);
            return true;
        }
    }

    recStack[v] = false;
    return false;
}

void Graph::printCycle(int start, int end, const vector<int>& parent) {
    stack<int> cycleStack;
    int current = end;

    cycleStack.push(start);
    while (current != start) {
        cycleStack.push(current);
        current = parent[current];
    }

    // Print the cycle
    cout << "Cycle detected: ";
    while (!cycleStack.empty()) {
        cout << cycleStack.top() << " ";
        cout<<"Connection Details: "<<connectionsVector[cycleStack.top()].u<<" "<<connectionsVector[cycleStack.top()].v<<" "<<connectionsVector[cycleStack.top()].departure<<" "<<connectionsVector[cycleStack.top()].duration<<endl;
        cycleStack.pop();
    }
    cout << start << endl;  // Complete the cycle by printing the start vertex again


}

void Graph::topologicalSort() 
{
    vector<int> indegree(nCount, 0); // Initialize indegree vector to 0

    level.assign(nCount, -1);     // Initialize level vector to -1 for all vertices
    processed.assign(nCount, false); // Initialize processed vector to false for all vertices
    // Calculate indegrees for each vertex
    for (int i = 0; i < nCount; i++) {
        for (auto it = adj[i].begin(); it != adj[i].end(); it++) {
            indegree[*it]++;
        }
    }

    cout << "Indegree Calculated Successfully" << endl;

    // // Print indegree values
    // for (int i = 0; i < nCount; i++) {
    //     cout << indegree[i] << " ";    
    // }
    // cout << endl;

   queue<int> q;

    for (int i = 0; i < nCount; i++) {
        if (indegree[i] == 0) {
            q.push(i);
            level[i] = 0; // Initialize level for source nodes
            processed[i] = true; // Mark as processed
        }
    }
    
	//////////////////////////////////////////////////////////////////////////

    /*vector<bool> reachable(nCount, false);
	queue<int> q;
	for (int i = 0; i < nCount; i++) {
	    if (indegree[i] == 0) {
	        q.push(i);
	        reachable[i] = true;
	    }
	}

	while (!q.empty()) {
	    int u = q.front(); q.pop();
	    for (int v : adj[u]) {
	        if (!reachable[v]) {
	            reachable[v] = true;
	            q.push(v);
	        }
	    }
	}

	// Check for unreachable vertices
	for (int i = 0; i < nCount; i++) {
	    if (!reachable[i]) {
	        cout << "Vertex " << i << " is not reachable from any source node." << endl;
	    }
	}




    //////////////////////////////////////////////////////////////////////*/
    int processedCount = 0; // To track the number of processed vertices

    while (!q.empty()) {
        int u = q.front(); 
        q.pop();
        processedCount++; // Increment processed count
        
        for (auto it = adj[u].begin(); it != adj[u].end(); it++) {
            int v = *it;
            indegree[v]--;

            level[v] = max(level[v], level[u] + 1);

            if (indegree[v] == 0) {
                q.push(v);
                processed[v] = true; // Mark as processed
            }
        }
    }

    // Check if the topological sort was successful
    if (processedCount != nCount) {
        cerr << "Error: Graph contains a cycle or is not connected." << endl;
        // throw runtime_error("Topological sorting failed due to cycle or unprocessed vertices.");
    }

    for (int i = 0; i < nCount; i++) {
        if (!processed[i]) {
            cerr << "Vertex " << i << " was not processed. Indegree: " << indegree[i] << endl;
        }
    }

    cout << "Topological Sort completed successfully." << endl;
    // cout << "Levels of vertices: ";
    // for (int i = 0; i < nCount; i++) {
    //     cout << level[i] << " ";
    // }
    // cout << endl;
}


/*
void Graph::topologicalSort()
{
	int *indegree=new int[nCount];
	// level = new int[nCount];

	//iterate over all nodes
	for(int i=0;i<nCount;i++)
	{
	    indegree[i]=0;
	    level.push_back(-1);
	}

	for(int i=0;i<nCount;i++)
	{
	    for(auto it=adj[i].begin();it!=adj[i].end();it++)
	    {
	        indegree[*it]++;
	    }    
	}

	cout<<"Indegree Caluculated Sucessfully"<<endl;

	for(int i=0;i<nCount;i++)
	{
	    cout<<indegree[i]<<" ";    
	}
	cout<<endl;

	queue<int> q;
	for(int i=0;i<nCount;i++)
	{
	    if(indegree[i]==0)
	    {
	        q.push(i);
	        level[i] = 0;
	    }
	}

	int u;

	while(!q.empty())
	{
	    u= q.front();q.pop();
	    
	    for(auto it=adj[u].begin();it!=adj[u].end();it++)
	    {
	        int v = *it;
	        indegree[v]--;

	        level[v] = max(level[v], level[u]+1);

	        if(indegree[v]==0)
	        {
	            q.push(v);
	        }
	    }
	}//while loop is over
	
}// topological sort is over
*/

void Graph::printLevelNumbers(ofstream &connectionsWithLevels)
{
    // cout<<"level numbers: "<< endl;
    /*for(int i=0; i<nCount; i++)
    {
        cout<<level[i]<< " ";
        // connectionsWithLevels<<
    }*/
    connectionsWithLevels<<vCount<<" "<<nCount<<" "<<totalDependencies<<" "<<*max_element(level.begin(), level.end())+1<<endl;

    
    // for(int i=0;i<=nCount;i++)
	// {
	// 	connectionsWithLevels<<offset[i]<<" ";
	// }
	// connectionsWithLevels<<endl;

	// for(int i=0;i<dependency.size();i++)
	// {
	// 	connectionsWithLevels<<dependency[i]<<" ";
	// }
	// connectionsWithLevels<<endl;

    vector<array<int, 5>> connectionsWithLevelsVector;
    for(int i=0; i<nCount; i++)
    {
        conn edge = connectionsVector[i];
        array<int, 5> connectionWithLevel = {edge.u, edge.v, edge.departure, edge.duration, level[i]};
        connectionsWithLevelsVector.push_back(connectionWithLevel);
    }
    sort(connectionsWithLevelsVector.begin(), connectionsWithLevelsVector.end(), [](const array<int, 5>& a, const array<int, 5>& b) {
        return a[4] < b[4]; // Sort by level
    });
    int i= 0;
    for(const auto& connection : connectionsWithLevelsVector) {
        connectionsWithLevels << i << " " << connection[0] << " " << connection[1] << " " << connection[2] << " " << connection[3] << " " << connection[4] << endl;
        i++;
    }


}

int main(int argc, char *argv[])
{
	ifstream connections; 
   	connections.open(argv[1]); 

   	int lvCount,lnCount, ltotalDependencies;

   	connections>>lvCount>>lnCount>>ltotalDependencies;

   	cout<<"In main"<<lvCount<<" "<<lnCount<<" "<<ltotalDependencies<<endl;


   	Graph g(lvCount,lnCount,ltotalDependencies);

   	g.readESDGFromFile(connections);

   	cout<<"read connections done"<<endl;

   	
	ofstream connectionsWithLevels;
   	connectionsWithLevels.open(argv[2]);
   	cout<<"outfile created "<<endl;  

   	/*if (g.isCyclic()) {
    	cout << "Error: The graph contains a cycle. Topological sort is not possible." << endl;
    	// return 1;
	}*/
 	
 	/*if (g.detectCycle()) {
        cout << "Cycle found!" << endl;
    } else {
        cout << "No cycle found." << endl;
    }*/

    g.topologicalSort();
    g.printLevelNumbers(connectionsWithLevels);

	return 0;
}