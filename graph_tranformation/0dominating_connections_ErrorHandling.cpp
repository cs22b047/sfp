/*
g++ -o dom.out 0dominating_connections.cpp
./dom.out connections_input.txt connections_sorted_t_output.txt
./dom.out ../../../Data/conn_sorted_t/connections_sample1.txt ../../../Data/dominated_conn_sort_t/sample1_connections_dom.txt

*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <queue>
#include <algorithm>

using namespace std;

typedef struct details
{
	int departure;
	int duration;

	bool operator()(details a,details b)
	{
		if(a.departure!=b.departure)
			return a.departure<b.departure;
		else 
			return a.duration>b.duration;
	}
}details;


typedef struct connection
{
	int u,v,departure,duration;

	bool operator()(connection a,connection b)
	{
		if(a.departure!=b.departure)
			return a.departure<b.departure;
		else 
			if(a.duration!=b.duration)
				return a.duration<b.duration;
			else
				return (a.v == b.u);
	}

}connection;


unordered_map<string,priority_queue<details,vector<details>,details>> getMap(ifstream &oldConnection,int totalVertices,int totalConnections)
{
	int u,v,t,lambda;

	unordered_map<string,priority_queue<details,vector<details>,details>> m;

	for(int i=0;i<totalConnections;i++)
	{
		oldConnection>>u>>v>>t>>lambda;

		details d;
		d.departure = t;
		d.duration = lambda;

		string key = to_string(u)+" "+to_string(v);
		auto it = m.find(key);
		
		if(it!=m.end())
		{
			it->second.push(d);
		}
		else
		{
			priority_queue<details,vector<details>,details> q;
			q.push(d);
			m.insert({key,q});
		}
	}

	
/*
	for(auto it=m.begin(); it!=m.end();it++)
	{
		cout<<"connections between "+it->first<<endl;
		while(!it->second.empty())
		{
			details d;
			d = it->second.top();
			it->second.pop();
			cout<<d.departure<<" "<<d.duration<<endl;
		}
		cout<<endl;
	}*/

	return m;
}

void getNewConnections(ifstream &oldConnection,ofstream &newConnection)
{
	int totalVertices,totalConnections;
	oldConnection>>totalVertices;
	oldConnection>>totalConnections;
	/*cout<<"Enter totalVertices: "<<endl;
	cin>>totalVertices;
	cout<<" Enter totalConnections "<<endl;
	cin>>totalConnections;*/
	// totalVertices=2601977;
	// totalConnections=63497050;
	// totalVertices=986;
	// totalConnections=332334;

	unordered_map<string,priority_queue<details,vector<details>,details>> m;
	m = getMap(oldConnection,totalVertices,totalConnections);

	//verified the content of m

	vector<connection> nonDominatedConnections;
	int u,v;
	string key;


	for(auto it = m.begin();it!=m.end();it++)
	{
		connection c;
		details currDetails;

		currDetails = it->second.top(); it->second.pop();
		key = it->first;

		stringstream keyStream(key);
		keyStream>>u>>v;

		c.u = u;
		c.v = v;
		c.departure = currDetails.departure;
		c.duration = currDetails.duration;

		//cout<<"keystream u "<<u<<" v "<<v<<endl;

		nonDominatedConnections.push_back(c);

		while(!it->second.empty())
		{

			currDetails = it->second.top(); it->second.pop();
			if((currDetails.departure+currDetails.duration) < (nonDominatedConnections.back().departure +nonDominatedConnections.back().duration))
			{
				c.u = u;
				c.v = v;
				c.departure = currDetails.departure;
				c.duration = currDetails.duration;

				nonDominatedConnections.push_back(c);
			}

		}



	}

	sort(nonDominatedConnections.begin(), nonDominatedConnections.end(), connection());

	newConnection<<totalVertices<<" "<<nonDominatedConnections.size()<<endl;
	for(int i=0;i<nonDominatedConnections.size();i++)
	{
		connection curr = nonDominatedConnections[i];
		newConnection<<i<<" "<<curr.u<<" "<<curr.v<<" "<<curr.departure<<" "<<curr.duration<<endl;
	}
	
}


int main( int argc, char *argv[])
{
	ifstream oldConnection; 
   	oldConnection.open(argv[1]); 
	
	ofstream newConnection;
   	newConnection.open(argv[2]);

   	getNewConnections(oldConnection,newConnection);
}