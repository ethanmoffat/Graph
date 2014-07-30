/*
 * Ethan Moffat
 * Dijkstra - Shortest Path First Lab
 * 03-05-2014
 * main.cpp
 */

#include <iostream>
#include "graph.h"

using namespace std;

//a zero represents no connection between vertices
static const int GRAPH_SIZE = 9;
const int graphData[GRAPH_SIZE][GRAPH_SIZE] = {
	{ 0, 4, 0, 0, 0, 0, 0, 8, 0 },
	{ 4, 0, 8, 0, 0, 0, 0, 11, 0 },
	{ 0, 8, 0, 7, 0, 4, 0, 0, 2 },
	{ 0, 0, 7, 0, 9, 14, 0, 0, 0 },
	{ 0, 0, 0, 9, 0, 10, 0, 0, 0 },
	{ 0, 0, 4, 14, 10, 0, 2, 0, 0 },
	{ 0, 0, 0, 0, 0, 2, 0, 1, 6 },
	{ 8, 11, 0, 0, 0, 0, 1, 0, 7 },
	{ 0, 0, 2, 0, 0, 0, 6, 7, 0 }
	}; //this is the representation of the graph from the assignment

const int slides[GRAPH_SIZE][GRAPH_SIZE] = {
	{ 0,1,0,0,0,1,0,0,1 },
	{ 1,0,1,0,1,0,0,0,0 },
	{ 0,1,0,1,1,0,0,0,0 },
	{ 0,0,1,0,0,0,1,1,0 },
	{ 0,1,1,0,0,0,1,0,0 },
	{ 1,0,0,0,0,0,1,0,0 },
	{ 0,0,0,1,1,1,0,0,0 },
	{ 0,0,0,1,0,0,0,0,0 },
	{ 1,0,0,0,0,0,0,0,0 }
	}; //this is a representation of the graph from the chapter 3 slides (for testing DFS/BFS)

int printVert(int vertex)
{ //generic print function. prints the node.
	cout << "Visited node: " << vertex << endl;
	return 0;
}

int printVertAlpha(int vertex)
{ //converts the node to a letter representation starting from 'a'.
	cout << "Visited node: " << (char)(vertex + 'a') << endl;
	return 0;
}

int printPlus1(int vertex)
{ //prints the node as a 1-based index (for MST testing)
	cout << "Visited node: " << vertex + 1 << endl;
	return 0;
}

void CreateMSTTest(Graph &g) //making an adjacency matrix for this would have been hell. use the methods in Graph class to generate it.
{// see http://www.me.utexas.edu/~jensen/exercises/mst_spt/mst_demo/mst1.html for more info
	if (g.NumVerts() != 11)
		g = Graph(11); //we have 11 nodes in the graph

	g.AddEdge(0,1,8); g.AddEdge(0,6,9); g.AddEdge(0,7,10); g.AddEdge(0,8,6); g.AddEdge(0,9,12); g.AddEdge(0,10,3);
	g.AddEdge(1,2,10); g.AddEdge(1,4,2); g.AddEdge(1,10,7);
	g.AddEdge(2,3,9); g.AddEdge(2,10,5);
	g.AddEdge(3,4,13); g.AddEdge(3,5,12);
	g.AddEdge(4,5,10); g.AddEdge(4,6,6);
	g.AddEdge(5,6,8);
	g.AddEdge(6,7,7);
	g.AddEdge(7,8,3);
	g.AddEdge(8,9,10);
	g.AddEdge(9,10,8);
}

int main(int c, char * argv[])
{
	Graph * g = NULL;
	g = Graph::CreateFromArray<GRAPH_SIZE>(graphData);

	for (int i = 0; i < g->NumVerts(); ++i)
	{ //run Dijkstra's algorithm from every starting vertex
		cout << "Running Dijkstra's algorithm starting from: " << i << endl;
		cout << left << setw(8) << "Vertex" << setw(10) << "Distance" << "Path" << endl;
		g->Dijkstra(i, [](int vert, int dist, list<int> pathInfo) //lambda expression for visit function
						{
							cout << setw(8) << left << vert << setw(10) << dist; //vertex/distance
							list<int>::iterator iter = pathInfo.begin();
							while (iter != pathInfo.end()) //iterate over the path list and print the values
							{
								cout << *iter;
								if (++iter != pathInfo.end())
									cout << "-->"; //separate with the arrow
							}
							cout << endl; //endline after the visit function is done
							return 0;
						});

		cout << endl << endl; //spacing between calls to g->Dijkstra
	}

	delete g; //the factory function CreateFromArray dynamically allocates memory.

	//The following is a demo of the MST algorithm
	//Graph mstStart(11), mst;
	//CreateMSTTest(mstStart);
	//mstStart.GetMinimumSpanningTree(mst);
	//cout << endl;
	//cout << "starting graph " << (mstStart.IsBipartite() ? "is" : "is not") << " bipartite." << endl;
	//cout << "MST      graph " << (mst.IsBipartite() ? "is" : "is not") << " bipartite." << endl;
	//cout << "Edge count: " << mst.NumEdges() << " | Vertex count: " << mst.NumVerts() << endl;
	//cout << mst;
	return 0;
}