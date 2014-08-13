/*
 * ETHAN MOFFAT
 * Dijkstra shortest path first lab
 * graph.h
 */

#ifndef GRAPH_DIJKSTRA_H
#define GRAPH_DIJKSTRA_H

#include <ostream>
#include <iomanip>
#include <functional>
#include <queue>
#include <stack>
#include <list>

/*
Considering organization as follows
	0  1  2
0 | 0  1  1
1 | 0  0  0
2 | 0  1  0

1--2
| /
3

So, 2-dim array where each element is type
struct gEntry
{
	int distance;
	T value;
};

With associated T * vertVals
[0]: 5
[1]: 7
[2]: 1
That contains values for the vertices (lookup table)

enhanced:
	use adjList[1][2]->distance to get the edge distance between vert 1 and vert 2 (instead of adjList[1][2])
new:
	use adjList[1][2]->value to get the value of the edge between vert 1 and vert 2 (templated)
	use vertVals[1] to get the value of vertex 1 (templated)
*/

class Graph
{
public:
	//(factory function) create undirected graph from 2D array
	//This is designed this way because otherwise I would have had to template the entire
	//	class with vertCount parameter. This way, I only template the static factory function
	//	method and we keep the class looking cleaner. It also provides a more verbose construction
	//	technique when creating graphs from adjacency matrices.
	template<int vertCount>
	static Graph* CreateFromArray(const int vertData[vertCount][vertCount]);

public: /* -- CONSTRUCTORS/DESTRUCTORS/ASSIGNMENT OPERATORS -- */
	//default constructor: create an empty, undirected graph
	Graph() : size(0), edgeCount(0), directed(false), adjList(NULL) { }
	//creates a graph with the specified number of vertices but no edges
	Graph(int numVerts, bool directed = false);
	//copy from other graph
	Graph(const Graph& rhs);
	//assignment operator (copy from right hand side)
	Graph& operator=(const Graph &rhs);
	~Graph(); //clean up memory

public: /* -- MODIFIER METHODS -- */
	//Adds a vertex to the graph (resizes internal array container). Returns the new vertex number.
	int AddVertex();
	//removes a vertex (resizes array) and any edges referencing it.
	void RemoveVertex(int vertex);
	//adds an edge. Adds a bidirectional edge to an undirected graph (A->B and B->A)
	void AddEdge(int vertexA, int vertexB, int weight = 1);
	//removes an edge. removes both A->B and B->A in an undirected graph
	void RemoveEdge(int vertexA, int vertexB);
	//Removes all the edges. Used in GetShortestPathTree and GetMinimumSpanningTree
	inline void RemoveAllEdges() //resets all vertices to have no edges
	{
		for (int i = 0; i < size; ++i)
		{
			memset(adjList[i], 0, sizeof(int)* size);
		}
		edgeCount = 0;
	}

public: /* -- TRAVERSAL METHODS -- */
	//General parameters for these: specify the starting vertex and a 'visit' function
	//The visit function is called on the vertex when the vertex is visited by the algorithm

	//performs Breadth-first traversal (queue) on vertices
	void BFS(const int startVert, const std::function<void(int)> &visit = [](int){}) const;
	//performs Depth-first traversal (stack) on vertices
	void DFS(const int startVert, const std::function<void(int)> &visit = [](int){}) const;

	//Calculates Dijkstra's shortest-path-first algorithm from a starting vertex
	//parameters for Dijkstra visit function are as follows: VertexNumber, DistancefromStart, list PathToVertex
	void Dijkstra(const int startVert, 
		const std::function<void(int, int, std::list<int>)> &visit = [](int, int, std::list<int>){}) const;
	//uses Dijkstra method to return a new graph containing a shortest path tree (returns as parameter)
	void GetShortestPathTree(const int startVert, Graph &retGraph) const;
	//Uses Prim's algorithm to find and return the minimum spanning tree graph. Returns by reference.
	void GetMinimumSpanningTree(Graph &retGraph) const;
	
public: /* -- OTHER OPERATORS -- */
	//the << operator overload does a lot of fancy stuff formatting wise
	friend std::ostream& operator<<(std::ostream &str, const Graph &g);
	//this is how the user can interface directly with graph data
	const int * const operator[](int vertex) const;

	/* -- ACCESSOR METHODS -- */
	//First three simply access a data member
	inline int NumVerts() const { return this->size; }
	inline bool IsDirectedGraph() const { return this->directed; }
	inline int NumEdges() const { return this->edgeCount; }
	//This method uses red-blue coloring to determine if the graph is bipartite
	bool IsBipartite() const;

private:
	int **adjList; //adjacency matrix of edges
	int size, edgeCount; //size is number of vertices, edgeCount is number of edges
	bool directed; //flag for whether or not the graph is directed; important only when AddEdge/RemoveEdge are called

	//helper function for the dfs traversal: recursively called
	void dfsHelper(const int startVert, bool * const visited, const std::function<void(int)> &visit) const;
};

//definition of static factory function
template<int vertCount>
static Graph* Graph::CreateFromArray(const int vertData[vertCount][vertCount])
{
	Graph * ret = new Graph(vertCount);
	for (int i = 0; i < vertCount; ++i)
	{
		//we only need to diagonally fill half the array, since the graph is undirected
		//	it will automatically add edge B,A when edge A,B is created
		for (int j = i; j < vertCount; ++j)
		{
			ret->AddEdge(i, j, vertData[i][j]);
		}
	}
	return ret;
}
#endif