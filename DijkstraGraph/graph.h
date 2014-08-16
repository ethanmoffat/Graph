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
#include <map>

template<typename T, typename U>
class Graph
{
private:
	struct edge; //fwd declaration

public:
	//(factory function) create undirected graph from 2D array
	//This is designed this way because otherwise I would have had to template the entire
	//	class with vertCount parameter. This way, I only template the static factory function
	//	method and we keep the class looking cleaner. It also provides a more verbose construction
	//	technique when creating graphs from adjacency matrices.
	template<int vertCount>
	static Graph* CreateFromArray(const int vertData[vertCount][vertCount]);

public: //typedefs for visit functions
	typedef std::function<void(int, const T&, int, std::list<std::pair<int, T>>)> SPFVisitFunction;
	typedef std::function<void(int, int)> TraversalFunction;

public: /* -- CONSTRUCTORS/DESTRUCTORS/ASSIGNMENT OPERATORS -- */
	//default constructor: create an empty, undirected graph
	Graph() : size(0), edgeCount(0), directed(false), adjList(NULL), vertValues(NULL) { }
	//creates a graph with the specified number of vertices but no edges
	Graph(int numVerts, bool directed = false);
	//copy from other graph
	Graph(const Graph& rhs);
	//assignment operator (copy from right hand side)
	Graph& operator=(const Graph &rhs);
	~Graph(); //clean up memory

public: /* -- MODIFIER METHODS -- */
	//Adds a vertex to the graph (resizes internal array container). Returns the new vertex number.
	int AddVertex(const T * const value = NULL);
	//removes a vertex (resizes array) and any edges referencing it.
	void RemoveVertex(int vertex);
	//removes the first vertex (resizes array) with the value specified
	void RemoveVertex(const T& value);
	//adds an edge. Adds a bidirectional edge to an undirected graph (A->B and B->A)
	void AddEdge(int vertexA, int vertexB, int distance = 1, const U * const value = 0);
	//removes an edge. removes both A->B and B->A in an undirected graph
	void RemoveEdge(int vertexA, int vertexB);
	//Removes all the edges. Used in GetShortestPathTree and GetMinimumSpanningTree
	inline void RemoveAllEdges() //resets all vertices to have no edges
	{
		for (int i = 0; i < size; ++i)
		{
			memset(adjList[i], 0, sizeof(edge)* size);
		}
		edgeCount = 0;
	}

	//add values for edges/vertices (useful when using CreateFromArray factory function)
	void AddVertexValue(int vertValue, const T& value);
	void AddEdgeValue(int vertAIndex, int vertBIndex, const U& value);

public: /* -- TRAVERSAL METHODS -- */
	//General parameters for these: specify the starting vertex and a 'visit' function
	//The visit function is called on the vertex when the vertex is visited by the algorithm

	//performs Breadth-first traversal (queue) on vertices
	void BFS(const int startVert, const std::function<void(int, const T&)> &visit = [](int, T){}) const;
	//performs Depth-first traversal (stack) on vertices
	void DFS(const int startVert, const std::function<void(int, const T&)> &visit = [](int, T){}) const;

	//Calculates Dijkstra's shortest-path-first algorithm from a starting vertex to every other vertex in the graph
	//parameters for Dijkstra visit function are as follows: VertexNumber, DistancefromStart, list PathToVertex
	void Dijkstra(const int startVert, const SPFVisitFunction &visit = [](int, const T&, int, std::list<std::pair<int,T>>){}) const;
	//A* search implementation. A* is an algorithm which, unlike Dijkstra's algorithm, solves single-pair rather than single-source
	std::list<std::pair<int, T>> AStar(const int startVert, const int endVert) const;
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

	//the following structs are internal structures for managing graph data
private:
	T * vertValues; //this is a lookup table that is sizeof(numVerts). Each T is the value of it's corresponding vertex
	struct edge
	{
		int distance;
		U value;

		edge() : distance(0) { }
		inline edge& operator=(const edge &rhs)
		{
			if (this == &rhs)
				return *this;

			this->distance = rhs.distance;
			this->value = U(rhs.value);//invoke copy constructor
			
			return *this;
		}
	};

private:
	edge **adjList; //adjacency matrix of edges
	int size, edgeCount; //size is number of vertices, edgeCount is number of edges
	bool directed; //flag for whether or not the graph is directed; important only when AddEdge/RemoveEdge are called

	//helper function for the dfs traversal: recursively called
	void dfsHelper(const int startVert, bool * const visited, const std::function<void(int, const T&)> &visit) const;

	//helper function to reconstruct the path for the return value of A*
	std::list<std::pair<int, T>> reconstruct_path(const std::map<int, int> &came_from, int vert) const;
};

//definition of static factory function
template<typename T, typename U>
template<int vertCount>
static Graph<T, U>* Graph<T, U>::CreateFromArray(const int vertData[vertCount][vertCount])
{
	Graph<T, U> * ret = new Graph<T, U>(vertCount);
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

template<typename T, typename U>
Graph<T, U>::Graph(int numVerts, bool directed)
{
	this->size = numVerts;
	this->directed = directed;
	this->edgeCount = 0;

	if (size > 0)
	{
		this->adjList = new edge*[this->size];
		this->vertValues = new T[this->size];

		for (int i = 0; i < this->size; ++i)
		{
			this->adjList[i] = new edge[this->size];
			this->vertValues[i] = T();
			memset(this->adjList[i], 0, sizeof(edge)* this->size);
		}
	}
	else
	{
		adjList = NULL;
		vertValues = NULL;
	}
}

template<typename T, typename U>
Graph<T, U>::Graph(const Graph<T, U> &rhs)
{
	edgeCount = rhs.NumEdges();
	size = rhs.NumVerts();
	directed = rhs.IsDirectedGraph();

	if (size > 0)
	{
		adjList = new edge*[size];
		vertValues = new T[size];
		for (int i = 0; i < size; ++i)
		{
			adjList[i] = new edge[size];
			memcpy(adjList[i], rhs.adjList[i], sizeof(edge)* size);
			vertValues[i] = T(rhs.vertValues[i]);
		}
	}
	else
	{
		adjList = NULL;
		vertValues = NULL;
	}
}

template<typename T, typename U>
Graph<T, U>& Graph<T, U>::operator=(const Graph<T, U> &rhs)
{
	if (&rhs == this)
		return *this;

	if (adjList != NULL)
	{ //clear out this graph if there is data
		for (int i = 0; i < size; ++i)
			delete[] adjList[i];
		delete[] adjList;
		adjList = NULL;
	}
	if (vertValues != NULL)
	{
		delete[]vertValues;
		vertValues = NULL;
	}

	edgeCount = rhs.NumEdges();
	size = rhs.NumVerts();
	directed = rhs.IsDirectedGraph();

	if (size > 0)
	{
		adjList = new edge*[size];
		vertValues = new T[size];
		for (int i = 0; i < size; ++i)
		{
			adjList[i] = new edge[size];
			memcpy(adjList[i], rhs.adjList[i], sizeof(edge)* size);
			vertValues[i] = T(rhs.vertValues[i]);
		}
	}
	else
	{
		adjList = NULL;
		vertValues = NULL;
	}

	return *this;
}

template<typename T, typename U>
Graph<T, U>::~Graph()
{
	if (adjList != NULL)
	{
		for (int i = 0; i < size; ++i)
		{
			//free each row
			delete[] adjList[i];
			adjList[i] = NULL;
		}
		//free the whole thing
		delete[] adjList;
		adjList = NULL;
	}

	if (vertValues != NULL)
	{
		delete[] vertValues;
		vertValues = NULL;
	}
}

template<typename T, typename U>
int Graph<T, U>::AddVertex(const T * const value)
{
	edge ** oldAdjList = NULL;
	T * oldVertValues = NULL;

	if (size > 0)
	{//pointers to old data
		oldAdjList = adjList; 
		oldVertValues = vertValues;
	}
	size++;
	adjList = new edge*[size]; //allocate new memory
	vertValues = new T[size];

	for (int i = 0; i < size; ++i)
	{
		adjList[i] = new edge[size]; //allocate each row
		memset(adjList[i], 0, sizeof(edge)* size); //set to 0
	}

	if (oldAdjList != NULL)
	{
		for (int i = 0; i < size - 1; ++i)
		{
			memcpy(adjList[i], oldAdjList[i], sizeof(edge)* (size - 1)); //copy the memory over
			delete[] oldAdjList[i]; //delete the old memory
			oldAdjList[i] = NULL;
		}
		delete[] oldAdjList; //free the whole of the old memory
		oldAdjList = NULL;
	}

	if (oldVertValues != NULL)
	{
		memcpy(vertValues, oldVertValues, sizeof(T)* (size - 1));
		delete[] oldVertValues;
		oldVertValues = NULL;
	}

	if (value != NULL)
	{
		vertValues[size - 1] = *value;
	}
	else
	{
		vertValues[size - 1] = T();
	}

	return size - 1; //new highest index is the old size
}

template<typename T, typename U>
void Graph<T, U>::RemoveVertex(int vertex)
{
	if (vertex < 0 || vertex >= size || size == 0)
		return;

	edge ** oldData = adjList;
	T * oldVertValues = vertValues;
	size--;

	if (size > 0)
	{
		adjList = new edge*[size];
		vertValues = new T[size];
		for (int src = 0, dst = 0; dst < size; ++dst, ++src)
		{ //foreach row: skip over vertex row, copy otherwise
			if (src == vertex)
			{//skip source row, but still copy over for this iteration (for dest)
				delete[]oldData[src];
				oldData[src] = NULL;
				src++;
			}
			adjList[dst] = new edge[size];

			for (int row = 0, off = 0; (row - off) < size; ++row)
			{ //for this row: skip over the vertex, copy otherwise
				if (row == vertex)
				{
					if (oldData[src][row].distance > 0)
						edgeCount--;
					off = 1;
					row++;
				}
				adjList[dst][row - off] = oldData[src][row];
			}
			delete[] oldData[src];
			oldData[src] = NULL;

			vertValues[dst] = T(oldVertValues[src]);
		}
		delete[] oldVertValues;
		oldVertValues = NULL;
	}
	else
	{ //if we're deleting the last vertex things are very easy
		delete[] adjList[0];
		delete[] adjList;
		delete[] vertValues;
		adjList = NULL;
		vertValues = NULL;
	}
}

template<typename T, typename U>
void Graph<T, U>::RemoveVertex(const T &value)
{
	for (int i = 0; i < size; ++i)
	{
		if (value == vertValues[i])
		{
			//This is an easy wrapper using a sequential search to find the
			//	first matching vertex value
			RemoveVertex(i);
			break;
		}
	}
}

template<typename T, typename U>
void Graph<T, U>::AddEdge(int vertexA, int vertexB, int distance, const U * const value)
{
	//The following condition implies that A and B must both be created vertices
	//	in order to add an edge between them. Adding an edge DOES NOT manipulate the 
	//	number of vertices in the graph.
	if (vertexA < 0 || vertexA >= size || //vertex A must be between 0 and size - 1
		vertexB < 0 || vertexB >= size || //vertex B must be between 0 and size - 1
		distance < 1 || //weight must be greater than 0 (adding an edge implies we are not removing it - a zero would remove it
		vertexA == vertexB) //don't allow edges to point from A->A (I don't want to support this case)
		return;

	if (adjList[vertexA][vertexB].distance == 0)
		edgeCount++;

	edge newEdge;
	newEdge.distance = distance;
	if (value != NULL)
		newEdge.value = *value;

	adjList[vertexA][vertexB] = newEdge;
	if (!directed)
		adjList[vertexB][vertexA] = newEdge;
}

template<typename T, typename U>
void Graph<T, U>::RemoveEdge(int vertexA, int vertexB)
{
	if (vertexA < 0 || vertexA >= size ||
		vertexB < 0 || vertexB >= size || 
		vertexA == vertexB || edgeCount == 0) //similar restraints for removing the edge
		return;

	adjList[vertexA][vertexB] = edge();
	edgeCount--;
	if (!directed)
		adjList[vertexB][vertexA] = edge();
}

template<typename T, typename U>
void Graph<T, U>::AddVertexValue(int vertIndex, const T& value)
{
	if (vertIndex < 0 || vertIndex >= size)
		return;

	vertValues[vertIndex] = value;
}

template<typename T, typename U>
void Graph<T, U>::AddEdgeValue(int vertAIndex, int vertBIndex, const U& value)
{
	if (vertAIndex < 0 || vertAIndex >= size || vertBIndex < 0 || vertBIndex >= size)
		return;

	adjList[vertAIndex][vertBIndex].value = value;
	if (!directed)
		adjList[vertBIndex][vertAIndex].value = value;
}

template<typename T, typename U>
void Graph<T, U>::BFS(const int startVert, const std::function<void(int, const T&)> &visit) const
{ //BFS traversal algorithm as described in slides
	if (startVert < 0 || startVert >= size || size == 0 || edgeCount == 0)
		return;

	std::queue<int> toVisit; //use queue to track vertices we want to visit
	bool * visited = new bool[size];
	memset(visited, 0, sizeof(bool)* size);

	toVisit.push(startVert);
	visited[startVert] = true;
	while (toVisit.size() > 0)
	{
		int curVert = toVisit.front();
		toVisit.pop();

		for (int i = 0; i < size; ++i)
		{
			if (adjList[curVert][i].distance > 0 && !visited[i])
			{
				toVisit.push(i);
				visited[i] = true;
			}
		}

		visit(curVert, vertValues[curVert]);
	}
	delete[] visited;
}

template<typename T, typename U>
void Graph<T, U>::dfsHelper(const int startVert, bool * const visited, const std::function<void(int, const T&)> &visit) const
{
	visited[startVert] = true;
	visit(startVert, vertValues[startVert]);

	for (int i = 0; i < size; ++i)
	{
		if (adjList[startVert][i].distance > 0 && !visited[i])
			dfsHelper(i, visited, visit);
	}
}

template<typename T, typename U>
void Graph<T, U>::DFS(const int startVert, const std::function<void(int, const T&)> &visit) const
{ //DFS traversal implemented as discussed in class
	if (startVert < 0 || startVert >= size || size == 0 || edgeCount == 0)
		return;

	bool * const visited = new bool[size];
	memset(visited, 0, sizeof(bool)* size);
	dfsHelper(startVert, visited, visit); //use the call stack as a stack container for DFS (recursion)
	delete[] visited;
}

template<typename T, typename U>
void Graph<T, U>::Dijkstra(const int startVert, const SPFVisitFunction &visit) const
{
	if (startVert < 0 || startVert >= size || size == 0 || edgeCount == 0)
		return;

	//set up distance array
	int * dist = new int[size];
	for (int i = 0; i < size; ++i) dist[i] = 0x7FFFFFFF;//distance of 0x7FFFFFFF is INFINITE, same as MAX_INT
	dist[startVert] = 0; //initial vertex has distance of 0
	
	bool * spt = new bool[size]; //spt is the list of vertices that are finalized. Set to true when done
	memset(spt, false, sizeof(bool)* size); //initially they are all false
	
	std::list<std::pair<int, T>> * pathList = new std::list<std::pair<int, T>>[size]; //create an array of lists: each list contains the path
	for (int i = 0; i < size; ++i)
	{
		pathList[i].push_back(pair<int, T>(startVert, vertValues[startVert])); //each list starts at startVert
	}

	for (int i = 0; i < size - 1; ++i) //iterate over all vertices
	{
		int curVert, min = 0x7FFFFFFF; //min is used only the the next loop to search sequentially for the vertex with the minimum distance
		for (int j = 0; j < size; ++j) //find the minimum unchecked vertex
		{
			if (!spt[j] && dist[j] <= min)
			{
				min = dist[j];
				curVert = j;
			}
		}

		//curVert always equals startVert on the first iteration through
		spt[curVert] = true;
		for (int j = 0; j < size; j++)
		{
			if (adjList[curVert][j].distance == 0) //j is not an adjacent vertex
				continue;

			if (!spt[j] && dist[curVert] != 0x7FFFFFFF && 
				dist[curVert] + adjList[curVert][j].distance < dist[j])
			{
				//update the distance and the paths for vertex J
				dist[j] = dist[curVert] + adjList[curVert][j].distance;
				
				pathList[j].clear();
				for (std::list<std::pair<int, T>>::iterator iter = pathList[curVert].begin(); iter != pathList[curVert].end(); ++iter)
					pathList[j].push_back(*iter);
				pathList[j].push_back(std::pair<int, T>(j, vertValues[j])); //push the current vertex
			}
		}
	}

	for (int i = 0; i < size; ++i) visit(i, vertValues[i], dist[i], pathList[i]); //visit each vertex
	delete[] dist;
	delete[] spt;
	delete[] pathList;
}

template<typename T, typename U>
std::list<std::pair<int, T>> Graph<T, U>::reconstruct_path(const std::map<int, int> &came_from, int curVert) const
{
	auto path = std::list<std::pair<int, T>>();
	//reconstruct path following this pseudo-code:
	/*
	function reconstruct_path(came_from, current_node)
		if current_node in came_from
			p : = reconstruct_path(came_from, came_from[current_node])
			return (p + current_node)
		else
			return current_node
			*/
	//path.push_back(std::pair<int, T>(startVert, vertValues[startVert]));
	return path;
}

template<typename T, typename U>
std::list<std::pair<int, T>> Graph<T, U>::AStar(const int startVert, const int endVert) const
{
	//this follows the implementation provided on the page http://en.wikipedia.org/wiki/A*_search_algorithm

	//the set of nodes already evaluated
	bool * closedSet = new bool[size];
	memset(closedSet, false, sizeof(bool)* size);
	
	//the set of nodes to be evaluated
	//(default less function will sort on first int value of pair, which is used to represent priority)
	std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>> openset;
	openset.push(std::pair<int, int>(size - 1, startVert));

	//helper array providing a method of tracking whether or not a vertex is in openset
	bool *inOpenSet = new bool[size];
	memset(inOpenSet, false, sizeof(bool)*size);
	inOpenSet[startVert] = true;

	int * g_score = new int[size];
	memset(g_score, 0, sizeof(int)*size);

	std::map<int, int> came_from;

	while (!openset.empty())
	{
		//next vertex is the top of the priority queue
		int nextVert = openset.top().second;
		openset.pop();
		inOpenSet[nextVert] = false;
		if (nextVert == endVert)
			break;

		for (int i = 0; i < size; ++i)
		{
			if (adjList[nextVert][i].distance && closedSet[i])
				continue;

			int heuristic_cost = adjList[nextVert][i].distance;
			int tentative_g_score = g_score[nextVert] + heuristic_cost;
			if (!inOpenSet[i] || tentative_g_score < g_score[i])
			{
				came_from[i] = nextVert;
				g_score[i] = tentative_g_score;

				if (!inOpenSet[i])
				{
					openset.push(std::pair<int, int>(g_score[i] + heuristic_cost, i));
					inOpenSet[i] = true;
				}
			}
		}
	}

	delete[] closedSet;
	delete[] g_score;

	return reconstruct_path(came_from, endVert);
}
#endif