/*
* ETHAN MOFFAT
* Dijkstra shortest path first lab
* graph.h
*/

#include "graph.h"

//Graph::Graph(int numVerts, bool directed)
//{
//	//create the correct number of vertices
//	//don't add any edges though
//	this->size = numVerts;
//	this->directed = directed;
//	this->edgeCount = 0;
//	this->adjList = new int*[this->size];
//	for (int i = 0; i < this->size; ++i)
//	{
//		this->adjList[i] = new int[this->size];
//		memset(this->adjList[i], 0, sizeof(int)* this->size);
//	}
//}

//Graph::Graph(const Graph& rhs)
//{
//	edgeCount = rhs.NumEdges();
//	size = rhs.NumVerts();
//	directed = rhs.IsDirectedGraph();
//
//	if (size > 0)
//	{
//		adjList = new int*[size];
//		for (int i = 0; i < size; ++i)
//		{
//			adjList[i] = new int[size];
//			for (int j = 0; j < size; ++j)
//			{
//				adjList[i][j] = rhs.adjList[i][j];
//			}
//		}
//	}
//}

//Graph& Graph::operator=(const Graph &rhs)
//{
//	if (&rhs == this)
//		return *this;
//
//	if (adjList != NULL)
//	{ //clear out this graph if there is data
//		for (int i = 0; i < size; ++i)
//			delete[] adjList[i];
//		delete[] adjList;
//	}
//
//	edgeCount = rhs.NumEdges();
//	size = rhs.NumVerts();
//	directed = rhs.IsDirectedGraph();
//
//	if (size > 0)
//	{
//		adjList = new int*[size];
//		for (int i = 0; i < size; ++i)
//		{
//			adjList[i] = new int[size];
//			for (int j = 0; j < size; ++j)
//			{
//				adjList[i][j] = rhs.adjList[i][j];
//			}
//		}
//	}
//
//	return *this;
//}

//Graph::~Graph()
//{
//	for (int i = 0; i < size; ++i)
//	{
//		//free each row
//		delete[] adjList[i];
//	}
//	//free the whole thing
//	delete[] adjList;
//}

//int Graph::AddVertex()
//{
//	int ** oldData = NULL;
//	if(size > 0)
//		oldData = adjList; //pointer to old data
//	size++;
//	adjList = new int*[size]; //allocate new memory
//	for (int i = 0; i < size; ++i)
//	{
//		adjList[i] = new int[size]; //allocate each row
//		memset(adjList[i], 0, sizeof(int)* size); //set to 0
//	}
//	if (oldData != NULL)
//	{
//		for (int i = 0; i < size - 1; ++i)
//		{
//			memcpy(adjList[i], oldData[i], sizeof(int) * (size - 1)); //copy the memory over
//			delete[] oldData[i]; //delete the old memory
//		}
//		delete[] oldData; //free the whole of the old memory
//	}
//	return size - 1; //new highest index is the old size
//}
//
//void Graph::RemoveVertex(int vertex)
//{
//	if (vertex < 0 || vertex >= size || size == 0)
//		return;
//
//	int ** oldData = adjList;
//	size--;
//	if (size > 0)
//	{
//		adjList = new int*[size];
//		for (int src = 0, dst = 0; dst < size; ++dst, ++src)
//		{ //foreach row: skip over vertex row, copy otherwise
//			if (src == vertex)
//			{//skip source row, but still copy over for this iteration (for dest)
//				delete[]oldData[src];
//				src++;
//			}
//			adjList[dst] = new int[size];
//
//			for (int row = 0, off = 0; (row - off) < size; ++row)
//			{ //for this row: skip over the vertex, copy otherwise
//				if (row == vertex)
//				{
//					if (oldData[src][row] > 0)
//						edgeCount--;
//					off = 1;
//					row++;
//				}
//				adjList[dst][row - off] = oldData[src][row];
//			}
//			delete[] oldData[src];
//		}
//	}
//	else
//	{ //if we're deleting the last vertex things are very easy
//		delete[] oldData[0];
//		delete[] oldData;
//	}
//}
//
//void Graph::AddEdge(int vertexA, int vertexB, int weight)
//{
//	//The following condition implies that A and B must both be created vertices
//	//	in order to add an edge between them. Adding an edge DOES NOT manipulate the 
//	//	number of vertices in the graph.
//	if (vertexA < 0 || vertexA >= size || //vertex A must be between 0 and size - 1
//		vertexB < 0 || vertexB >= size || //vertex B must be between 0 and size - 1
//		weight < 1 || //weight must be greater than 0 (adding an edge implies we are not removing it - a zero would remove it
//		vertexA == vertexB) //don't allow edges to point from A->A (I don't want to support this case)
//		return;
//
//	if (adjList[vertexA][vertexB] == 0)
//		edgeCount++;
//
//	adjList[vertexA][vertexB] = weight;
//	if (!directed)
//		adjList[vertexB][vertexA] = weight;
//}
//
//void Graph::RemoveEdge(int vertexA, int vertexB)
//{
//	if (vertexA < 0 || vertexA >= size ||
//		vertexB < 0 || vertexB >= size || 
//		vertexA == vertexB || edgeCount == 0) //similar restraints for removing the edge
//		return;
//
//	adjList[vertexA][vertexB] = 0;
//	edgeCount--;
//	if (!directed)
//		adjList[vertexB][vertexA] = 0;
//}
//
//void Graph::BFS(const int startVert, const std::function<void(int)> &visit) const
//{ //BFS traversal algorithm as described in slides
//	if (startVert < 0 || startVert >= size || size == 0 || edgeCount == 0)
//		return;
//
//	std::queue<int> toVisit; //use queue to track vertices we want to visit
//	bool * visited = new bool[size];
//	memset(visited, 0, sizeof(bool)* size);
//
//	toVisit.push(startVert);
//	visited[startVert] = true;
//	while (toVisit.size() > 0)
//	{
//		int curVert = toVisit.front();
//		toVisit.pop();
//
//		for (int i = 0; i < size; ++i)
//		{
//			if (adjList[curVert][i] > 0 && !visited[i])
//			{
//				toVisit.push(i);
//				visited[i] = true;
//			}
//		}
//
//		visit(curVert);
//	}
//	delete[] visited;
//}
//
//void Graph::dfsHelper(const int startVert, bool * const visited, const std::function<void(int)> &visit) const
//{
//	visited[startVert] = true;
//	visit(startVert);
//
//	for (int i = 0; i < size; ++i)
//	{
//		if (adjList[startVert][i] > 0 && !visited[i])
//			dfsHelper(i, visited, visit);
//	}
//}
//
//void Graph::DFS(const int startVert, const std::function<void(int)> &visit) const
//{ //DFS traversal implemented as discussed in class
//	if (startVert < 0 || startVert >= size || size == 0 || edgeCount == 0)
//		return;
//
//	bool * const visited = new bool[size];
//	memset(visited, 0, sizeof(bool)* size);
//	dfsHelper(startVert, visited, visit); //the recursive traversal is the stack container
//	delete[] visited;
//}
//
//void Graph::Dijkstra(const int startVert, const std::function<void(int, int, std::list<int>)> &visit) const
//{
//	if (startVert < 0 || startVert >= size || size == 0 || edgeCount == 0)
//		return;
//
//	//set up distance array
//	int * dist = new int[size];
//	for (int i = 0; i < size; ++i) dist[i] = 0x7FFFFFFF;//distance of 0x7FFFFFFF is INFINITE, same as MAX_INT
//	dist[startVert] = 0; //initial vertex has distance of 0
//	
//	bool * spt = new bool[size]; //spt is the list of vertices that are finalized. Set to true when done
//	memset(spt, false, sizeof(bool)* size); //initially they are all false
//	
//	std::list<int> * pathList = new std::list<int>[size]; //create an array of lists: each list contains the path
//	for (int i = 0; i < size; ++i) pathList[i].push_back(startVert); //each list starts at startVert
//
//	for (int i = 0; i < size - 1; ++i) //iterate over all vertices
//	{
//		int curVert, min = 0x7FFFFFFF; //min is used only the the next loop to search sequentially for the vertex with the minimum distance
//		for (int j = 0; j < size; ++j) //find the minimum unchecked vertex
//		{
//			if (!spt[j] && dist[j] <= min)
//			{
//				min = dist[j];
//				curVert = j;
//			}
//		}
//
//		//curVert always equals startVert on the first iteration through
//		spt[curVert] = true;
//		for (int j = 0; j < size; j++)
//		{
//			if (adjList[curVert][j] == 0) //j is not an adjacent vertex
//				continue;
//
//			if (!spt[j] && dist[curVert] != 0x7FFFFFFF && 
//				dist[curVert] + adjList[curVert][j] < dist[j])
//			{
//				//update the distance and the paths for vertex J
//				dist[j] = dist[curVert] + adjList[curVert][j];
//				
//				pathList[j].clear();
//				for (std::list<int>::iterator iter = pathList[curVert].begin(); iter != pathList[curVert].end(); ++iter)
//					pathList[j].push_back(*iter);
//				pathList[j].push_back(j); //push the current vertex
//			}
//		}
//	}
//
//	for (int i = 0; i < size; ++i) visit(i, dist[i], pathList[i]); //visit each vertex
//	delete[] dist;
//	delete[] spt;
//	delete[] pathList;
//}
//
//void Graph::GetShortestPathTree(const int startVert, Graph &retGraph) const
//{
//	retGraph = Graph(this->size, this->directed); //copy from this object
//
//	//use the dijkstra algorithm method to create a tree based on the paths contained in each visited vertex
//	this->Dijkstra(startVert, 
//		//lambda format: [&] means local variables accessible by reference
//		[&](int vert, int dist, std::list<int> path)
//		{
//			std::list<int>::iterator iter = path.begin();
//			while (iter != path.end())
//			{
//				int prev = *iter;
//				++iter;
//				if (iter == path.end())
//					break;
//				retGraph.AddEdge(prev, *iter, adjList[prev][*iter]); //edge between the two vertices in the path
//				//the weight of the new edge is stored conveniently in the adjList of this object
//			}
//		}
//	);
//}
//
//void Graph::GetMinimumSpanningTree(Graph &retGraph) const
//{
//	retGraph = Graph(*this);
//	retGraph.RemoveAllEdges();
//
//	bool * checkedVerts = new bool[size]; //mark vertices that are checked
//	memset(checkedVerts, 0, sizeof(bool)* size);
//	checkedVerts[0] = true;
//
//	while (retGraph.NumEdges() < retGraph.NumVerts() - 1) //not sure if this is necessary
//	{
//		int minVertA, minVertB, minDist = INT_MAX; //store minimum A, minimum B, and the minimum distance
//		bool foundVert = false;
//		for (int i = 0; i < size; ++i)
//		{ //find the next vert that is in S1 (checkedVerts)
//			if (checkedVerts[i])
//			{
//				for (int j = 0; j < size; ++j)
//				{ //make sure it has an edge connecting to a vertex in S2
//					if (!checkedVerts[j] && adjList[i][j] > 0 && adjList[i][j] < minDist)
//					{ //haven't added vert J to tree yet, and it is the new minimum edge
//						foundVert = true;
//						minVertA = i;
//						minVertB = j;
//						minDist = adjList[i][j];
//					}
//				}
//			}
//		}
//		if (!foundVert) //all vertices have been checked, terminate
//			break;
//
//		retGraph.AddEdge(minVertA, minVertB, minDist); //add the found edge to the MST
//		checkedVerts[minVertB] = true; //mark B as checked (A is already marked)
//	}
//
//	delete[] checkedVerts;
//}
//
//std::ostream& operator<<(std::ostream &str, const Graph &g)
//{
//	//fancy formatting
//	str << "       ";
//	for (int i = 0; i < g.size; ++i)
//	{
//		str << std::left << std::setw(4) << i;
//	}
//	str << std::endl;
//	str << "-------";
//	for (int i = 0; i < g.size; ++i)
//	{
//		str << "----";
//	}
//	str << std::endl;
//
//	//print rownum : [ contents ] for each vertex
//	for (int i = 0; i < g.size; ++i)
//	{
//		str << std::setw(2) << i << " : [ ";
//		for (int j = 0; j < g.size; ++j)
//		{
//			str << std::setw(2) << g.adjList[i][j];
//			if (j != g.size - 1)
//			{
//				str << ", ";
//			}
//		}
//		str << " ]" << std::endl;
//	}
//	return str;
//}
//
//const int * const Graph::operator[](int vertex) const
//{
//	if (vertex < 0 || vertex >= size)
//		return NULL;
//
//	return adjList[vertex];
//}
//
//bool Graph::IsBipartite() const
//{
//	if (size == 0 || edgeCount == 0) //don't bother doing work if we don't need to
//		return false;
//
//	bool * blueVerts = new bool[size]; //blue and red vertex arrays
//	bool *  redVerts = new bool[size];
//	memset(blueVerts, false, sizeof(bool)* size);
//	memset(redVerts, false, sizeof(bool)* size);
//
//	bool colorBlue = false; //colorBlue is set to true when we're coloring a vertex BLUE
//
//	blueVerts[0] = true; //color vertex 0 blue
//
//	int markedVerts = 1; //we've colored a vertex
//	for (int curVertex = 0; curVertex < size; ++curVertex)
//	{ //iterate over all vertices
//		for (int i = 0; i < size; ++i)
//		{
//			if (adjList[curVertex][i]) //make sure there is an edge
//			{
//				if (!blueVerts[i] && !redVerts[i]) //it is an unmarked vertex
//				{
//					markedVerts++;
//					if (colorBlue)
//						blueVerts[i] = true;
//					else
//						redVerts[i] = true; //mark it
//				}
//				else if (blueVerts[i] && !colorBlue) //trying to color a blue vertex red: not bipartite
//				{
//					return false;
//				}
//				else if (redVerts[i] && colorBlue) //trying to color a red vertex blue: not bipartite
//				{
//					return false;
//				}
//			}
//		}
//		colorBlue = !colorBlue; //switch color for each vertex we're going at
//		if (markedVerts == size) //break the loop if we've marked every vertex (trying to mark twice causes problems)
//			break;
//	}
//
//	delete[] blueVerts;
//	delete[] redVerts;
//	
//	//if we get here we were able to mark each vertex without overlap. This means it is bipartite.
//	return true; 
//}