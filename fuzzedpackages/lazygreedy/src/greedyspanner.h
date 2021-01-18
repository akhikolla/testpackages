#ifndef GREEDYSPANNER_H
#define GREEDYSPANNER_H

#include "structures.h"
#include "binaryheap.h"
#include <cmath>

inline double distance(const vertex &a, const double x, const double y)
{
		return hypot(a.x-x,a.y-y);
}

inline double distance(const vertex &a, const vertex &b)
{
		return distance(a, b.x, b.y);
}

template<class Heap>
void DoDijkstra(unsigned int j,
				const pointset& vertices,
				std::vector< std::pair<unsigned int, double> >* myEdges,
				Heap &myHeap,
				double t,
				std::vector<int>& QContent,
				Heap& Q,
				bool* dirtyBits) {
	unsigned int N = vertices.size();
	QContent[j] = -1;
	dirtyBits[j] = false;
	myHeap.clear();
	for (unsigned int i = 0; i < N; i++)
		myHeap.insert(i, std::numeric_limits<double>::infinity());
	myHeap.decreaseKey(j, 0);
	while (myHeap.getCount() > 0) {
		std::pair<unsigned int, double> pair = myHeap.getMin();
		double directDist = distance(vertices[pair.first], vertices[j]);
		if (pair.second > t * directDist)
		{
			if (vertices[pair.first].x >= vertices[j].x) {
				if (QContent[j] == -1 || directDist < distance(vertices[j], vertices[QContent[j]]))
					QContent[j] = pair.first;
			}
		}
		myHeap.extractMin();
		for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
			std::pair<unsigned int, double> edge = myEdges[pair.first][i];
			double alt = pair.second + edge.second;
			if (alt < myHeap.getValue(edge.first))
				myHeap.decreaseKey(edge.first, alt);
		}
	}
	if (QContent[j] != -1) {
		if (!Q.contains(j))
			Q.insert(j, distance(vertices[QContent[j]], vertices[j]));
		else {
			if (Q.getValue(j) < distance(vertices[QContent[j]], vertices[j]))
				Q.increaseKey(j, distance(vertices[QContent[j]], vertices[j]));
			else
				Q.decreaseKey(j, distance(vertices[QContent[j]], vertices[j]));
		}
	}
	else if (Q.contains(j))
		Q.remove(j);
}

template<class Heap>
int GreedyLinspace3(const pointset &vertices, double t, edgelist &edges) {
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	t = (t < 2.0) ? t : 2.0;
	Heap myHeap(N, 0); //

	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	std::vector<int> QContent(N, -1);
	Heap Q(N, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[N];
	for (unsigned int j = 0; j < N; j++)
		DoDijkstra<Heap>(j, vertices, myEdges, myHeap, t, QContent, Q, dirtyBits);

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne)
		{
			minQ = Q.getMin();
			if (dirtyBits[minQ.first])
				DoDijkstra<Heap>(minQ.first, vertices, myEdges, myHeap, t, QContent, Q, dirtyBits);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		int minRepFrom = minQ.first;
		int minRepTo = QContent[minQ.first];
		double minRepDist = minQ.second;

		std::pair<unsigned int, double> pair(minRepFrom, minRepDist);
		myEdges[minRepTo].push_back(pair);
		pair.first = minRepTo;
		myEdges[minRepFrom].push_back(pair);

		edges.push_back(edge(minRepFrom,minRepTo));

		edgeCount++;

		for (int j = 0; j < N; j++)
			if (QContent[j] != -1)
				dirtyBits[j] = true;
	}
end:

	delete[] dirtyBits;
	delete[] myEdges;
	return edgeCount;
}

#endif
