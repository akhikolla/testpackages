#ifndef STRUCTURES_H
#define STRUCTURES_H

// #define nullptr NULL //tmp fix omdat ik hier nog een oude gcc heb :P

#include <vector>

struct vertex {
	double x,y;
	vertex() : x(0), y(0) {}
	vertex(double x, double y) : x(x), y(y) {}
	bool operator==(const vertex &v) const {return (v.x==x && v.y==y);}
};

struct edge {
	unsigned int x,y;
	edge() : x(0), y(0) {}
	edge(unsigned int x, unsigned int y) : x(x), y(y) {}
	bool operator==(const edge &e) const {return (e.x==x && e.y==y) || (e.x==y && e.y==x);}
    inline bool operator < (const edge &o) const {
        if (x != o.x) return x < o.x;
		else return y < o.y;
    }
};

typedef std::vector<vertex> pointset;
typedef std::vector<edge> edgelist;


#endif
