#ifndef BSP_H
#define BSP_H

#include <iostream>
#include <cstddef>
#include <vector>    
#include <map>
#include <cfloat>    /* for DBL_MAX */
#include <cstring>	 /* for string */
#include <cmath>

#include "type.h"
#include "util.h"
#include "bspconst.h"

using std::fabs;

/* a proposal for splitting, used in binarySplit() */
struct splitProp {
	double lclnArea;
	double rclnArea;
	int nlc;
	int nrc;
};

struct LLAsplitProp {
	vector<vector<int> > desc_numpts;
	vector<vector<double> > desc_lnarea;
};


/* forward declare bspTree */
class bspTree; 

class bspNode {
public:      						// only make sense when the node is a leaf node
	bspNode* prev;					// previous leaf node
	bspNode* next;					// next leaf node
private:
	vector<uint> idx;				// indices of the points that the node contains
	matrix range;	
	vector<int> clabels;
private:

	bspNode* parent;				
	bspNode* leftChild;
	bspNode* rightChild;

	int nodeID;
	int depth;

	bool isleaf;					// initialized to be true
	bool splitFurthur;				// initialized to be true
				
	float lnArea;					// log of the region size
	float lnMass;					// log of the density mass i.e. theta in Lu (2013) 
	uint dim;	

	/* Caching sis intermediate results */
	vector<bool> tested;
	vector<double> lclnArea;
	vector<double> rclnArea;
	vector<int> nlc;
	vector<int> nrc;

	/* spliting history */			
	vector<int> splitdims;			// keep track of dimensions that have been split along
	vector<int> dimHist;
	vector<int> luHist;
	int spliton;

public: 
	/* Constructors */
	bspNode();
	bspNode(double dummyval);
	bspNode(const matrix &range); 	// for LLA
	bspNode(bspNode *parent_);
	bspNode(const matrix &data, 
			const vector<double> &mmax, 
			const vector<double> &mmin);

	/* Update Functions */
	void updatelnMass();
	void updatelnArea();
	void updateIdx(bspNode *lc, bspNode *rc, const matrix *data, vector<uint> &parentIdx);
	
	/* Get Functions */
	inline double 	   getlnMass()       const { return lnMass; }
	inline double 	   getlnArea()       const { return lnArea; }
	inline int 	       getDepth()		 const { return depth; }
	inline bool 	   getSplitFurthur() const { return splitFurthur; 		 }
	inline uint 	   getNumPts() 	     const { return idx.size();  }
	inline double 	   getDensity() 	 const { return exp(lnMass-lnArea); }
	inline bspNode*    getLeftChild() 	 const { return leftChild;   }
	inline bspNode*    getRightChild()   const { return rightChild;  }
	inline int 	       getNodeID()		 const { return nodeID; }
	inline int 		   getDim()	      	 const { return dim; }
	inline bool 	   getIsleaf() 	  	 const { return isleaf; }
	inline bspNode*    getNextLeaf() 	 const { return next; }
	inline int 	       getSpliton() 	 const { return spliton; }
	inline vector<int> getSplitDims() 	 const { return splitdims; }
	inline int 		   getIdxAt(int i)   const { return idx[i]; }
		   int 	       getParentID()	 const;
		   int 	       getLcID()		 const; 
		   int 	       getRcID()		 const;
	
	
	/* Set Functions */
    
    inline void setIsleaf		(bool flag) 	   { isleaf = flag; 	  }
	inline void setSplitFurther (const bool flag)  { splitFurthur = flag; }	
	inline void setSpliton      (int dim)          { spliton = dim; }	
	inline void setNodeID	    (const int id)	   { nodeID = id; }	
	inline void setDepth	    (const int dp)	   { depth = dp; }	    
		   void setRange		(const matrix &r)  { range = r; }
		   void setSplitDims	(const uint dim)   { splitdims[dim] = 1; }
		   void updateSplitHist (const uint d, const int lu);
		   void addChildren	    (bspNode *lc, bspNode *rc);    
	
	/* Discrepancy */ 
	uint getMaxGapIdx 	  (const bspTree &T, const uint nCut) const;
	void discrepancySplit (const uint dim, const uint ptr, 
					  	   const uint nCut, bspTree &T,
					  	   bspNode *lc, bspNode *rc, int totalLeaves, double theta);	
	bool within 	 (const vector<double> &data) const;
	 int nwithinchild(const matrix &lcrange, const matrix *data);
	 int nwithinchild(const matrix &lcrange, const matrix *data, bspNode *child);

	/* SIS */
	void binarySplit(int cutdim, bspTree &T, bspNode* lc, bspNode* rc, int totalNodes);
	void whatIfSplit(const int cutdim, const bspTree &T, splitProp &pp);
	void LLAwhatIfSplit(const int cutdim, const bspTree &T, LLAsplitProp &pp);
	 int LLAsampleBinaryCut(const int level, bspTree &T, const double lnqprev);

};


class bspTree {
private:
	matrix *data; 
	bspNode root;
	bspNode *head;
	uint dim;	
	uint nleaves;
	uint nzleaves;
	double alpha;
	int nnodes;
	/* output */
	string ofilename; 
public:
	/* bsp-means */
	matrix leafctr; 
public:
	/* Constructors */
	bspTree() {};
	bspTree(matrix& data_,
   		    const vector<double> &mmax,
			const vector<double> &mmin);

	/* Destructor */
	~bspTree();
	void deallocateTree(const bspNode *n);

	/* Get Functions */
	inline const uint 	getNleaves() 	   const { return nleaves; }
	inline const uint 	getNZleaves() 	   const { return nzleaves; }
	inline const uint 	getDim() 		   const { return dim; }	
	inline const double getDataAt(uint row, uint col) const { return (*data)[row][col]; }
	inline matrix*		getDataPtr() 	   const { return data; }
	inline const uint 	getNNodes() 	   const { return nnodes; }
	inline const uint 	getNumData() 	   const { return (*data).size(); }
	inline bspNode* 	getFirstleaf() 	   const { return head; }
	inline double 		getAlpha()		   const { return alpha; }
	
	void appendleaf(bspNode *node);
	void removeleaf(bspNode *leaf);

	/* output */
	inline void setOutputFilename(const string& ofilename_) { ofilename = ofilename_; };

	/* Discrepancy */
	void dsp(const uint nCut, const uint maxlevel, double theta = 1.0);	

	/* SIS */
	void lla(const int maxlevel, const int minpts = 50);
	double logBPscore();

	/* bsp-means */
	void CalculateLeafCenter();
};

#endif /* BSP_H */
