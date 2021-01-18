#include <time.h>
#include <climits>
#include "bsp.h"

bspNode::bspNode() : prev(NULL)
, next(NULL)
, parent(NULL)
, leftChild(NULL) 
, rightChild(NULL)
, isleaf(true)
, splitFurthur(true)
{};

bspNode::bspNode(const matrix &range_) : prev(NULL)
, next(NULL)
, parent(NULL)
, leftChild(NULL) 
, rightChild(NULL)
, isleaf(true)
, splitFurthur(true)
{
	this->dim = range_.size();
	this->range = range_;
	this->tested.resize(range_.size());
	for(uint i = 0; i < dim; i++) {
		tested[i] = false;
	}
	lclnArea.resize(dim);
	rclnArea.resize(dim);
	nlc.resize(dim);
	nrc.resize(dim);	
};

bspNode::bspNode(double dummyval) : prev(NULL)
, next(NULL)
, parent(NULL)
, leftChild(NULL) 
, rightChild(NULL)
, isleaf(true)
, splitFurthur(true)
{
	dim = dummyval;
};

bspNode::bspNode(bspNode *parent_) : prev(NULL)
, next(NULL)
, leftChild(NULL) 
, rightChild(NULL)
, isleaf(true)
, splitFurthur(true)
{	
	parent = parent_;
	dim = (parent_)->dim;
	splitdims = (parent_)->splitdims;
	dimHist = (parent_)->dimHist;
	luHist = (parent_)->luHist;
	tested.resize(dim);
	for(uint i = 0; i < dim; i++) {
		tested[i] = false;
	}
	lclnArea.resize(dim);
	rclnArea.resize(dim);
	nlc.resize(dim);
	nrc.resize(dim);
	spliton = -1;
}

/*  data: data set used to build the tree
 *  mmax
 *  mmin: range of the data along each dimension
 *  No longer assumes the data lies within the unit range
 *  range of the node is initialized correspondingly based on mmax and mmin
 *  Used only to create root node
 */
 bspNode::bspNode(const matrix &data, 
 	const vector<double> &mmax, 
 	const vector<double> &mmin): prev(NULL)
 , next(NULL)
 , parent(NULL)
 , leftChild(NULL)
 , rightChild(NULL)
 , isleaf(true)
 , splitFurthur(true)
 , lnMass(0)

 {
 	uint npts = data.size();
 	dim = data[0].size();

 	for (uint i = 0; i < dim; i++) {
 		vector<double> one_d_range(1, mmin[i]);
 		one_d_range.push_back(mmax[i]);
 		range.push_back(one_d_range);
 	}

 	lnArea = reclnArea(range);
 	idx.resize(npts);
 	clabels.resize(npts);
	
 	for(uint i = 0; i < npts; i++) {	// root includes all the data points
 		idx[i] = i;
 		clabels[i] = 0;
 	}

 	splitdims.resize(dim);
 	for(uint i = 0; i < dim; i++) {
 		splitdims[i] = 0;
 	}	
 	tested.resize(dim);
 	for(uint i = 0; i < dim; i++) {
 		tested[i] = false;
 	}
 	lclnArea.resize(dim);
 	rclnArea.resize(dim);
 	nlc.resize(dim);
 	nrc.resize(dim);
 	spliton = -1;
 	return;
 }

 void bspNode::updatelnMass() {
 	this->lnMass = this->parent->getlnMass() + log((LAPLACE + this->getNumPts())/(2*LAPLACE + this->parent->getNumPts()));
 }

 void bspNode::updatelnArea() {
 	double la = 0;
 	for(uint i = 0; i < this->dim; i++) {
 		la += log(range[i][1] - range[i][0]);
 	}
 	lnArea = la;
 }

 void bspNode::updateIdx(bspNode *lc, bspNode *rc, const matrix *data, vector<uint> &parentIdx) {
 	uint n = parentIdx.size();
 	for(uint i = 0; i < n; i++) {
 		int dataidx = parentIdx[i];
 		if (lc->within((*data)[dataidx])) {
 			lc->idx.push_back(dataidx);
 		} else {
 			rc->idx.push_back(dataidx);
 		}
 	}
 	lc->clabels.resize(idx.size());
 	for(uint i = 0; i < idx.size(); i++) {
 		lc->clabels[i] = 0;
 	}
 	rc->clabels.resize(idx.size());
 	for(uint i = 0; i < idx.size(); i++) {
 		rc->clabels[i] = 0;
 	}
 	return;
 }

 void bspNode::updateSplitHist(const uint d, const int lu) {
 	this->dimHist.push_back((int)d);
 	switch (lu) {
 		case 0: this->luHist.push_back(0); break;
 		case 1: this->luHist.push_back(1); break;
 	}
 	return;
 }

/* determines whether a data point lies within the node */
 bool bspNode::within(const vector<double> &data) const {
 	for(uint d = 0; d < data.size(); d++) {
 		if (data[d] < this->range[d][0]) {
 			return false;
 		} else if (data[d] > this->range[d][1]) {
 			return false;
 		}
 	}
 	return true;
 }


/* returns the index of the maximal discrepancy */
 uint bspNode::getMaxGapIdx(const bspTree &T, const uint nCut) const {
  	uint dim = T.getDim();
 	uint n = this->getNumPts();
 	vector<double> gap(dim*(nCut-1), 0);

 	double incre = 1.0/n;
 	for(uint j = 0; j < dim; j++) {
 		double g = (range[j][1] - range[j][0])/ nCut;
		if (g < 1e-8) {
			return -1; 
		}
 		vector<double> pCount(nCut, 0);
 		for(uint i = 0; i < n; i++) {
 			int index = std::min(floor(double((T.getDataAt(idx[i], j) 
							- range[j][0])) / g) + 1, double(nCut)) - 1;

	

 			pCount[index] += incre;
 		}

 		for(uint k = 0; k < nCut - 1; k++) {
 			gap[j*(nCut-1)+k] = fabs(vecPartialSum(pCount, k) - (k+1.0)/nCut);
 		}
 	}

 	uint maxGapIndex = whichMax(gap);
 	return maxGapIndex;
 }

 void bspNode::discrepancySplit(const uint dim, const uint ptr, 
							 	const uint nCut, bspTree &T,
							 	bspNode *lc, bspNode *rc, 
							 	int totalNodes, double theta){
 	double gap = (range[dim][1] - range[dim][0]) / (nCut+1);
 	matrix lrange = range;
 	matrix rrange = range;

 	lrange[dim][1] = gap * ptr + lrange[dim][0];
 	rrange[dim][0] = gap * ptr + rrange[dim][0];

 	lc->setRange(lrange); 
 	rc->setRange(rrange);	

 	lc->updatelnArea();
 	rc->updatelnArea();

 	lc->updateSplitHist(dim, 0);
 	rc->updateSplitHist(dim, 1);

 	this->updateIdx(lc, rc, T.getDataPtr(), this->idx);
 	lc->updatelnMass();
 	rc->updatelnMass();

 	if (lc->getNumPts() < MAX_NPTS) {
 	//if (lc->getNumPts() < MAX_NPTS && lc->computeDiscrepancy(T) < THETA ) {
 		lc->setSplitFurther(false);
 	} else {
 		lc->setSplitFurther(true);
 	}
	if (rc->getNumPts() < MAX_NPTS) {
 	//if (rc->getNumPts() < MAX_NPTS && rc->computeDiscrepancy(T) < THETA)  {
 		rc->setSplitFurther(false);
 	} else {
 		rc->setSplitFurther(true);
 	}

 	lc->setNodeID(totalNodes + 1);
 	rc->setNodeID(totalNodes + 2);

 	lc->setDepth(this->getDepth() + 1);
 	rc->setDepth(this->getDepth() + 2);


	this->addChildren(lc, rc);
	return;
}

/* add two children to the node */
void bspNode::addChildren(bspNode *lc, bspNode *rc) {
	this->leftChild = lc;
	this->rightChild = rc;
	return;
}

int bspNode::getParentID() const {
	if (parent == NULL) {
		return -1;
	}
	return parent->getNodeID();
}

int bspNode::getLcID() const {
	if (leftChild == NULL) {
		return -1;
	}
	return leftChild->getNodeID();
}

int bspNode::getRcID() const {
	if (rightChild == NULL) {
		return -1;
	}
	return rightChild->getNodeID();
}

bspTree::bspTree(matrix& data_,
	const vector<double> &mmax,
	const vector<double> &mmin)
: nleaves(1)
, nzleaves(1)
, alpha(ALPHA)
, nnodes(1)
{
	data = &data_;
	dim = data_[0].size();
	root = bspNode(data_, mmax, mmin);
	root.setNodeID(0);
	root.setDepth(0);
	head = &root;
}

bspTree::~bspTree(){
	this->deallocateTree(root.getLeftChild());
	this->deallocateTree(root.getRightChild());
}

/* A recursive function to deallocate the binary tree for the nodes */
void bspTree::deallocateTree(const bspNode *n) {
	if (n == NULL) return;
	deallocateTree(n->getLeftChild());
	deallocateTree(n->getRightChild());
	delete n;
}

/* prepend node to the front */
void bspTree::appendleaf(bspNode *node){
	if (head == NULL) {
		head = node;
	} else {
		head->prev = node;
		node->next = head;
		head = node;
	}
	if (node->getNumPts() > 0) {
		nzleaves++;
	}
	nleaves++;
	nnodes++;
	return;
}

void bspTree::removeleaf(bspNode *leaf){
	if (leaf->prev != NULL) {
		leaf->prev->next = leaf->next;
	} else {
		head = leaf->next;
	}

	if (leaf->next != NULL) {
		leaf->next->prev = leaf->prev;
	} 

	leaf->prev = NULL;
	leaf->next = NULL;
	leaf->setIsleaf(false);
	if (leaf->getNumPts() > 0) {
		nzleaves--;
	}
	nleaves--;
	return;
}

/* quite version, only partition the space */
void bspTree::dsp(const uint nCut, const uint maxlevel, double theta) {
	while (this->getNZleaves() < maxlevel) {
		bspNode *leafptr = this->head;
		bool added = false;

		while (leafptr != NULL) {
			if (this->getNZleaves() >= maxlevel) break;
			if (leafptr->getSplitFurthur()) {	

				int maxIndex = leafptr->getMaxGapIdx(*this, nCut);
				
				if (maxIndex == -1) {
					leafptr->setSplitFurther(false);
					leafptr = leafptr->getNextLeaf();
					continue;
				}
				uint cutDim = floor(double(maxIndex/(nCut-1)));
				uint cutPtr = maxIndex % (nCut-1) + 1;
				leafptr->setSplitDims(cutDim);
				// creating two new nodes on the heap
				bspNode *lc = new bspNode(leafptr);
				bspNode *rc = new bspNode(leafptr);
				leafptr->discrepancySplit(cutDim, cutPtr, nCut,
										  *this, lc, rc,
										  this->getNNodes()-1, theta);
				leafptr->setSpliton(cutDim);
				this->appendleaf(rc);
				this->appendleaf(lc);
				vector<int> sd;
				sd = leafptr->getSplitDims();
				bspNode *next = leafptr->getNextLeaf();
				this->removeleaf(leafptr);
				leafptr = next;
				added = true;
			} else {
				leafptr = leafptr->getNextLeaf();
			}

		}
		if (!added)	break;

	}
	return;
}

int bspNode::nwithinchild(const matrix &lcrange, const matrix *data) {
	int count = 0;
	int n = (this->idx).size();
	for(int i = 0; i < n; i++){
		int dataidx = (this->idx)[i];
		if (::within(lcrange, (*data)[dataidx]) ) {
			count++;
		}
	} 
	return count;
}

int bspNode::nwithinchild(const matrix &lcrange, const matrix *data, bspNode *child) {
	int count = 0;
	int n = (this->idx).size();
	for(int i = 0; i < n; i++){
		int dataidx = (this->idx)[i];
		if (::within(lcrange, (*data)[dataidx]) ) {
			child->idx.push_back(dataidx);
			count++;
		}
	}
	
	return count;
}

/* Compute the number of points and log region size of the two children
 * The result is stored in splitProp
 */
 void bspNode::whatIfSplit(const int cutdim, const bspTree &T, splitProp &pp){
 	if (!tested[cutdim]) {
 		int n = this->idx.size();
 		matrix lcrange = this->range;
 		matrix rcrange = this->range;
 		
 		double cutptr = 0.5*(this->range[cutdim][1] + this->range[cutdim][0]);
 		lcrange[cutdim][1] = cutptr;
 		rcrange[cutdim][0] = cutptr;
 		pp.lclnArea = reclnArea(lcrange);
 		pp.rclnArea = reclnArea(rcrange);

 		pp.nlc = this->nwithinchild(lcrange, T.getDataPtr());
 		pp.nrc = n - pp.nlc;
 		this->lclnArea[cutdim] = pp.lclnArea;
 		this->rclnArea[cutdim] = pp.rclnArea;
 		this->nlc[cutdim] = pp.nlc;
 		this->nrc[cutdim] = pp.nrc;
 		tested[cutdim] = true;
 	} else {
 		pp.lclnArea = this->lclnArea[cutdim];
 		pp.rclnArea = this->rclnArea[cutdim];
 		pp.nlc = this->nlc[cutdim];
 		pp.nrc = this->nrc[cutdim];
 		
 	}
 	return;
}

 int bspNode::LLAsampleBinaryCut(const int level, bspTree &T, const double lnqprev) {

 	int dim = this->getDim();
 	vector<double> lnh(dim, 0);
 	vector<double> p(dim);
 	int N = T.getNumData();

 	for(int d = 0; d < dim; d++) {
 		LLAsplitProp pp;
 		this->LLAwhatIfSplit(d, T, pp);	
 		int n = this->getNumPts();
 		vector<double> des_score(2*dim);
 		for (int i = 0; i < 2*dim; i++) {
 			double alpha = T.getAlpha();
 			double u1 = n*this->getlnArea() - pp.desc_lnarea[i][0] -  pp.desc_lnarea[i][1] -  pp.desc_lnarea[i][2]; 
 			double u2 = lgamma(pp.desc_numpts[i][0] + alpha) + lgamma(pp.desc_numpts[i][1] + alpha) + lgamma(pp.desc_numpts[i][2] + alpha) - lgamma(n + alpha);
 			double u3 = lgamma(N + (level - 1)*alpha) - lgamma(N + (level+1)*alpha);
 			double u4 = lgamma((level+1)*alpha) - lgamma((level-1)*alpha) - lgamma(alpha)*2;
 			des_score[i] = lnqprev + u1 + u2 + u3 + u4;

 		}

 		double maxscore = vecmax(des_score); 
 		for(uint i = 0; i < des_score.size(); i++) {
 			des_score[i] = des_score[i] - maxscore;
 			lnh[d] += des_score[i];
 		}
 	}

	double maxlnh = vecmax(lnh); 
	for(uint i = 0; i < lnh.size(); i++) {
		lnh[i] = lnh[i] - maxlnh;
	}

	double sumlnh = 0;
	for(uint i = 0; i < lnh.size(); i++) {
		double tmp = exp(lnh[i]);
		p[i] = tmp;
		sumlnh += tmp;
	}
	
	for(uint i = 0; i < lnh.size(); i++) {
		p[i] /= sumlnh;
	}

	int cutDim = randsample(0, dim, p); 

	return cutDim;
}

void bspNode::LLAwhatIfSplit(const int cutdim, const bspTree &T, LLAsplitProp &pp) {
	int D = T.getDim();
	pp.desc_numpts.resize(D*2);					// for one step look ahead, this should be 2d x 3
	pp.desc_lnarea.resize(D*2);

	matrix lcrange = this->range;
	matrix rcrange = this->range;
	double cutptr = 0.5*(this->range[cutdim][1] + this->range[cutdim][0]);
	lcrange[cutdim][1] = cutptr + 1e-5;
	rcrange[cutdim][0] = cutptr + 1e-5;

	bspNode lcdummy(lcrange);	
	bspNode rcdummy(rcrange);

	double lclnArea = reclnArea(lcrange);
	double rclnArea = reclnArea(rcrange);

	int nlc = this->nwithinchild(lcrange, T.getDataPtr(), &lcdummy);
 	int nrc = this->nwithinchild(rcrange, T.getDataPtr(), &rcdummy); 
 	
 	for(int i = 0; i < D; i++) {
 		splitProp lcpp;
 		pp.desc_numpts[i].resize(3);
 		pp.desc_lnarea[i].resize(3);
 		lcdummy.whatIfSplit(i, T, lcpp);

 		pp.desc_numpts[i][0] = nrc;
 		pp.desc_numpts[i][1] = lcpp.nlc;
 		pp.desc_numpts[i][2] = lcpp.nrc;

 		pp.desc_lnarea[i][0] = rclnArea;
 		pp.desc_lnarea[i][1] = lcpp.lclnArea;
 		pp.desc_lnarea[i][2] = lcpp.rclnArea;

 		splitProp rcpp;
 		pp.desc_numpts[i + D].resize(3);
 		pp.desc_lnarea[i + D].resize(3);
 		rcdummy.whatIfSplit(i, T, rcpp);
 		pp.desc_numpts[i + D][0] = nlc;
 		pp.desc_numpts[i + D][1] = rcpp.nlc;
 		pp.desc_numpts[i + D][2] = rcpp.nrc;

 		pp.desc_lnarea[i + D][0] = lclnArea;
 		pp.desc_lnarea[i + D][1] = rcpp.lclnArea;
 		pp.desc_lnarea[i + D][2] = rcpp.rclnArea;
 	}	
 	return;
 }

void bspTree::lla(const int maxlevel, const int minpts) {

	/* Initialization */
	double lnqprev = 0;
	int level = 1;
	
	while (level < maxlevel) {
		
		bool added = false;
		bspNode *leafptr = this->head;
		while(leafptr != NULL) {
			if ((int)this->getNleaves() >= maxlevel) break;

			if (leafptr->getNumPts() > (uint)minpts) {	
				level++;
				
				int cutDim = leafptr->LLAsampleBinaryCut(level, *this, lnqprev); // TODO
			 	leafptr->setSplitDims(cutDim);
				leafptr->setSpliton(cutDim);
				bspNode *lc = new bspNode(leafptr);
				bspNode *rc = new bspNode(leafptr);

				leafptr->binarySplit(cutDim, *this, lc, rc, this->getNNodes());
				this->appendleaf(rc);
				this->appendleaf(lc);

				bspNode *next = leafptr->getNextLeaf();

				this->removeleaf(leafptr);
				leafptr = next;
				added = true;
				lnqprev = logBPscore();
			} else {
				leafptr = leafptr->getNextLeaf();
			}
			
		}

		if (!added) {
			break;
		} 

	}
	return;
}

double bspTree::logBPscore() {
	bspNode *leafptr = this->head;
	double score = 0;
	int NL = this->getNleaves();
	int N = 0;
	while (leafptr != NULL) {
		score += lgamma(leafptr->getNumPts() + this->alpha);
		score -= leafptr->getNumPts()*leafptr->getlnArea();
		N += leafptr->getNumPts();
		leafptr = leafptr->getNextLeaf();

	}
	score -= lgamma(NL * this->alpha + N);
	score -= NL*lgamma(this->alpha);
	score -= lgamma(NL*this->alpha);

	return score;
}

void bspNode::binarySplit(int cutdim, bspTree &T, bspNode* lc, bspNode* rc, int totalNodes) {

	this->setSplitDims(cutdim);
	matrix lcrange = this->range;
	matrix rcrange = this->range;

	double cutptr = 0.5*(this->range[cutdim][1] + this->range[cutdim][0]);

	lcrange[cutdim][1] = cutptr;
	rcrange[cutdim][0] = cutptr;

	lc->setRange(lcrange);
	rc->setRange(rcrange);

	lc->updatelnArea();
	rc->updatelnArea();

	lc->updateSplitHist(cutdim, 0);
	rc->updateSplitHist(cutdim, 1);

	this->updateIdx(lc, rc, T.getDataPtr(), this->idx);

	lc->updatelnMass();
	rc->updatelnMass();

	this->addChildren(lc, rc);

	if (lc->getNumPts() < 200) {
		lc->setSplitFurther(false);
	}

	if (rc->getNumPts() < 200) {
		rc->setSplitFurther(false);
	}

	lc->setNodeID(totalNodes );
	rc->setNodeID(totalNodes + 1);

	lc->setDepth(this->getDepth() + 1);
	rc->setDepth(this->getDepth() + 2);

	return;
}

void bspTree::CalculateLeafCenter() {
	int nl = this->getNZleaves();
	leafctr.resize(nl);

	bspNode *leafptr = this->getFirstleaf();
	int idx = 0;
	while (leafptr != NULL) {
		if (leafptr->getNumPts() == 0) {
			leafptr = leafptr->getNextLeaf();
			continue;
		}
		leafctr[idx].resize(this->getDim());
		for (uint i = 0; i < leafptr->getNumPts(); i++) {
			int data_idx = leafptr->getIdxAt(i);
			for(uint d = 0; d < this->getDim(); d++) {
				leafctr[idx][d] += this->getDataAt(data_idx, d) / leafptr->getNumPts();
			}
		}
		idx++;
		leafptr = leafptr->getNextLeaf();
	}

}

