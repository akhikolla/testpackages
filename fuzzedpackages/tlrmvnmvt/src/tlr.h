#ifndef TLR_H
#define TLR_H
#include <RcppEigen.h>
struct TLRNode
{
        Eigen::MatrixXd U;
        Eigen::MatrixXd V;
        int maxColNum;
        int crtColNum;

        TLRNode& operator=(TLRNode &&refNode){
                U = std::move(refNode.U);
                V = std::move(refNode.V);
                maxColNum = refNode.maxColNum;
                crtColNum = refNode.crtColNum;
                return *this;
        }

	TLRNode(const TLRNode &refNode){
		U = refNode.U;
		V = refNode.V;
		maxColNum = refNode.maxColNum;
		crtColNum = refNode.crtColNum;
	}

	TLRNode(){}
};
#endif
