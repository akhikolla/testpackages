#include <Rcpp.h>
#include <map>
#include <vector>
#include <string>

using namespace Rcpp;
using std::string;
using std::vector;
using std::map;


void findNames(vector<string>& name, List x);
CharacterVector names(List x);
void dendIndex(map<string, int>& dendNameInx,
               CharacterVector& dendNames, CharacterVector& newNames);
IntegerVector clusterIndex(List x, map<string, int>& dendNameInx);
void runDend(map< int, IntegerVector>& hierarchy, map< int, vector<int> >& subsets, 
             List x, 
             int& counter, int superset, NumericVector& height, int heightInx, 
             map<string, int>& dendNameInx);
List  dend2hier(List x, NumericVector height, CharacterVector newNames);


//' @title Create a hierarchy from a dendrogram
//' 
//' @description A function which can be called from R. It is creating a 
//' hierarchy from a dendrogram.
//' 
//' @param x A dendrogram S3 R object.
//' @param height A vector of heights at which nodes are grouped.
//' @param newNames Names of the variabels which should be part of the 
//' hierarchy.
//'
//' @keywords internal
// [[Rcpp::export]]
List dend2hier(List x, NumericVector height, CharacterVector newNames)
{
    CharacterVector dendNames = names(x);
    map<string, int> dendNameInx; 
    dendIndex(dendNameInx, dendNames, newNames);
    // initial cluster
    IntegerVector cluster = wrap(dendNameInx);
    cluster.attr("height") = height[0];
    vector<int> nullSubset;
    cluster.attr("subset") = nullSubset;
    // output 
    map<int, vector<int> > subsets;
    map<int, IntegerVector> out;
    out[0] = cluster;
    int counter = 1;
    // run through the dendrogram
    runDend(out, subsets, x, counter, 0, height, 0, dendNameInx);
    // add subset attr to clusters
    for (map<int, vector<int> >::iterator i = subsets.begin(); i != subsets.end(); i++) {
        out[i->first].attr("subset") = i->second;
    }
    return wrap(out);
}


// @title Creat a node of a hierarchy
// 
// @description A function which recursively is called to generate all nodes 
// in the hierarchy. Call only from within a C++ functions!
// 
// @param hierarchy A map to which the node is added.
// @param subsets A map of indixces of subset of a set.
// @param x A dendrogram S3 R object.
// @param counter A interger for the position of the node in the 
// hierarchy.
// @param superset A integer giving the position of the next higher node.
// @param height A vector of heights at which nodes are grouped.
// @param heightInx A integer of the heightInx of heights.
// @param dendNameInx Name index to add to node.
// 
// @keywords internal
void runDend(map< int, IntegerVector>& hierarchy, map< int, vector<int> >& subsets, 
             List x, int& counter, int superset, NumericVector& height, int heightInx, 
             map<string, int>& dendNameInx)
{
    RObject node;
    double nodeHeight = 0;
    int newHeightInx = 0;
    for (int i = 0; i < x.size(); ++i) {
        node = as<RObject>(x[i]);
        for (int j = 0; j < height.size(); ++j) {
            nodeHeight = node.attr("height");
            if (nodeHeight <= height[j])
                newHeightInx = j;
            else
                break;
        }
        if (newHeightInx > heightInx) {
            IntegerVector cluster;
            RObject superNode = as<RObject>(hierarchy[superset]);
            if (node.hasAttribute("leaf") && node.attr("leaf")) {
                // make cluster
                string name = as<string>(node.attr("label"));
                cluster = dendNameInx[name];
                cluster.attr("height") = height[newHeightInx];
                cluster.attr("superset") = superset +1;
                hierarchy[counter] = cluster;
                // modify super cluster
                subsets[superset].push_back(counter +1);
                ++ counter;
            }
            else {
                // make cluster
                cluster = clusterIndex(x[i], dendNameInx);
                cluster.attr("height") = height[newHeightInx];
                cluster.attr("superset") = superset +1;
                hierarchy[counter] = cluster;
                // modify super cluster
                subsets[superset].push_back(counter +1);
                int newSuperset = counter;
                ++ counter;
                runDend(hierarchy, subsets, x[i], counter, newSuperset, 
                        height, newHeightInx, dendNameInx);
            }
        }
        else if (!node.hasAttribute("leaf"))
            runDend(hierarchy, subsets, x[i], counter, superset, 
                    height, heightInx, dendNameInx);
    }
}


// @title dendrogram names in order to the new names
//
// @description A function which creates a nameIndex. Call only from within a 
// C++ functions!
//
// @param dendNameInx A map to which the dendrogram names index is written.
// @param dendNames The names of variables which are part of the dendrogram.
// @param newNames Names of the variabels which should be part of the 
// hierarchy.
// 
// @keywords internal
void dendIndex(map<string, int>& dendNameInx, CharacterVector& dendNames, 
               CharacterVector& newNames)
{
    for (int i = 0; i < dendNames.size(); ++i) {
        for (int j = 0; j < newNames.size(); ++j) {
            if (dendNames[i] == newNames[j]) {
                dendNameInx[as<string>(dendNames[i])] = j + 1;
                break;
            }
        }
    }
}


// @title Index of cluster
// 
// @description A function which creates a Index of clusters. Call only from 
// within a C++ functions!
//
// @param x A dendrogram S3 R object.
// @papam dendNameInx  A map of name indexes.
//
// @keywords internal
IntegerVector clusterIndex(List x, map<string, int>& dendNameInx)
{
    CharacterVector name = names(x);
    int n = name.size();
    IntegerVector out(n);
    for (int i = 0; i < n; ++i)
        out[i] = dendNameInx[as<string>(name[i])];
    return out.sort();
}


// @titile Names of cluster
// 
// @description A function which creates a vector of names.
// 
// @param x A dendrogram S3 R object.
//
// @keywords internal
CharacterVector names(List x)
{
    vector<string> out;
    findNames(out, x);
    return wrap(out);
}


// @title Find names of dendrogram
// 
// @description A function which finds names by recursively calling it self. 
// Call only from within a C++ functions!
//
// @param name A vector in which the names are written.
// @param x A dendrogram S3 R object.
//
// @keywords internal
void findNames(vector<string>& name, List x) 
{
    for (int i = 0; i < x.size(); ++i) {
        RObject node = as<RObject>(x[i]);
        if (node.hasAttribute("label"))
            name.push_back(as<string>(node.attr("label")));
        else
            findNames(name, x[i]);
    }
}
