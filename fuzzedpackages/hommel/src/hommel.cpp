#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
using std::vector;


// Implementation of Fortune 1989
std::vector<int> findhull (
    int             m,              // length of p
    vector<double>  &p)             // p-values (sorted!)
{
  // intialize output length
  int r;
  vector<int> hull(1);
  hull.push_back(1);
  bool notconvex;
  
  // find the hull
  for (int i = 2; i <= m; i++)
  {
    if (i == m || (m-1) * (p[i-1] - p[0]) < (i-1) * (p[m-1] - p[0])) {
      do {
        r = hull.size()-1;
        if (r>1) {
          notconvex = ((i-hull[r-1]) * (p[hull[r]-1] - p[hull[r-1]-1]) >=
            (hull[r] - hull[r-1]) * (p[i-1] - p[hull[r-1]-1])); 
        } else {
          if (r == 1) {
            notconvex = (i * p[hull[1]-1] >= hull[1] * p[i-1]); 
          } else {
            notconvex = false;
          }
        }
        if (notconvex) hull.resize(r);
      } while (notconvex);
      hull.push_back(i);
    }
  }
  
  return hull;
}

// Finds the jumps of h(alpha)
// [[Rcpp::export]]
std::vector<double> findalpha( 
    std::vector<double>   &p,            // vector of p-values (sorted!)
    int                   m,             // length of p
    std::vector<double>   &simesfactor,  // denominator of local test
    bool                  simes )        // assume simes yes or no
{
  // initialize output
  std::vector<double> alpha(m);
  
  // initialize
  vector<int> hull = findhull(m, p);
  double Dk=0;
  int k=hull.size()-1;
  int i=1;
  
  // algorithm for alpha*
  while (i <= m) {
    if (k > 1)
      Dk = p[hull[k-1]-1] * (hull[k] - m + i) - p[hull[k]-1] * (hull[k-1] - m +i);
    if (k > 1 && Dk < 0) {
      k--;
    } else {
      alpha[i-1] = simesfactor[i] * p[hull[k]-1] / (hull[k] - m + i);
      i++;
    }
  }
  
  // Bound alpha by 1
  if (!simes) {
    for (int i=m-1; i>=0; i--) 
      if (alpha[i] > 1)
        alpha[i] = 1;
  }
  
  // Find the cumulative maximum
  if (!simes) {
    for (int i=m-2; i>=0; i--) {
      if (alpha[i] < alpha[i+1])
        alpha[i] = alpha[i+1];
    }
  }
  
  // add (m+1)st element to alpha
  alpha.push_back(0);  
  
  return alpha;
}


// Old linearithmic procedure
// // Implementation of Procedure 1 in the paper
// void findjumps (
//     int             minrow,         // row an column range 
//     int             maxrow,         // note: index between 1 and m 
//     int             mincol,         // (not 0 and m-1)
//     int             maxcol,
//     int             m,              // length of p
//     vector<double>  &p,             // p-values (sorted!)
//     vector<double>  &simesfactor,   // denominator of simes test
//     vector<double>  &alpha )        // alpha to be filled
// {
//   // find the middle column
//   int j = (mincol + maxcol)/2;      
//   
//   // Find the minimum in the jth column
//   int start = std::max(minrow, j);    // because matrix is lower triangular
//   double Mij = 2;                     // initialize: greater than any Mij
//   int minloc = 0;                     // location of the minimum
//   double minvalue = Mij;              // value of the minimum
//   for (int i = start; i <= maxrow; i++)
//   {
//     Mij = p[i-1] / (i-j+1);
//     if (Mij < minvalue) {
//       minvalue = Mij;
//       minloc = i;
//     }
//   }
//   
//   // Fill the middle value for alpha    
//   alpha[m-j] = simesfactor[m-j+1] * minvalue;
//   
//   // Recursion to submatrices left and right
//   if (j > mincol)
//     findjumps(minrow, minloc, mincol, j-1, m, p, simesfactor, alpha);
//   if (j < maxcol)
//     findjumps(minloc, maxrow, j+1, maxcol, m, p, simesfactor, alpha);
// }
// 
// // Finds the jumps of h(alpha)
// // [[Rcpp::export]]
// std::vector<double> findalpha( 
//     std::vector<double>   &p,            // vector of p-values (sorted!)
//     int                   m,             // length of p
//     std::vector<double>   &simesfactor,  // denominator of local test
//     bool                  simes )        // assume simes yes or no
// {
//   // initialize output
//   std::vector<double> alpha(m);
// 
//   // Recursively fill alpha
//   findjumps(1, m, 1, m, m, p, simesfactor, alpha);
// 
//   // Bound alpha by 1
//   if (!simes)
//     for (int i=m-1; i>=0; i--) 
//       if (alpha[i] > 1)
//         alpha[i] = 1;
//   
//   // Find the cumulative maximum
//   if (!simes)
//     for (int i=m-2; i>=0; i--) 
//       if (alpha[i] < alpha[i+1])
//         alpha[i] = alpha[i+1];
// 
//   // add (m+1)st element to alpha
//   alpha.push_back(0);  
//       
//   return alpha;
// }

// Calculates the denominator of the local test
// [[Rcpp::export]]
std::vector<double> findsimesfactor(
    bool            simes,        // assume simes yes or no
    int             m )           // number of p-values
{
  vector<double> simesfactor(m+1);  // denominator of simes test
  double multiplier = 0;            // extra multiplier term needed for non-simes 
  simesfactor[0] = 0;
  if (simes) 
    for (int i=1; i<=m; i++)
      simesfactor[i] = i;
  else 
    for (int i=1; i<=m; i++)
    {
      multiplier += 1.0/i;
      simesfactor[i] = i * multiplier;
    }

  return simesfactor;
}

// Calculate adjusted p-values for all elementary hypotheses
// [[Rcpp::export]]
std::vector<double> adjustedElementary( 
    std::vector<double>   &p,             // vector of p-values (sorted!)
    std::vector<double>   &alpha,         // jumps of the function h
    int                   m,              // length of p
    std::vector<double>   &simesfactor )  // denominator of local test
{
  // adjust p
  vector<double> adjusted(m, 0);
  int i = 1;
  int j = m+1;
  while (i <= m) 
  {
    if (simesfactor[j-1] * p[i-1] <= alpha[j-1]) 
    {
      adjusted[i-1] = std::min(simesfactor[j] * p[i-1], alpha[j-1]);
      i++;
    }
    else
      j--;
  }
  
  return adjusted;
}
  
// Calculate the adjusted p-value of an intersection hypotheses
// [[Rcpp::export]]
double adjustedIntersection( 
    double                pI,             // p-value
    std::vector<double>   &alpha,         // jumps of the function h
    int                   m,              // length of p
    std::vector<double>   &simesfactor )  // denominator of local test
{
  // intitialize 
  int lower = 1;
  int upper = m+2;
  int mid = 0;
  
  // bisection algorithm
  // we find the largest value such that condition (*) holds
  // throughout (*) holds for lower, not for upper
  // in the loop we converge lower and upper
  while (lower < upper -1) 
  {
    mid = (lower + upper)/2;
    if (simesfactor[mid-1] * pI <= alpha[mid-1])  // (*)
    {
      lower = mid;
    }
    else
    {
      upper = mid;
    }
  }
  
  // we found lower as the largest value for which (*) holds 
  pI = std::min(simesfactor[lower] * pI, alpha[lower-1]);
  
  return pI;
} 

// Calculate the value of h(alpha) for a given alpha
// [[Rcpp::export]]
int findHalpha (std::vector<double> &jumpalpha,   // points where h jumps
                double              alpha,        // alpha where h is to be evaluated
                int                 m)            // size of the multiple testing problem
{
  // We use bisection
  int lower = 0;
  int upper = m+1;
  int mid = 0;
  while (lower+1 < upper) {
    mid = (lower + upper + 1)/2;
    if (jumpalpha[mid-1] > alpha)
    {
      lower = mid;
    }
    else
    {
      upper = mid;
    }
  }
  return lower;
}
 
// Calculates the size of the concentration set at a fixed alpha
// [[Rcpp::export]]
int findConcentration(std::vector<double> &p,           // vector of p-values (sorted!)
                      double              simesfactor,  // simesfactor at h(alpha)
                      int                 h,            // h(alpha)
                      double              alpha,        // alpha itself
                      int                 m)            // size of the problem
  
{
  // from m-h we increase z until we fulfil the condition
  int z=m-h;
  if (z > 0)  // h=m implies z=0
  {
    while ((z < m) & (simesfactor * p[z-1] > (z - m + h + 1) * alpha))
    {
      z++;
    }
  }
  return z;
}

// Find function for disjoint set data structure
// (1) Old recursive version
// int Find(int x,
//          vector<int> &parent)
// {
//   if (parent[x] != x)
//   {
//     parent[x] = Find(parent[x], parent);
//   }
//   
//   return parent[x];
// }
// (2) iterative version (more stable)
int Find(int x,
         vector<int> &parent)
{
  while (parent[x] != x)
  {
    parent[x] = parent[parent[x]];
    x         = parent[x];
  }
  
  return x;
}



// Union function for disjoint set data structure
// Extra: we keep track of the lowest entry of each disjoint set
// That way we can find the lower set to merge with
void Union(int x,
           int y,
           vector<int> &parent,
           vector<int> &lowest,
           vector<int> &rank)
{
  int xRoot = Find(x, parent);
  int yRoot = Find(y, parent);
  
  // if x and y are already in the same set (i.e., have the same root or representative)
  if (xRoot == yRoot) return; // Note: this never happens in our case
  
  // x and y are not in same set, so we merge
  if (rank[xRoot] < rank[yRoot])
  {
    parent[xRoot] = yRoot;
    lowest[yRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
  else if (rank[xRoot] > rank[yRoot])
  {
    parent[yRoot] = xRoot;
    lowest[xRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
  else
  {
    parent[yRoot] = xRoot;
    rank[xRoot]++;
    lowest[xRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
}

// Calculate the category for each p-value
int getCategory(double   p,               // p-value for which we need the category
                double   simesfactor,     // simesfactor at h(alpha)
                double   alpha,           // alpha itself
                int      m)               // size of the problem
{
  if (p==0 || simesfactor==0)
    return 1;
  else
    if (alpha == 0)
      return m+1;
    else
    {
      double cat = (simesfactor / alpha) * p;
      return static_cast<int> (std::ceil(cat));
    }
}

// Calculates the lower bound to the number of false hypotheses
// Implements the algorithm based on the disjoint set structure
// [[Rcpp::export]]
std::vector<int> findDiscoveries(std::vector<double>   &p,           // p-values in set I (not sorted)
                                 std::vector<double>   &allp,        // all p-values (sorted!)
                                 double                simesfactor,  // simesfactor at h(alpha)
                                 int                   h,            // h(alpha)
                                 double                alpha,        // alpha
                                 int                   k,            // size of I
                                 int                   m)            // size of the problem
{
  // calculate categories for the p-values
  vector<int> cats;
  for (int i=0; i<k; i++)
  {
    cats.push_back(getCategory(p[i], simesfactor, alpha, m));
  }
  
  // find the maximum category needed
  int z = findConcentration(allp, simesfactor, h, alpha, m);
  int maxcat = std::min(z-m+h+1, k);
  int maxcatI = 0;
  for (int i=k-1; i >= 0; i--)
  {
    if (cats[i] > maxcatI) {
      maxcatI = cats[i];
      if (maxcatI >= maxcat) break; 
    }
  }
  maxcat = std::min(maxcat, maxcatI);
  
  // prepare disjoint set data structure
  vector<int> parent;
  vector<int> lowest;
  vector<int> rank;
  for (int i=0; i <= maxcat; i++)
  {
    parent.push_back(i);
    lowest.push_back(i);
    rank.push_back(0);
  }

  // The algorithm proper. See pseudocode in paper
  vector<int> discoveries(1,0);
  int lowestInPi;
  for (int i=0; i < k; i++)
  {
    if (cats[i] <= maxcat)
    {
      lowestInPi = lowest[Find(cats[i], parent)];
      if (lowestInPi == 1)
      {
        discoveries.push_back(discoveries.back()+1);
      }
      else
      {
        discoveries.push_back(discoveries.back());
        Union(lowestInPi-1, Find(cats[i], parent), parent, lowest, rank);
      }
    }
    else
    {
      discoveries.push_back(discoveries.back());
    }
  }
  
  return discoveries;
}

