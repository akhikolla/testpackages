#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <algorithm> 
#include <deque>

#include "constants.h"
#include "templates.h"

//Namespaces
using namespace Rcpp;
using namespace std;

#ifndef __structs__
#define __structs__

//Data Structures

//2D vector
struct Vec2 {
  Vec2() : x(0), y(0){}
  Vec2(double x, double y) : x(x), y(y){}
  double x, y;
  double gradient() const { return x == 0 ? D_MAX : y / x; }
  double length() const {   return pow(pow(x,2.0) + pow(y,2.0),0.5); }
  bool operator  ==(const Vec2& b) const{ return x == b.x && y == b.y; }
  bool operator  !=(const Vec2& b) const{ return !(*this == b); }
  Vec2& operator  =(const Vec2& b) { x = b.x; y = b.y; return *this;}
  Vec2  operator  +(const Vec2& b) const { return Vec2(x + b.x, y + b.y) ;}
  Vec2  operator  -(const Vec2& b) const { return Vec2(x - b.x, y - b.y) ;}
  Vec2  operator  *(const double& t) const { return Vec2(x*t,y*t); }
  Vec2  operator  /(const double& t) const { return Vec2(x/t,y/t); }
  Vec2& operator +=(const Vec2& b) {  x += b.x; y += b.y; return *this; }
  Vec2& operator -=(const Vec2& b) {  x -= b.x; y -= b.y; return *this; }
  Vec2& operator *=(const double& t){ x *= t; y *= t; return *this; }
  Vec2& operator /=(const double& t){ x /= t; y /= t; return *this; }
  friend std::ostream &operator<<(std::ostream& os, Vec2 const &b){  os << "[" << b.x << "," << b.y << "]";  return os;  }
  bool inLine(Vec2 p1, Vec2 p2, bool strictlyBetween = true){
    double xmn = min(p1.x,p2.x), xmx = max(p1.x,p2.x), ymn = min(p1.y,p2.y), ymx = max(p1.y,p2.y);
    bool isBetween = (!strictlyBetween) || (x >= xmn - D_TOL && x <= xmx + D_TOL && y >= ymn - D_TOL && y <= ymx + D_TOL);
    if(!isBetween){ return false; }
    if(abs(p1.x-p2.x) < D_TOL) return abs(p1.x-x) < D_TOL; //Vertical   Case
    if(abs(p1.y-p2.y) < D_TOL) return abs(p1.y-y) < D_TOL; //Horizontal Case
    double slope = (p2.y-p1.y)/(p2.x-p1.x);
    return (abs(slope*(x - p1.x) - y) < D_TOL);
  }
};

//2D Vector with ID Field
struct Node2: Vec2{
  Node2() : Vec2(), id(0){}
  Node2(int id, double x, double y) : Vec2(x,y),id(id){}
  Node2(int id, Vec2 vec) : Vec2(vec), id(id) {}
  int id;
  bool operator ==(const Node2& n) const { return Vec2::operator==(n) && id == n.id; }
  bool operator !=(const Node2& n) const { return !(*this == n); }
};

//3D vector, inherits fro 2D vector
struct Vec3 : Vec2 {
  Vec3(): Vec2(), z(0) {}
  Vec3(Vec2 b, double z) : Vec2(b),z(z){}
  Vec3(double x, double y, double z) : Vec2(x,y), z(z) {}
  double z;
  double gradient() const { return (x + y == 0) ? D_MAX : z / Vec2::length(); }
  double length() const {return pow(pow(x,2) + pow(y,2) + pow(x,2),0.5); }
  Vec2 xy(){return Vec2(x,y);} //Point representing the x-y projection
  Vec2 xz(){return Vec2(x,z);} //Point representing the x-z projection
  Vec2 yz(){return Vec2(y,z);} //Point representing the y-z projection
  bool operator  ==(const Vec3& b) const { return Vec2::operator==(b) && z == b.z; }
  bool operator  !=(const Vec3& b) const { return !(*this == b); }
  Vec3& operator  =(const Vec3& b) { Vec2::operator=(b); z = b.z; return *this;}
  Vec3  operator  +(const Vec3& b) const { return Vec3(x + b.x, y + b.y, z + b.z) ;}
  Vec3  operator  -(const Vec3& b) const { return Vec3(x - b.x, y - b.y, z - b.z) ;}
  Vec3  operator  *(const double& t) const {return Vec3(x*t,y*t,z*t); }
  Vec3  operator  /(const double& t) const {return Vec3(x/t,y/t,z/t); }
  Vec3& operator +=(const Vec3& b) {  Vec2::operator+=(b); z += b.z; return *this; }
  Vec3& operator -=(const Vec3& b) {  Vec2::operator-=(b); z -= b.z; return *this; }
  Vec3& operator *=(const double& t){ Vec2::operator*=(t); z *= t;   return *this; }
  Vec3& operator /=(const double& t){ Vec2::operator/=(t); z /= t;   return *this; }
  friend std::ostream &operator<<(std::ostream& os, Vec3 const &b){  os << "[" << b.x << "," << b.y << "," << b.z << "]";  return os; }
  bool inLine(Vec3 p1, Vec3 p2, bool strictlyBetween){
    return (*this).xy().inLine(p2.xy(),p2.xy(),strictlyBetween) &&
           (*this).xz().inLine(p2.xz(),p2.xz(),strictlyBetween) &&
           (*this).yz().inLine(p2.yz(),p2.yz(),strictlyBetween);
  }
};

//Data structure to store Node information
struct Node3 : Vec3{
  Node3(){};
  Node3(int id,Vec3 vec) : Vec3(vec),id(id){}
  Node3(int id, double x, double y, double z) : Vec3(x,y,z),id(id){}
  int id;
  bool operator ==(const Node3& rhs) const { return id == rhs.id; }
  bool operator !=(const Node3& rhs) const { return id != rhs.id; }
  bool operator < (const Node3& rhs) const { return id <  rhs.id; }
  bool operator > (const Node3& rhs) const { return id >  rhs.id; }
};

//Data Structure to enable sorting of deleaunay triangles.
struct Centroid{
  //Centroid(); //Intentionally declared but not defined
  Centroid() : index(0), x(0), y(0) { }
  Centroid(int index) : index(index), x(0), y(0){ }
  //Centroid(int index, double x, double y) : index(index), x(x), y(y) {}
  int index;
  double x;
  double y;
  bool operator<(const Centroid& rhs) const { 
    //return atan2(y,x) < atan2(rhs.y,rhs.x); 
    return (x < rhs.x) ? true : ((x == rhs.x) ? y < rhs.y : false);
  }
  Centroid& operator +=(const Vec2& b) {  x += b.x; y += b.y; return *this; }
};

//Data structure to store Edge information:
struct Edge{
  Edge(){}; //Intentionally declared but not defined
  Edge(Node3 A, Node3 B){
    if(A == B)
      throw std::invalid_argument("A and B must have unique indexes");
    if(A < B){a = A; b = B;}else{a = B; b = A;}
  }
  Node3 a, b;
  bool isBound(const double& z) const { return (z > a.z && z < b.z) || (z < a.z && z > b.z) ; } 
  bool isAbove(const double& z) const { return a.z > z && b.z > z; }
  bool isBelow(const double& z) const { return a.z < z && b.z < z; }
  Node3 highest(){ return a.z > b.z ? a : b; }
  Node3 lowest(){  return a.z < b.z ? a : b; }
  //Determines is an edge, e, is contained within this object in terms of its z level
  bool isBound(const Edge& e){ return isBound(e.a.z) && isBound(e.b.z); }
  double length(){return (a-b).length(); }
  bool operator ==(const Edge& rhs) const {
    //In the following equality check, the second operation *should* not be required, 
    // since a and b are put in order when Edge is constructed
    return (a == rhs.a && b == rhs.b) || (a == rhs.b && b == rhs.a);
  }
};

struct ContourData{ 
  ContourData(int levelID, int groupID, int pathID, double x, double y, double z) : 
    levelID(levelID), groupID(groupID), pathID(pathID), x(x), y(y), z(z){}
  NumericVector toNumericVector() const {
    NumericVector ret(6);
    ret(0)=levelID; ret(1)=groupID; ret(2)=pathID; ret(3)=x; ret(4)=y; ret(5)=z;
    return ret;
  }
  private:
      int levelID,
          groupID,
          pathID;
      double x, y, z;
};

//Forward Declaration.
Vec3 interpolateZVal(Vec3 a, Vec3 b, double z);

//Data structure for storing a deleaunay triangle, constructed by three (3) unique Nodes,
//internally it will store three (3) Edges...
struct Del{
  Del(); //Intentionally declared but not defined
  //Constructor, being made up of three nodes
  Del(Node3 &A, Node3 &B, Node3 &C) {
    if(A == B || A == C || B == C){
      throw std::invalid_argument("The Nodes A, B and C must be unique");
    }
    //Put the edges in order
    vector<Node3> n = vector<Node3>(); n.push_back(A); n.push_back(B); n.push_back(C);
    std::sort(n.begin(),n.end());
    for(int i = 0; i < 3; i++){
      sisters[i] = NULL;
      edges[i] = Edge(n[i],n[ (i<2?i+1:0) ]);
    }
    aspectRatio = calcAspectRatio();
  }
  public:
    //Public Members
    Edge edges[3]; 
    Del* sisters[3];
    vector<Del*> cousins;
    
    //reset the drawn levels
    void reset(){ 
      levelsDrawn.clear(); 
    }
    
    //Function to recursively assign pointers in two adjacent and joining Dels.
    bool makeSisters(Del* sister){
      if(sister == NULL) 
        return false;
      if(this->isFull())
        return sister == this->sisters[0] || sister == this->sisters[1] || sister == this->sisters[2]; 
      for(int i = 0; i < 3; i++){
        if(this->sisters[i] == sister)
          return true;
        for(int j = 0; j < 3; j++){
          if(sister->edges[j] == edges[i]){
            this->sisters[i] = sister; 
            return sister->makeSisters(this); //Network the counterpart
          } 
        }
      }; return false;
    }
    
    bool makeCousins(Del* cousin){
      if(cousin == NULL){return false;}
      if(isIn<Del*>(cousins,cousin)){return true;}
      for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
          if(cousin->edges[j].a == edges[i].a){
            cousins.push_back(cousin);
            return cousin->makeCousins(this);
          }
        }
      }; return false;
    }
    
    bool isBound(const double& z) const {
      double zarr[4] = {edges[0].a.z,edges[0].b.z,
                        edges[1].a.z,edges[1].b.z};
      double  zMin = zarr[0],
              zMax = zarr[0];
      for(int i = 1; i < 4; i++){
        zMin = min(zMin,zarr[i]);
        zMax = max(zMax,zarr[i]);
      }
      if(zMin < z && zMax > z) return true;
      return false;
    }      
    
    //Has this Del been fully networked, ie not on the hull
    bool isFull() const { return numberOfSisters() == 3; }
    
    //Does the del touch the hull
    bool touchesHull(const vector<Node3>& chull) const{
      for(int i = 0; i < 2; i++){
        if(isIn<Node3>(chull,edges[i].a) || isIn<Node3>(chull,edges[i].b))
          return true;
      }
      return false;
    }
    
    //Draw Contour Recursively
    void drawContour( const double& level, 
                      const vector<Node3>& chull, 
                      vector<Vec3>& contours, 
                      Del* from, 
                      vector<Del>& delLibrary,
                      const double criticalRatio){  
      
      //Check candidature and availability
      if(isIn(levelsDrawn,level) || !isBound(level) || (aspectRatio >= abs(criticalRatio) && touchesHull(chull)))
        return;
      
      for(unsigned int i = 0, cnt = 0; cnt < 2 && i < 3; i++){
        Del* sister  = sisters[i];
        Edge* edge   = &edges[i]; 
        
        //Must be a different Neighbour
        if((from == NULL || sister != from)  && edge->isBound(level)) {
        
          //Reverse
          if(cnt > 0 && from == NULL)
            std::reverse(contours.begin(),contours.end());
            
          //Make this level unavailable anymore...
          levelsDrawn.push_back(double(level));
          
          //Push back the result for this level.
          contours.push_back(interpolateZVal(edge->a,edge->b,level));
          
          //If Share an Edge, Execute Recursively, Referencing the 'from' pointer as 'this'
          if(sister != NULL)
            sister->drawContour(level,chull,contours,this,delLibrary,criticalRatio);
          //Increment the Counter
          cnt++;
        }
      }
    }
  
  //Private Members and functions
  private:
    vector<double> levelsDrawn;
    double aspectRatio;
    double calcAspectRatio(){
      double   a = edges[0].length(), 
               b = edges[1].length(), 
               c = edges[2].length();
      double   s = (a + b + c) / 2.0;
      return a*b*c/(8*(s-a)*(s-b)*(s-c));
    }
    int  numberOfSisters() const {
      int n = 0;
      for(int i = 0; i < 3; i++)
        n += this->sisters[i] != NULL ? 1 : 0;
      return n;
    }
};

#endif

