#ifndef BINTREE__H
#define BINTREE__H

#include <stack>

#include <Rversion.h>
#if defined(R_VERSION) && R_VERSION >= R_Version(3, 3, 0)
  /* nothing */
#else
  #define NO_C_HEADERS true // disables the including of stdlib.h, stdio.h, limits.h, math.h by R.h
  #include <cstdlib>        // manually loading of cpp versions of disabled headers
  #include <cstdio>
  #include <climits>
  #include <cmath>
  #include <cstddef>
  using std::size_t;
#endif

#include <R.h>
#include <Rinternals.h>

/***************
* class BinTree
* implements a binary tree as a template
* including functions to move through the tree
* Thomas Hotz, 2007-2008
***************/

template <class T>
class BinTree {
  private:
    class Node {
      public:
        T value; // the value
        Node *left; // the left node
        Node *right; // the right node
        bool isRight; // is this a right node? true for root
//         Node(T v, bool r); // constructor
        Node();
        Node(T v, bool r) : value(v), left(NULL), right(NULL), isRight(r) {}
        Node(Node* l, Node* r, bool i) : left(l), right(r), isRight(i) {}
    };

    mutable std::stack<Node*> s; // the stack records where we are in the tree
    Node *root; // the root
    int numNodes; // number of nodes

  public:
//     BinTree(const T& value); // constructor turning the root into a leaf
    BinTree(T value); // constructor turning the root into a leaf

    bool isLeaf(); // is the current node a leaf?
//     const T& getValue(); // get the current node's value
    T getValue(); // get the current node's value
//     void setValue(const T& value); // set the current node's value
    void setValue(T value); // set the current node's value
    int size(); // get the number of nodes

    int depth(); // get the depth of the current node

    void up(); // go up
    void left(); // go down and left
    void right(); // go down and right
    bool previous(); // go the previous leaf to the left or return false if there is none
    bool next(); // go to the next leaf to the right or return false if there is none
    void first(); // go to the first leaf
    void last(); // go to the last leaf

//     void addLeft(const T& left); // add left leaf
    void addLeft(T value); // add left leaf
//     void addRight(const T& right); // add right leaf
    void addRight(T value); // add right leaf

    BinTree(); // default constructor may not be called
};

template <class T>
BinTree<T>::BinTree() {
  error("BinTree needs to contain at least one node!");
}

template <class T>
// BinTree<T>::BinTree(const T& value) {
BinTree<T>::BinTree(T value) {
  Node *root = (Node*) R_alloc(1, sizeof(Node));
  *root = Node(value, TRUE); // new root is right leaf
  s.push(root);
  numNodes = 1;
}

// template <class T>
// BinTree<T>::Node::Node(T v, bool r) : value(v), left(NULL), right(NULL), isRight(r) {}

template <class T>
bool BinTree<T>::isLeaf() {
  return s.top()->left == NULL && s.top()->right == NULL;
}

template <class T>
// const T& BinTree<T>::getValue() {
T BinTree<T>::getValue() {
  return s.top()->value;
}

template <class T>
// void BinTree<T>::setValue(const T& value) {
void BinTree<T>::setValue(T value) {
  s.top()->value = value;
}

template <class T>
int BinTree<T>::size() {
  return numNodes;
}

template <class T>
int BinTree<T>::depth() {
  return s.size();
}

template <class T>
void BinTree<T>::up() {
  if(s.size() == 1) {
    error("There is no element above the root!");
  } else {
    s.pop();
  }
}

template <class T>
void BinTree<T>::left() {
  if(isLeaf()) {
    error("There is no element below a leaf!");
  } else {
    s.push(s.top()->left);
  }
}

template <class T>
void BinTree<T>::right() {
  if(isLeaf()) {
    error("There is no element below a leaf!");
  } else {
    s.push(s.top()->right);
  }
}

template <class T>
bool BinTree<T>::previous() {
  if(isLeaf()) {
    // go as far up as needed
    while(s.size() > 1 && !s.top()->isRight) {
      s.pop();
    }
  }

  if(s.size() > 1) {
    // turn left
    s.pop();
    left();
    // go as far right as possible
    while(!isLeaf()) {
      right();
    }
    return TRUE;
  } else {
    // we were at the very left
    first();
    return FALSE;
  }
}

template <class T>
bool BinTree<T>::next() {
  if(isLeaf()) {
    // go as far up as needed
    while(s.size() > 1 && s.top()->isRight) {
      s.pop();
    }
  }

  if(s.size() > 1) {
    // turn right
    s.pop();
    right();
    // go as far left as possible
    while(!isLeaf()) {
      left();
    }
    return TRUE;
  } else {
    // we were at the very right
    last();
    return FALSE;
  }
}

template <class T>
void BinTree<T>::first() {
  // go to root
  while(s.size() > 1) {
    s.pop();
  }
  // go to the very left
  while(!isLeaf()) {
    left();
  }
}

template <class T>
void BinTree<T>::last() {
  // go to root
  while(s.size() > 1) {
    s.pop();
  }
  // go to the very right
  while(!isLeaf()) {
    right();
  }
}

template <class T>
// void BinTree<T>::addLeft(const T& value) {
void BinTree<T>::addLeft(T value) {
  if(isLeaf()) {
    Node *r = s.top(); // current will become right node
    s.pop();

    Node *l = (Node*) R_alloc(1, sizeof(Node));
    *l = Node(value, FALSE); // new left node

    Node *m = (Node*) R_alloc(1, sizeof(Node));
    *m = Node(l, r, r->isRight); // new middle node
    r->isRight = TRUE;

    if(s.size() > 0) { // if we aren't at the top
      Node* top = s.top(); // the node to hold the new middle node
      if(m->isRight) {
        top->right = m;
      } else {
        top->left = m;
      }
    } else {
      root = m;
    }

    s.push(m);
    numNodes += 1;
  } else {
    error("Cannot add element to non-leaf!");
  }
}

template <class T>
// void BinTree<T>::addRight(const T& value) {
void BinTree<T>::addRight(T value) {
  if(isLeaf()) {
    Node *l = s.top(); // current will become left node
    s.pop();

    Node *r = (Node*) R_alloc(1, sizeof(Node));
    *r = Node(value, TRUE); // new right node

    Node *m = (Node*) R_alloc(1, sizeof(Node));
    *m = Node(l, r, l->isRight); // new middle node
    l->isRight = FALSE;

    if(s.size() > 0) { // if we aren't at the top
      Node *top = s.top(); // the node to hold the new middle node
      if(m->isRight) {
        top->right = m;
      } else {
        top->left = m;
      }
    } else {
      root = m;
    }

    s.push(m);
    numNodes += 1;
  } else {
    error("Cannot add element to non-leaf!");
  }
}

#endif
