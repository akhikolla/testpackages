/***************
* class TriArray
* implements a triangular array as a template
* first index may never exceed second index!
* second index running fastest!
* Thomas Hotz, 2007-2015
***************/

template <class T>
class TriArray {
  public:
    unsigned int size;
    T *array;
    
    TriArray(unsigned int n); // constructor creating the array, doesn't initialize!
    TriArray(unsigned int n, T* arr); // constructor accepting an array, assuming it already is of the right size!
    TriArray(); // default constructor may not be called
    virtual ~TriArray() {}; // default destructor, doesn't free the array as this is done by R
    
    virtual T& operator()(unsigned int i, unsigned int j); // index element
    virtual const T& operator()(unsigned int i, unsigned int j) const; // index element of const object
    
};

template <class T>
TriArray<T>::TriArray() {
  error("TriArray needs a size!");
}

template <class T>
TriArray<T>::TriArray(unsigned int n) {
  if(n <= 0) error("TriArray needs a postive size!");
  size = n;
  array = (T*) R_alloc(( n * ( n + 1 )) / 2, sizeof(T)); // allocate storage, freed by R
}

template <class T>
TriArray<T>::TriArray(unsigned int n, T* arr) {
  if(n <= 0) error("TriArray needs a postive size!");
  size = n;
  array = arr;
}

template <class T>
T& TriArray<T>::operator()(unsigned int i, unsigned int j) {
  if(i >= size) error("First index out of bound!");
  if(j >= size) error("Second index out of bound!");
  if(i > j) error("First index may not exceed second index!");
  return array[i * size - (i * ( i - 1 )) / 2 + ( j - i )]; // we skip a triangle of size i-1 before this element
}

template <class T>
const T& TriArray<T>::operator()(unsigned int i, unsigned int j) const {
  if(i >= size) error("First index out of bound!");
  if(j >= size) error("Second index out of bound!");
  if(i > j) error("First index may not exceed second index!");
  return array[size * i - (i * ( i - 1 )) / 2 + ( j - i )]; // we skip a triangle of size i-1 before this element
}


/***************
* class TriArrayFF
* sub-class of TriArray, but:
* first index running fastest!
* Thomas Hotz, 2007
***************/

template <class T>
class TriArrayFF : public TriArray<T> {
  public:
    TriArrayFF(unsigned int n) : TriArray<T>(n) {}; // constructor creating the array, doesn't initialize!
    TriArrayFF(unsigned int n, T* arr) : TriArray<T>(n, arr) {}; // constructor accepting an array, assuming it already is of the right size!
    TriArrayFF(); // default constructor may not be called
    virtual ~TriArrayFF() {}; // default destructor, doesn't free the array as this is done by R
    
    virtual T& operator()(unsigned int i, unsigned int j); // index element
    virtual const T& operator()(unsigned int i, unsigned int j) const; // index element of const object
};

template <class T>
T& TriArrayFF<T>::operator()(unsigned int i, unsigned int j) {
  if(i >= this->size) error("First index out of bound!");
  if(j >= this->size) error("Second index out of bound!");
  if(i > j) error("First index may not exceed second index!");
  return this->array[(j * ( j + 1 )) / 2 + i]; // there is a triangle of size j before this element
}

template <class T>
const T& TriArrayFF<T>::operator()(unsigned int i, unsigned int j) const {
  if(i >= this->size) error("First index out of bound!");
  if(j >= this->size) error("Second index out of bound!");
  if(i > j) error("First index may not exceed second index!");
  return this->array[(j * ( j + 1 )) / 2 + i]; // there is a triangle of size j before this element
}
