/* $Id: oja_geometry.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef OJA_GEOMETRY_H
#define OJA_GEOMETRY_H

#include "error.h"
#include "data.h"
#include "hyperplane.h"
#include "line.h"
#include "misc.h"

using namespace std; //df

class OjaData;
class OjaLine;
class Lattice;

enum OjaMedianMethod
{
	EVAL_ALL_POINTS, GRADIENT_DESCENT, BRUTE_FORCE, FOLLOW_INTERSECTION_LINES, FOLLOW_INTERSECTION_LINES_BOUNDED, FOLLOW_INTERSECTION_LINES_BOUNDED_APPROX,
	LATTICE_APPROX,LATTICE_APPROX2,LATTICE_APPROX3,BOOTSTRAP,SIMPLEX_APPROX
};

enum LatticeMeasure
{
	LM_DIAMETER,LM_EDGE_LENGTH_AVG,LM_EDGE_LENGTH_MAX
};

class OjaPoint:private Point 
{    
    IndexSet id;
    const OjaData* data;
	
  public:
    
    OjaPoint() : Point() {data = 0;}
    OjaPoint(const OjaData& d) : Point()
		{data = &d;}
    OjaPoint(const OjaData* d) : Point()
		{data = d;}
    OjaPoint(const OjaPoint& p) : Point(p)
		{id = p.id; data = p.data;}
    OjaPoint& operator=(const OjaPoint& p)
		{Point::operator=(p); id = p.id; data = p.data; return *this;}

    void get(const IndexSet& I);
    void set_index(const IndexSet&);
    void set_index(const IndexSet&,const Index&);
    void set_location(const Point& p)
		{Point::operator=(p);}
    
    IndexSet index() const
		{return id;}
    Point location() const
		{return (Point)*this;}
	bool is_nil() const
		{return Point::is_nil();}
	const OjaData* data_set() const
		{return data;}

	OjaLine scan_all_routes(OjaPoint& best,double& ojafn,double hi_score=-1.0);
	bool better_route(OjaLine& L,OjaPoint& P,double &ojafn);
	
    // BUG: n�m� pit�isi m��ritell� attribuuteiksi luokan sis�lle ja
    //      palauttaa valmiiksi laskettu arvo mik�li se on ajan tasalla.
    bool is_data() const
		{return id.how_many_common_digits()==1;}
    int data_index() const; // BUG: yhten�ist� seuraavan kanssa
	bool splits_line(int& p1,int& p2,Index& H) const;
	bool splits_line() const
		{int dummy1,dummy2; Index I; return splits_line(dummy1,dummy2,I);}


	// BUG: n�m� toimii periaatteessa v��rin (erottelee saman pisteen eri
	// esitykset)
    bool operator==(const OjaPoint& p2) const
		{return id==p2.id;}
    bool operator!=(const OjaPoint& p2) const
		{return !(operator==(p2));}
};

class OjaData: public Data
{
    OjaMedianMethod method;
	LatticeMeasure lattice_measure;
    HyperplaneSet* plane;
    int planes;
	set<int> includedPlanes;	// Planes included in bounded exact search algorithm
    Index* planeindex;
    Point norm_offset,norm_scale;
	Point exact_median; // Tarkka mediaani, jos tiedossa (tiedostosta)
	Point boundedMin;	Point boundedMax;  // min || max defined by bounds
	int original_size;

	// Parametereja eri menetelmille
	double epsilon;
	double chi2_limit;
	int set_size;
	int max_searchlines;
	double volume;

    double confidence_size(Lattice* L);
	
    OjaPoint medianAtDataPoints(double& value_at_min);
    OjaPoint medianEvalAllPoints();
    OjaPoint medianGradientDescent();
    OjaPoint medianFollowIntersectionLines();
	OjaPoint medianFollowIntersectionLinesBounded();
	OjaPoint medianBruteForceSearch();
	OjaPoint medianLatticeApprox();
	OjaPoint medianLatticeApprox2();
	OjaPoint medianLatticeApprox3(list<Hyperplane>* store=0,list<Index>* idxstore=0,matrix* Cov=0);
	OjaPoint medianBootstrap(const list<Hyperplane>& store,const list<Index>& idxstore,const matrix& CovInv);
	OjaPoint medianSimplexApprox();
	void brute_force_search(OjaLine& L,OjaPoint& p,double hi_score);
	
  public:
#ifdef _MSC_VER 
#ifdef DEEPDEBUG
	static OjaData S;
#endif
#endif
	Point h; double h0;

    OjaData();
    OjaData(const char* filename);  
    OjaData(int vecdim,int size);
    ~OjaData();
    
    void set_median_method(OjaMedianMethod new_method)
		{method=new_method;}
	void set_lattice_measure(LatticeMeasure m)
		{lattice_measure=m;}
	void set_epsilon(double eps)
		{epsilon=eps;}
	void set_chi2_limit(double val)
		{chi2_limit=val;}
	void set_set_size(int sz)
		{set_size=sz;}
	void set_max_searchlines(int ml)
		{max_searchlines=ml;}
	void set_volume(double vol)
		{volume = vol;}

    OjaPoint median();
	Point exact() const
		{return exact_median;}
	OjaPoint medianBootstrap(matrix& covariance,int how_many);
    double oja(const Point& x) const;
    double oja(const OjaPoint& x) const
		{return oja(x.location());}
    Point gradient(const Point& x) const;
	Point oja_rank(const Point& x) const;
	Hyperplane getBoundingHyperplane(Point& point, Point& grad);
	Hyperplane getBoundingHyperplane(Point& point);
    void get_oja_and_gradient(const Point& x,double& oja,Point& grad) const;
    bool derivative(const Point& x,const Point& direction,double& Dpos,double& Dneg) const; 
    bool derivative(const Point& x,const Line& L,double& Dpos,double& Dneg) const 
		{return derivative(x,L.dir(),Dpos,Dneg);}
	
    void generate_hyperplanes();
    void regenerate_hyperplanes();
    const Hyperplane& hyperplane(int i) const;
    Index hyperplaneindex(int i) const; 
    const HyperplaneSet& hyperplaneset() const; 
    int points() const
		{return size();}
    int hyperplanes() const
		{return planes ? planes : choices(size(),dim());}
	set<int> get_includedPlanes() const
		{ return includedPlanes; }
    int crossingpoints() const
		{return choices(hyperplanes(),dim());}
	
    void scale();
    Point scaled(Point x) const;
    Point descaled(Point x) const;

	Point min() const;
	Point max() const;
	void set_bounded_min_max(const Point bmin, const Point bmax);
	void add_bound_points(const vector<Point> & crossing_points);
	void add_bound(const Hyperplane& bound, const set<int>& crossings);
	int get_original_size() const 
		{return original_size;}
};

class OjaLine:public Line 
{
    const OjaData *data;
    IndexSet id;
	    
    OjaLine& operator=(const Line&);
	
  public:
    
    // BUG: tarvitaan kunnollinen is_nil k�sittely koko luokalle
    OjaLine();
    ~OjaLine() {}
    
    OjaLine(const OjaData& D);
    OjaLine(const OjaData* D);
    OjaLine(const OjaLine&);
    OjaLine& operator=(const OjaLine&);
    void set(const IndexSet& I,const Line& L);

	int dim() const;
	inline const OjaData* Data() const
		{return data;}
    Line line() const
		{return (Line)*this;}
 	bool is_nil() const
 		{return data==0 || id.dim()==0;} 
	    	
    IndexSet index() const {return id;}
    
    void get(const IndexSet& I);
    void get_random_through(int index);
    void get_random_through(int index1,int index2);
	OjaPoint crossing_point(const Index& hyperplane) const;
	
    OjaPoint min(double& ojafn) const; /* New */
    OjaPoint min() const
		{double dummy; return min(dummy);}
	
    bool operator==(const OjaLine& L2)
		{return id==L2.id;}
    bool operator!=(const OjaLine& L2)
		{return !(id==L2.id);}
    bool is_same(const OjaLine& L2);
};

class DotSet
{
	const OjaLine* line;
	
	bool sorted;
	list<pair<double,int> > dotlist; // K�ytet��n, jos sorted=false
	vector<pair<double,int> > *dotarray; // K�ytet��n, jos sorted=true
	
	Point h; // Ulkopuolelta kulkevien hypertasojen kertoimet
	double h0;

	void generate_dots();
	void generate_dots_bounded();
	void get_common_coefs(Point& h, double& h0);

	set<int> find_valid_bounds(set<int>& includedPlanes, const OjaData* data, Point& x);
	
public:

	DotSet(const OjaLine& L);
	DotSet(const OjaLine* L);
	DotSet(const OjaLine* L, Point& h, double& h0);
	~DotSet();
	
	void sort();
	
	int size() const
		{return sorted ? dotarray->size() : dotlist.size();}
	int dim() const
		{return line->dim();}
	pair<double,int>& dot(int index) const;
	
	inline const OjaData* Data() const
		{return line->Data();}
	Point point_at(double t) const
		{return line->line().at(t);}
	Point point(const pair<double,int>& tt) const
		{return point_at(tt.first);}
	Point point(int idx) const
		{return point(dot(idx));}
	const Hyperplane& hyperplane(const pair<double,int>& tt) const
		{return Data()->hyperplane(tt.second);}
	const Hyperplane& hyperplane(int idx) const
		{return hyperplane(dot(idx));}

	OjaPoint min_evaluate_all(double &ojafn);
	OjaPoint min(double &ojafn); /* New */
	OjaPoint min()
		{double dummy; return min(dummy);}
	
	double oja(double at) const;
	double oja(const pair<double,int>& tt)
		{return oja(tt.first);}
};

class OjaLineSet 
{
    list<OjaLine> line;
    list<OjaLine>::iterator last_access;
    int last_access_index;
    const OjaData* data;
    
    OjaLineSet()
		{last_access_index=-1; data=0;}
    
  public:
    
    OjaLineSet(const OjaData& D)
		{data=&D; last_access_index=-1;}
    OjaLineSet(const OjaData* D)
		{data=D; last_access_index=-1;}
    ~OjaLineSet()
		{}
    
    void clear()
		{line.clear(); last_access_index=-1;}
    void add(const OjaLine& L)
		{line.push_back(L);}
    
    int size()const
		{return line.size();}
    OjaLine& operator[](int i);
    OjaLine best_at(const OjaPoint& z);
    
	void make_all_combinations(int with);
    void make_data_combinations(int with);
    void make_combinations(const IndexSet& I);
    void make_split_line_combinations(int p1,int p2,const Index& I);
	void make_combinations(const OjaPoint& p);
	void make_full_combinations(const OjaPoint& p);
	
	int scan_all_routes(OjaPoint& best,double& ojafn);
};

class OjaLineIndexIterator
{
	enum {IT_SPLIT,IT_POINTOPOINT,IT_DATA,IT_NORMAL} mode;

	int i,j;
	Index I;
	IndexSet S;

	bool overflow;

	void refresh_value();

	IndexSet current;
	
  public:
	
	OjaLineIndexIterator(const OjaPoint& p,
	  bool cheat_with_data=false,bool cheat_whit_split=false); 

	const IndexSet& index() const
		{return current;}
	
	OjaLineIndexIterator& operator++(int);
	operator void*() const
		{return overflow ? (void*) 0 : (void *) -1;}
};

class OjaLineIterator: public OjaLineIndexIterator
{
	OjaLine L;
	
  public:

	OjaLineIterator(const OjaPoint& p);
	OjaLineIterator(const OjaPoint& p,bool cheat_with_data
	  ,bool cheat_whit_split);

	const OjaLine& line();
};

ostream& operator <<(ostream& os,const OjaPoint& L);
ostream& operator <<(ostream& os,const OjaLine& L);
ostream& operator <<(ostream& os,const OjaLineSet& S);

// T�m� tarvitaan, jotta leda_list rakennetta voi k�ytt�� OjaLine luokalle
istream& operator >>(istream& is,OjaLine& L);

#endif
