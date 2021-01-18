/*
*  bounded_search.cpp : Algorithms for computing exact Oja Median by bounded search.
*  Copyright (C) 2015 Oleksii Pokotylo
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include <math.h>
#include <set>
#include <algorithm>
#include "lattice.h"
#include "random.h"
#include "data.h"
#include "simplex.h"
#include "line.h"
#include "oja_geometry.h"
#include "global.h"
#include "matrix_wrapper.h"
#include "stl_tools.h"
#include <R.h>

#include "bitset"

using namespace std;

#include <ctime>
#include <fstream>
#include <sstream>

double inline dist(Point p1, Point p2)
{
	errif(p1.dim() != p2.dim(), "distance: points have different dimensions");
	double sum = 0;
	for (int i = 0; i < p1.dim(); i++){
		sum += pow(p1.coord(i) - p2.coord(i), 2);
	}
	return sqrt(sum);
}

double getVolume(Point& min, Point& max)
{
	int dim = min.dim();
	double volume = 1;
	{
		Point diff = max - min;
		for (int i = 0; i < dim; i++){
			volume *= diff.coord(i);
		}
	}
	return volume;
}

void bounded_min_max(vector<Point>& crossing_points, set<int>& used_crossing_points, Point& bmin, Point& bmax, Point& bmid, int dim){
	bmin = crossing_points[*(used_crossing_points.begin())]; bmax = bmin; bmid = 0*bmin;
	for (set<int>::iterator i = used_crossing_points.begin(); i != used_crossing_points.end(); i++)
	for (int j = 0; j < dim; j++){
		bmid[j] += crossing_points[*i][j];
		if (bmin[j] > crossing_points[*i][j]) bmin[j] = crossing_points[*i][j];
		if (bmax[j] < crossing_points[*i][j]) bmax[j] = crossing_points[*i][j];
	}
	bmid *= 1.0/used_crossing_points.size();
}

void addBound(Hyperplane& plane, bool initial, vector<Hyperplane>& bounds, vector<Point>& crossing_points, set<int>& used_crossing_points, vector<set<int> >& bounds_crossing_indexes, unsigned int dim){
	HyperplaneSet hs(dim);
	vector<Hyperplane> hv(dim);

	// check if it is crossing the region
	if (!initial){
		int side = 0;
		for (set<int>::iterator i = used_crossing_points.begin(); i != used_crossing_points.end(); i++){
			int s = plane.side(crossing_points[*i]);
			if (s == 0) continue;
			if (side == 0) {
				side = s; continue;
			}
			if (s != side){ // at least 2 points lie strictly on different sides
				side = 10; break;
			}
		}
		if (side != 10) return;	// the new bound is useless
	}

	plane.isBound = true;
	bounds.push_back(plane);
	bounds_crossing_indexes.push_back(set<int>());
	int bind = bounds.size() - 1;
	if (bounds.size() < dim) return;		// there can be no crossings

	// find crossings
	hv[0] = plane;
	set<int> crossed_bounds;
	for (Index index = Index(dim - 1, bounds.size() - 1); index; index++){
		for (int k = 0; k < dim - 1; k++){
			hv[k + 1] = bounds[index[k]];
		}
		hs.get(hv);
		Point cp = hs.crossing_point();
		if (cp.is_nil())
			continue;

		bool acceptPoint = true;
		for (int i = 0; i < bounds.size() - 1; i++){
			if (index.has(i)) continue;
			if (bounds[i].side(cp) == -1){
				acceptPoint = false;
				break;
			}
		}
		if (!acceptPoint)
			continue;

		crossing_points.push_back(cp);
		int cind = crossing_points.size() - 1;
		used_crossing_points.insert(cind);
		bounds_crossing_indexes[bind].insert(cind);
		for (int k = 0; k < dim - 1; k++){
			bounds_crossing_indexes[index[k]].insert(cind);
			crossed_bounds.insert(index[k]);
		}
	}

	if (initial) return;

	set<int> remove_cp;
	for (set<int>::iterator c = used_crossing_points.begin(); c != used_crossing_points.end(); c++){
		if (!bounds_crossing_indexes[bind].count(*c)/*not my own*/){
			if (plane.side(crossing_points[*c]) != 1)
				remove_cp.insert(*c);
		}
	}
	set<int> remove_b;
	for (int i = 0; i < bounds.size(); i++){
		if (crossed_bounds.count(i)){
			remove_subset(bounds_crossing_indexes[i], remove_cp);
		}
		else
		if (is_subset(remove_cp, bounds_crossing_indexes[i])){
			remove_b.insert(i);
		}  
	}

	remove_subset(used_crossing_points, remove_cp);
	for (set<int>::reverse_iterator b = remove_b.rbegin(); b != remove_b.rend(); b++){
		int e = *b;
		bounds.erase(bounds.begin() + e);
		bounds_crossing_indexes.erase(bounds_crossing_indexes.begin() + e);
	}
}

// remove all unused crossing points, refresh bounds-crossings links
void clearBounds(vector<Point>& crossing_points, set<int>& used_crossing_points, vector<set<int> >& bounds_crossing_indexes){
	vector<int> new_ind(crossing_points.size());
	int i = 0;
	for (set<int>::iterator cp = used_crossing_points.begin(); cp != used_crossing_points.end(); cp++, i++){
		if (i != (*cp))
			crossing_points[i] = crossing_points[*cp];
		new_ind[*cp] = i;
	}
	crossing_points.resize(i);
	for (i = 0; i < bounds_crossing_indexes.size(); i++){
		set<int> new_set;
		for (set<int>::iterator cp = bounds_crossing_indexes[i].begin(); cp != bounds_crossing_indexes[i].end(); cp++){
			new_set.insert(new_ind[*cp]);
		}
		bounds_crossing_indexes[i] = new_set;
	}
	used_crossing_points.clear();
	for (i = 0; i < crossing_points.size(); i++){ used_crossing_points.insert(i); }
}

#ifdef _MSC_VER 
#ifdef DEEPDEBUG
OjaData OjaData::S = OjaData();
#endif 
#endif

OjaPoint OjaData::medianFollowIntersectionLinesBounded()	//AP
{
	int fail_count = 0;
	int counter = 1;

	//if(verbose)
	//	cout << "Max. search lines " << max_searchlines << endl;

#ifdef _MSC_VER 
#ifdef DEEPDEBUG
	//  used to test if the objective function equals the real value in DotSet::min
//	LOG("Initializing Reference Data for tesing");
	OjaData::S = *this;
	for (int i = 0; i < dim(); i++)
	for (int j = 0; j < size(); j++)
	OjaData::S[j][i] = (*this)[j][i];
	OjaData::S.generate_hyperplanes();
	/**/
#endif 
#endif

#ifdef TO_FILE
	stringstream filename; 
	filename << "D:\\OjaExperiments\\" << dim() << " " << size() << " (" << (*this)[0][0] << ") " << volume << ".txt";
	ofstream fout(filename.str().c_str());	

	FLOG("d: " << dim()); FLOG("n: " << size()); FLOG("volume: " << volume); FLOG("");
#endif
	clock_t begin = clock();

	/* 1.*/
//	LOG("Generating hyperplanes");
	generate_hyperplanes();

	clock_t hp_generated = clock();

	LOGTIME("hp generated " << double(hp_generated - begin) / CLOCKS_PER_SEC);

	/* 1.1 Get bounds */

#pragma region add initial bounding box

	vector<Hyperplane> bounds;
	vector<Point> crossing_points;
	set<int> used_crossing_points;
	vector<set<int> > bounds_crossing_indexes;
	{
		valarray<double> v(dim()), s, g;	for (int i = 0; i < dim(); i++) v[i] = 0;
		Hyperplane h; h.isBound = true;
		Point sp, grad;
		for (int i = 0; i < dim(); i++){
			s = v; g = v; // clean
			s[i] = max()[i]; sp = Point(s);
			g[i] = -1;     grad = Point(g);
			h.get(sp, grad);
			addBound(h, true, bounds, crossing_points, used_crossing_points, bounds_crossing_indexes, dim());
			s[i] = min()[i]; sp = Point(s);
			g[i] = +1;     grad = Point(g);
			h.get(sp, grad);
			addBound(h, true, bounds, crossing_points, used_crossing_points, bounds_crossing_indexes, dim());
		}
	}
#pragma endregion

	Point bmin = min(), bmax = max(), bmid = 0.5*(bmin+bmax);
	double initvol = getVolume(bmin, bmax);
	double volume;

#ifdef REMOVED_PARTS

if (false)	// Approach B go along the gradients
{
#pragma region go along the gradients (any start point)
//	LOG("Building gradients");
	//Going along the gradients as long as they are getting smaller
	Point p = average();	// starting point
	double d = DOUBLE_XMAX;
	vector<Point> grads;
	grads.push_back(p);

	while (true){
		Point g = oja_rank(p);
		Point p1 = p - g;
		double d1 = dist(p, p1);

		if (d1 == 0){
			OjaPoint T(this);
			T.set_location(p);
//			LOG("Zero gradient found in point " << p);
			return T;
		}

		grads.push_back(p1);

		if (d1 >= d)
			break;

		p = p1;
		d = d1;

#pragma endregion

#pragma region get bounds		
//		LOG("Building bounds");
		Hyperplane b1;
		for (int i = grads.size() - 2; i >= 0; i--){	// here: add ADD found gradients as bounds. other variants are possible
			Point g = -(grads[i] - grads[i + 1]);

			b1.get(grads[i], g);	// hp through point with gradient
			b1.isBound = true;

//			LOG("Gradient from " << grads[i]);
			addBound(b1, false, bounds, crossing_points, used_crossing_points, bounds_crossing_indexes, dim());
		}
	}
#pragma endregion
}

if (false) //Approach D build bounds at once, reduce length
{
#pragma region go along the gradients (any start point)
	FLOG("Building gradients");
	//Going along the gradients as long as they are getting smaller
	Point p = average();	// starting point
	double d = DOUBLE_XMAX;
	Hyperplane b1;
	int k = 1;
	int gf = 0;

	do{
		Point g = oja_rank(p);
		Point p1 = p - g;
		double d1 = dist(p, p1);

		if (d1 == 0){
			OjaPoint T(this);
			T.set_location(p);
			FLOG("Zero gradient found in point " << p);
			return T;
		}

		if (d1 > d){
			gf++;
			if (gf == dim()) { d *= 0.7; gf = 0; }
			g *= d / d1;

			d1 = d;
			p1 = p - g;
		}
		else { gf = 0; }

		g = -g;
		b1.get(p, g);	// hp through point with gradient
		b1.isBound = true;

		//		LOG("Gradient from " << p);
		addBound(b1, false, bounds, crossing_points, used_crossing_points, bounds_crossing_indexes, dim());

		bounded_min_max(crossing_points, used_crossing_points, bmin, bmax, bmid, dim());
		volume = getVolume(bmin, bmax);

//		LOG(k << " Min " << bmin << "; Max " << bmax << "; size: " << bmax - bmin << " (" << volume / initvol << "); bounds: " << bounds.size());
		FLOG(k << " Min " << bmin << "; Max " << bmax << "; size: " << bmax - bmin << " (" << volume / initvol << "); bounds: " << bounds.size());
		k++;

		p = p1;
		d = d1;
	} while (volume / initvol > this->volume /*desired volume*/ && k < 100/*max cuts*/);
#pragma endregion
}


if (false) //Approach D build bounds at half of gradient
{
#pragma region go along the gradients (any start point)
	FLOG("Building gradients");
	//Going along the gradients as long as they are getting smaller
	Point p = average();	// starting point
	Hyperplane b1;
	int k = 1;
	int gf = 0;

	do{
		Point g = oja_rank(p);
		Point p1 = p - g;

		if (p1 == p){
			OjaPoint T(this);
			T.set_location(p);
			FLOG("Zero gradient found in point " << p);
			return T;
		}

		Line l; l.get_through(p, p1);

		p1 = p - 1e6*g;
		for (int i = 0; i < bounds.size(); i++){
			if (bounds[i].side(p1) != 1){
				p1 = l.intersect(bounds[i]);
			}
		}

		g = -g;
		b1.get(p, g);	// hp through point with gradient
		b1.isBound = true;
		p = 0.5*(p + p1);

		//		LOG("Gradient from " << p);
		addBound(b1, false, bounds, crossing_points, used_crossing_points, bounds_crossing_indexes, dim());

		bounded_min_max(crossing_points, used_crossing_points, bmin, bmax, bmid, dim());
		volume = getVolume(bmin, bmax);

//		LOG(k << " Min " << bmin << "; Max " << bmax << "; size: " << bmax - bmin << " (" << volume / initvol << "); bounds: " << bounds.size());
		FLOG(k << " Min " << bmin << "; Max " << bmax << "; size: " << bmax - bmin << " (" << volume / initvol << "); bounds: " << bounds.size());
		k++;
	} while (volume / initvol > this->volume /*desired volume*/ && k < 100/*max cuts*/);
#pragma endregion
}

#endif // REMOVED_PARTS

clearBounds(crossing_points, used_crossing_points, bounds_crossing_indexes);

if (true)	// Approach A. reducing procedure
{
#pragma region reduce the bounded area

	if (this->volume > 1 || this->volume <= 0)
		this->volume = 0.0001;

	int k = 1;
	unsigned count = INT_MAX;
	do{
		bounded_min_max(crossing_points, used_crossing_points, bmin, bmax, bmid, dim());
		volume = getVolume(bmin, bmax);

//		LOG(k << " Min " << bmin << "; Max " << bmax << "; size: " << bmax - bmin << " (" << volume / initvol << "); bounds: " << bounds.size());
		FLOG(k << " Min " << bmin << "; Max " << bmax << "; size: " << bmax - bmin << " (" << volume / initvol << "); bounds: " << bounds.size());


		Point g = -oja_rank(bmid);
		if (g.length() < 1e-15) {
//			LOG("The median is found : " << bmid);
			FLOG("The median is found : " << bmid);
			OjaPoint pp; pp.set_location(bmid);
			return pp;
		}

		Hyperplane b1;
		b1.get(bmid, g);	// hp through point with gradient
		addBound(b1, false, bounds, crossing_points, used_crossing_points, bounds_crossing_indexes, dim());

		k++;
	} while (volume / initvol > this->volume /*desired volume*/ && k < 100/*max cuts*/);

	//	LOG("Set coord bounds");
	clearBounds(crossing_points, used_crossing_points, bounds_crossing_indexes);
	bounded_min_max(crossing_points, used_crossing_points, bmin, bmax, bmid, dim());
	set_bounded_min_max(bmin, bmax);

	FLOG(k << " Min " << bmin << "; Max " << bmax << "; size: " << bmax - bmin << " (" << volume / initvol << "); bounds: " << bounds.size());


#pragma endregion
}

if (method == FOLLOW_INTERSECTION_LINES_BOUNDED_APPROX)
{
	OjaPoint pp; pp.set_location(bmid);
	return pp;
}

/*
#ifdef TO_FILE
fout.close();
#endif
OjaPoint pp; pp.set_location(min());
return pp;
*/

#pragma region Remove hyperplanes lying out of bounds, add bounds and crossing points to search
	/* 1.2 Remove hyperplanes lying out of bounds */
	//	LOG("Remove hyperplanes lying out of bounds");
	for (int i = 0; i < planes; i++){
		int s = (*plane)[i].side(crossing_points[0]);
		if (s == 0)
			includedPlanes.insert(i);
		else
		for (int j = 1; j < crossing_points.size(); j++){
			int sj = (*plane)[i].side(crossing_points[j]);
			if (s != sj){	// min 2 points of the bounds lie on different sides => hyperplane crosses the bounded region
				includedPlanes.insert(i);
				break;
			}
		}
	}
	//	LOG("Hyperplanes left: " << includedPlanes.size() << " out of  " << planes << " (" << round(100.0*includedPlanes.size() / planes, 0.01) << "%)");
	FLOG("Hyperplanes left: " << includedPlanes.size() << " out of  " << planes << " (" << round(100.0*includedPlanes.size() / planes, 0.01) << "%)");
	/*
	for (int i = 0; i < bounds_crossing_indexes.size(); i++)
	{
		LOG(bounds_crossing_indexes[i]);
	}
	*/
#pragma endregion

#pragma region Add bounds to hyperplanes
	{
		vector<Line> axes;

		valarray<double> s(dim()), g;	for (int i = 0; i < dim(); i++) s[i] = 0;
		Point sp = Point(s), grad; Line l;
		for (int i = 0; i < dim(); i++){	// pure axes
			g = s; // clean
			g[i] = 1;     grad = Point(g);
			l.get_toward(sp, grad);
			axes.push_back(l);
		}

		for (int len = dim(); len>1; len--){// diagonal axes. first the ones diagonal to the most of the original axes
			Int32 b = 0;	for (int i = 0; i < dim(); i++)	b |= (1 << i);
			for (; b > 1; b--){
				Int32 v = b - ((b >> 1) & 0x55555555); v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
				if (len != ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24) continue; // count set bits, check if the sequence is of desired length 
				
				g = s; // clean
				for (int i = 0; i < dim(); i++){
					if (check_bit(b, i))
						g[i] = 1;
				}
				grad = Point(g);
				l.get_toward(sp, grad);
				axes.push_back(l);
			}
		}

		vector<Point> planes_points(bounds.size()*dim()); int pc = 0;
		vector<set<int> > bci(bounds.size());
		Point ip; int os = get_original_size();
		for (int i = 0; i < bounds.size(); i++){

			for (int a = 0; bci[i].size() < dim(); a++){
				ip = (bounds[i].intersect(axes[a]));
				if (ip.is_nil()) continue;

				planes_points[pc] = ip;
				bci[i].insert(os + pc);
				pc++;
			}
		}

		bounds_crossing_indexes = bci;

		add_bound_points(planes_points);
		// Add bounds to the planes list, add indexes
		for (int i = 0, p = planes; i < bounds.size(); i++)
			add_bound(bounds[i], bci[i]);
	
	}

#pragma endregion

	clock_t bb_generated = clock();
	//	LOG("bounds generated " << double(bb_generated - begin) / CLOCKS_PER_SEC);
	LOGTIME("bounds generated " << double(bb_generated - hp_generated) / CLOCKS_PER_SEC);

	/* 2. */
	OjaLine L(*this);
	IndexIdentifier Lid, Tid;
	
#pragma region Choosing initial line

	//	LOG("Choosing initial line");


	/* 8. */
	OjaPoint T(*this), hatT(*this);
	double hatD, D;

	vector<set<int> > hv(dim() - 1);
	HyperplaneSet hss(dim() - 1);	// higher precision then taking a line just by indexes
	Index index = Index(dim() - 1, bounds.size() - 1);
step3:
	for (; index && h.is_nil(); index++){
		for (int k = 0; k < dim() - 1; k++){
			hv[k] = bounds_crossing_indexes[index[k]];
			hss[k] = bounds[index[k]];
		}

		Line l; l.get(hss);		
		L.set(IndexSet(dim(), size(), hv), l);

		if (L.is_nil())
			continue;

		errif(L.is_nil(), "Could not find an initial line")
		Lid.get(L.index());
//		LOG(counter++ << " Trying " << Lid);

		DotSet tds(&L, h, h0);	// initialize the common coefficients
		if (h.is_nil())
			continue; // The first line lies out of bounds

#ifdef GRAPHICS
		set_line(L.line());
		wait_if_pause();
#endif
	}
	errif(L.is_nil(), "Could not find an initial line")
	errif(h.is_nil(), "Could not find an initial line")

	/*	vector<set<int> > hi(3);
		hi[0].clear();
		hi[0].insert(0);
		hi[0].insert(1);
		hi[0].insert(2);
		hi[0].insert(18);

		hi[1].clear();
		hi[1].insert(42);
		hi[1].insert(43);
		hi[1].insert(45);
		hi[1].insert(47);

		hi[2].clear();
		hi[2].insert(44);
		hi[2].insert(45);
		hi[2].insert(46);
		hi[2].insert(47);

		L.get(IndexSet(dim(), size(), hi)); */ // use this block to test any particular line

#pragma endregion

	//	LOG((counter = 1) << " Minimizing " << Lid); counter++;
	hatT = L.min(hatD);

#ifdef GRAPHICS
	add_orbit(hatT.location());
#endif
	Tid.get(hatT.index());
	//	LOG("  Minimum " << Tid << " (object function " << hatD << ")");

	/* 9. */
	set<IndexIdentifier> calL;

	calL.insert(Lid);

	/* 10. */
step10:
	T = hatT;
	D = hatD;

	/* 11. */

	int nL = Tid.sup_objects();

	/* 12. */

	if (nL > max_searchlines)
	{
	  //	LOG("Too many line possibilities (" << nL << ") at " << Tid);
		FLOG("Too many line possibilities (" << nL << ") at " << Tid);
		fail_count++;
		goto step3;
	}

	/* 15. */
	//	LOG("Generating " << nL << " lines");

	set<IndexIdentifier> calLprime;
	Tid.put_sup_objects(calLprime, 1);

	/* 16. */
	set<IndexIdentifier>::const_iterator i;

	for (i = calL.begin(); i != calL.end(); i++)
		calLprime.erase(*i);

	/* 17. */

	for (i = calLprime.begin(); i != calLprime.end(); i++)
	{
	  //		LOG("     "<<*i);
	}

	while (calLprime.size())
	{

		/* 18. */
		Point R = calLprime.size() == 1 ? Point() : gradient(hatT.location());
		double proj, proj_max = -1.0;
		JokerIdentifier J;
		IndexSet I;
		OjaLine Ltmp(*this);

		int realparts = 0;
		for (i = calLprime.begin(); i != calLprime.end(); i++)
		{
			int rp = (*i).real_partitions(get_original_size());
			if (rp > realparts) realparts = rp;
			if (rp < realparts) continue;
			J.get(*i);
			J.put(I);
			Ltmp.get(I);

			proj = calLprime.size() == 1 ? 0 :
				fabs((R | Ltmp.line().dir())*(1.0 / Ltmp.line().dir().length()));
			if (proj_max == -1.0 || proj > proj_max)
			{
				proj_max = proj;
				L = Ltmp;
			}
		}
		Lid.get(L.index());
		//		LOG("  Best direction " << Lid << " (projection " << proj_max << ")");

		/* 19. */
		//		LOG(counter << " Minimizing " << Lid); counter++;
#ifdef GRAPHICS
		set_line(L.line());
		wait_if_pause();
#endif
		OjaPoint p = L.min(hatD);
		if (p.is_nil()){
		  //			LOG("Line lies out of bounds");
			calLprime.erase(Lid);
			calL.insert(Lid);
			continue;
		}
		hatT = p;
#ifdef GRAPHICS
		add_orbit(hatT.location());
#endif
		Tid.get(hatT.index());
		//	LOG("  Minimum " << Tid << " (object function " << hatD << ")");

		/* 20. */
		calLprime.erase(Lid);
		calL.insert(Lid);

		/* 21. */
		if (hatD < D)
			goto step10;

		/* 24. */
	}

	/* 25. */
	/*if(verbose)
	{
	cout << calL.size() << " lines minimized" << endl;
	cout << fail_count << " failures" << endl;
	}*/

	//	LOG("! Minimum " << Tid << " (" << hatT.location() << ") (object function " << hatD << ")");
	FLOG("! Minimum " << Tid << " (" << hatT.location() << ") (object function " << hatD << ")");

//#ifdef TO_FILE
	clock_t end = clock();
	double elapsed_secst = double(end - hp_generated) / CLOCKS_PER_SEC;
	double elapsed_secs = double(end - bb_generated) / CLOCKS_PER_SEC;

	LOGTIME("counter: " << counter-1 << "\nTime total: " << elapsed_secst << "\nTime count: " << elapsed_secs);

//	fout.close();
//#endif

	return T;
}

void Hyperplane::get(Point& point, Point& grad){
	set_dim(point.dim());
	for (int i = 1; i<cofs; i++)
		cof[i] = grad.coord(i-1);
	cof[0] = -(point | grad);
}

Hyperplane OjaData::getBoundingHyperplane(Point& point){
	Hyperplane h;
	Point grad = oja_rank(point);
	h.get(point, grad);
	return h;
}
