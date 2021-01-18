/* 
 *  algorithms.cpp : Algorithms for computing exact or approximate Oja Median.
 *  Copyright (C) 2000,2001 Tommi Ronkainen
 *  CVS $Id: algorithms.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "lattice.h"
#include "random.h"
#include "data.h"
#include "simplex.h"
#include "line.h"
#include "oja_geometry.h"
#include "global.h"
#include "matrix_wrapper.h"
#include "stl_tools.h"
#ifdef GRAPHICS
#include "graphics.h"
#endif

#include <ctime>
#include <fstream>
#include <sstream>

/*
defined in "global.h"
#define LOG(a) {if(debug)cout << a << endl;}
//#define TRACE(a) {if(trace)cerr << a << endl;}
//#define LOGIF(a,b) {if(debug && (a))cout << b << endl;}
//#define LOG(a) {}
#define TRACE(a) {}
#define LOGIF(a,b) {}
#define FLOG(a) {}//{if(0)fout << a << endl;}
*/

using namespace std;

bool unique_median;

//  Default values
//  ==============

#define LATTICE_GRIDSIZE 5
/* Kuinka monta hilapistett� alussa per dimensio */
#define SAMPLE_SIZE_MULTIPLIER 3
/* Hypertasojen m��r�n kasvukerroin */
#define SAMPLE_GIVEUP_LIMIT 10
/* Kuinka moneen hypertasosetin j�lkeen luovutetaan */

using namespace std;

// Median algorithms
// =================

OjaPoint OjaData::median()
{
	matrix dummy;

	switch(method)
    {
	  case EVAL_ALL_POINTS:
		  return medianEvalAllPoints();
	  case GRADIENT_DESCENT:
		  return medianGradientDescent();
	  case FOLLOW_INTERSECTION_LINES:
		  return medianFollowIntersectionLines();
	  case FOLLOW_INTERSECTION_LINES_BOUNDED:
	    return medianFollowIntersectionLinesBounded();  // This is a test, if this return fixed the valgrind message on this line! Check with Oleksii if the return is correct!!! (DF)
	  case FOLLOW_INTERSECTION_LINES_BOUNDED_APPROX:
		  return medianFollowIntersectionLinesBounded();
	  case BRUTE_FORCE:
		  return medianBruteForceSearch();
	  case LATTICE_APPROX:
		  return medianLatticeApprox();
	  case LATTICE_APPROX2:
		  return medianLatticeApprox2();
	  case LATTICE_APPROX3:
		  return medianLatticeApprox3();
	  case BOOTSTRAP:
		  return medianBootstrap(dummy,100);
	  case SIMPLEX_APPROX:
		  return medianSimplexApprox();
	}
    
    return OjaPoint(*this);
}

OjaPoint OjaData::medianAtDataPoints(double& min)
{
    errif(size()==0,"OjaData::medianAtDataPoints: no data");
    Point med;
    double o;
	
    unique_median = true;
    min = oja((*this)[0]);
    med = (*this)[0];
    for(int i=1; i<size(); i++)
    {
		TRACE((*this)[i]);
		
		o = oja((*this)[i]);
		if(o < min) 
		{
			min = o;
			med = (*this)[i];
			unique_median = true;
		}
		else if (o == min && (*this)[i] != med)
			unique_median = false;
    }	
    
	OjaPoint p(this);
	p.set_location(med);

    return p;
}

OjaPoint OjaData::medianEvalAllPoints()
{
    Point med;
    double o,min;
    
    med=medianAtDataPoints(min).location();
	
    IndexSet I(dim(),dim(),size());
    HyperplaneSet H(dim());
    Point P,mn,mx;
	long total=0,eval=0;
   
    mn = this->min();
    mx = this->max();
    while(I)
    {
		H.get(*this,I);// BUG: pit�isi k�ytt�� valmiita mahdollisuuksien mukaan

		P = H.crossing_point();
		total++;
		
		if(!P.is_nil() && P.in_box(mn,mx))
		{
			TRACE(P);
			o = oja(P);
			eval++;
			if(o < min) 
			{
				min = o;
				med = P;
				unique_median = true;
			}
			else if (o==min && P != med)
				unique_median = false;
		}
		
		I++;
    }
	
	//if(verbose)
	//	cout << "Evaluation done: " << size() << " data points and " << 
	//		eval << "/" << total << " crossing points" << endl;

    //if(!quiet && !unique_median)
	//	cout << "Median is not unique!" << endl;

	OjaPoint p(this);
	p.set_location(med);

    return p;
}

// ALGORITMI: Davidon-Fletcher-Powell, Kaisa Miettinen: "Optimointi", s.40
// BUG: Rankka kludge indeksien kanssa (huono luokkajako toimintojen osalta)
OjaPoint OjaData::medianGradientDescent()
{
	OjaPoint y(this);
    OjaLine L(*this);
	Line L0;
	Point G(dim()),oldG(dim()),oldY(dim()),p(dim()),q(dim());
	matrix D(dim(),dim());
	double diff;

	for(int i=0; i<dim(); i++)
		D(i,i)=1.0;
	
    generate_hyperplanes();
	
    L.get_random_through(center_index());
	oldY=(*this)[center_index()];
	oldG=gradient(oldY);
	
	for(int iter=1;;iter++)
	{
		LOG("Minimizing " << L.line());
		// BUG: minimointia voitaisiin nopeuttaa, koska tiedet��n minimin sijaitsevan
		//      negatiivisen gradientin suunnassa
		y = L.min();
		LOG("  Min. at " << y.location());
		
		G = gradient(y.location());
		LOG("  Gradient length " << G.length());
		diff=y.location().dist(oldY);
		LOG("  Result moved " << diff);
		
		if(diff < 0.0000000000001 && iter > 1) // BUG: kiinte� raja
			break;

		G=D*G;
		p=(y.location()-oldY);
		oldY=y.location();
		q=(G-oldG);
		oldG=G;
		
		D-= (D * ::covariance(q,q) * D) * (1.0/(q|(D*q)));		
		D+= ::covariance(p,p) * (1.0/(p|q));
		
		LOG("p="<<p);
		LOG("q="<<q);
		LOG("D="<<D);
		
		L0.get_toward(y.location(),D*G);
		
		L.set(L.index(),L0);
	}
	
    return y;
}

void OjaData::brute_force_search(OjaLine& L,OjaPoint& p,double hi_score)
{
	errif(p.is_nil(),"OjaData::brute_force_search: point is nil");
	
	OjaPoint q(*this);
	double oja;
	
	for(;;)
	{
		LOG("Searching at: " << p);
		LOG("  Hiscore is: " << hi_score);
		L = p.scan_all_routes(q,oja,hi_score);
		LOG("  Minimum on: " << L);
		if(oja < hi_score)
			hi_score=oja;
		
		LOG("  Minimum f(" << q << ") = " << oja);
 		if(L.is_nil())
 			break;
		p=q;
	}
}

OjaPoint OjaData::medianBruteForceSearch()
{
    OjaPoint p(*this);
    OjaLine L(*this);
	
    generate_hyperplanes();

    L.get_random_through(center_index());

	LOG("Starting at: " << L);
	p = L.min();

	brute_force_search(L,p,oja(p.location()));
	
	return p;
}

OjaPoint OjaData::medianFollowIntersectionLines()
{
	int fail_count=0;
	int counter = 1;
	//if(verbose)
	//	cout << "Max. search lines " << max_searchlines << endl;

	/*
	stringstream filename;
	filename << "D:\\OjaExperiments\\" << dim() << " " << size() << " (" << (*this)[0][0] << ").txt";
	ofstream fout(filename.str().c_str());
	*/
	clock_t begin = clock();
	FLOG("d: " << dim()); FLOG("n: " << size()); FLOG("volume: 100"); FLOG("");


	/* 1.*/
//XXX	LOG("Generating hyperplanes");
	generate_hyperplanes();

	clock_t hp_generated = clock();

	LOGTIME("hp generated " << double(hp_generated - begin) / CLOCKS_PER_SEC);
	LOGTIME("bounds generated 0");
	
	/* 2. */
    OjaLine L(*this);
	IndexIdentifier Lid,Tid;

  step3:
//XXX	LOG("Choosing initial line");
    L.get_random_through(center_index());
	Lid.get(L.index());
	if(Lid.dim() != 1)
	{
		LOG("Dimension of " << Lid << " is " << Lid.dim() << ", trying again");
		goto step3;
	}

	/* 8. */
	OjaPoint T(*this),hatT(*this);
	double hatD,D;

//XXX	LOG(counter << " Minimizing " << Lid);  counter++;
#ifdef GRAPHICS
	set_line(L.line());
	wait_if_pause();
#endif
	hatT=L.min(hatD);
#ifdef GRAPHICS
	add_orbit(hatT.location());
#endif
	Tid.get(hatT.index());
//XXX	LOG("  Minimum " << Tid << " (object function " << hatD << ")");

	/* 9. */
	set<IndexIdentifier> calL;

	calL.insert(Lid);

	/* 10. */
  step10:
	T=hatT;
	D=hatD;

	/* 11. */

	int nL;

	nL=Tid.sup_objects();

	/* 12. */

	if(nL > max_searchlines)
	{
//XXX		LOG("Too many line possibilities (" << nL << ") at " << Tid);
		fail_count++;
		goto step3;
	}

	/* 15. */
//XXX	LOG("Generating " << nL << " lines");
	
	set<IndexIdentifier> calLprime;
	Tid.put_sup_objects(calLprime,1);

	/* 16. */
	set<IndexIdentifier>::const_iterator i;

	for(i=calL.begin(); i!=calL.end(); i++)
		calLprime.erase(*i);

	/* 17. */

	while(calLprime.size())
	{
		
	/* 18. */
	Point R=gradient(hatT.location());
	double proj,proj_max=-1.0;
	JokerIdentifier J;
	IndexSet I;
	OjaLine Ltmp(*this);

	for(i=calLprime.begin(); i!=calLprime.end(); i++)
	{
		J.get(*i);
		J.put(I);
		Ltmp.get(I);
		
		proj=fabs((R|Ltmp.line().dir())*(1.0/Ltmp.line().dir().length()));
		if(proj_max==-1.0 || proj > proj_max)
		{
			proj_max=proj;
			L=Ltmp;
		}
	}
	Lid.get(L.index());
//XXX	LOG("  Best direction " << Lid << " (projection " << proj_max << ")");

	/* 19. */
//XXX	LOG(counter << " Minimizing " << Lid);  counter++;
#ifdef GRAPHICS
	set_line(L.line());
	wait_if_pause();
#endif
	hatT=L.min(hatD);
#ifdef GRAPHICS
	add_orbit(hatT.location());
#endif
	Tid.get(hatT.index());
//XXX	LOG("  Minimum " << Tid << " (object function " << hatD << ")");

	/* 20. */
	calLprime.erase(Lid);
	calL.insert(Lid);

	/* 21. */
	if(hatD < D)
		goto step10;
	
	/* 24. */
	}

	/* 25. */
	/*if(verbose)
	{
		cout << calL.size() << " lines minimized" << endl;
		cout << fail_count << " failures" << endl;
	}*/
	
//XXX	LOG("! Minimum " << Tid << " (" << hatT.location() << ") (object function " << hatD << ")");
	FLOG("! Minimum " << Tid << " (" << hatT.location() << ") (object function " << hatD << ")");

	clock_t end = clock();
	double elapsed_secst = double(end - hp_generated) / CLOCKS_PER_SEC;
	double elapsed_secs = double(end - hp_generated) / CLOCKS_PER_SEC;

	LOGTIME("counter: " << counter - 1 << "\nTime total: " << elapsed_secst << "\nTime count: " << elapsed_secs);

//	fout.close();

	return T;
}

double OjaData::confidence_size(Lattice* L)
{
	switch(lattice_measure)
	{
	  case LM_DIAMETER:
		  return L->volume_diameter_under(chi2_limit);
		  
	  case LM_EDGE_LENGTH_AVG:
		  return L->average_edge_length_under(chi2_limit);
		  
	  case LM_EDGE_LENGTH_MAX:
		  return L->max_edge_length_under(chi2_limit);
	}

	return 0.0;
}

OjaPoint OjaData::medianLatticeApprox()
{
	// Etsit��n ensimm�isen jaon tiheys
	double INITIAL_STEP;
	Point sz=max()-min();
	double sum=0.0;

	for(int i=0; i<dim(); i++)
		sum += sz[i]/LATTICE_GRIDSIZE;

	INITIAL_STEP=sum / dim();
	
	// Alustetaan muuttujat
	OjaPoint p(*this);
	Lattice L(min(),max(),INITIAL_STEP);
	Lattice *Lp=&L;
/*
	if(verbose)
	{
		cout << "Starting with " << L.points() << " points" << endl;
		cout << "First sample is " << set_size << " hyperplanes" << endl;
		cout << "Limit is " << chi2_limit << endl;
		cout << "Maximum lattice size at the end is " << epsilon << endl;
		cout << "Lattice measure " << (lattice_measure==1 ? "avg" : lattice_measure==2 ? "max" : "diam") << endl;
	}
*/	
	int sets=set_size; // Arvottavan setin koko
	int set=1; // Arvottavan setin j�rjestysnumero
	int n=0; // Hilaan lis�ttyjen hypertasojen m��r�
	int lattice_number=1; // Kuinka mones hila menossa.
	int lattice_points=L.points(); // Yhteenlaskettujen hilapisteiden m��r� kaikissa hiloissa.
	Index I(dim(),size()); // Arvotun hypertason indeksi
	Hyperplane H; // Arvottu hypertaso
	Point R(dim()); // V�liaikaismuuttuja gradientille
	matrix S(dim(),dim()); // Kovarianssimatriisi
	Point z; // V�liaikaismuuttja k�sitelt�v�lle hilapisteelle

	// Nollataan jokaisen solmun gradienttiestimaatti
	for(LatticeIterator j(L); j; j++)
		j.node().gradient=R;

#ifdef GRAPHICS
	set_lattice(&L);
#endif
	
	for(;;)
	{
		LOG("Sampling " << sets << " hyperplanes");
		LOG("Lattice have " << Lp->points() << " points");
		
#ifdef GRAPHICS
		message(sets);
#endif
		
		// Poimitaan satunnaisia hypertasoja 'SAMPLE_SIZE*4^set' kappaletta
		
		for(int sample=sets; sample; sample--)
		{
			// Arvontatapa riippuu siit�, ovatko hypertasot generoitu valmiiksi
			if(planes)
			{
				H=hyperplane(random(0,planes-1));
			}
			else
			{
				I.random();
				H.get(*this,I);
			}

			n++;
						
			// Lis�t��n hypertaso kerroksen gradienttiestimaatteihin
			for(LatticeIterator j(*Lp); j; j++)
			{
				z=j.point();
				R=j.node().gradient;
				
				R = (double(n-1)/double(n)) * R
					+((1.0/double(n)) * H.side(z)) * H.normal();
				
				j.node().gradient=R;
			}

			// Lis�t��n hypertaso kovarianssimatriisiin
 			S= S *(double(n-1)/double(n))
 				+ ::covariance(H.normal(),H.normal()) * (1.0/double(n));
			
// 			if(i > 2)
// 				for(LatticeIterator j(*Lp); j; j++)
// 				{
// 					R=j.node().gradient;	
// 					j.node().goodness=double(i) * (R | (S.inv() * R));
// 				}
		
// 			show_lattice();
// 			wait_if_pause();
		}

		// Lasketaan testisuure kaikille kerroksen hilapisteille
		matrix Sinv=S.inv();
		
		for(LatticeIterator j(*Lp); j; j++)
		{
  			R=j.node().gradient;
  			j.node().goodness=double(n) * (R | (Sinv * R));
		}

#ifdef GRAPHICS	
		show_lattice();
		wait_if_pause();
#endif
		
		// Tarkistetaan, joko luottamusalue on riitt�v�n pieni
		
		double boxsz;

		boxsz=confidence_size(Lp);

		LOG("We are now closer than " << boxsz);

		if(boxsz < epsilon)
			break;

		errif(set==SAMPLE_GIVEUP_LIMIT,"Giving up after " << set << " sets");

		// Tarkennetaan hilaa
		if(!Lp->focus_under(chi2_limit,true))
		{
			// Ei sopivia. Otetaan paras ja sen naapurit.
			Lp->focus_on(Lp->smallest_goodness(),Lp->smallest_goodness(),true);
		}
		
		// Muodostetaan alkuarvot tihennetylle hilalle
		Lp=&Lp->sub();
		Lp->update_from_parent();
#ifdef GRAPHICS	
		show_lattice();
		wait_if_pause();
#endif

		lattice_number++;
		lattice_points+=Lp->points();
		
		// M��ritell��n seuraavan kierroksen otoskoko
		sets *= SAMPLE_SIZE_MULTIPLIER;

		set++;
	}

#ifdef GRAPHICS	
	set_lattice(0);
#endif

	// Otetaan hilapiste ja sen l�hinaapurit, joiden avulla lasketaan
	// lineaarinen aproksimaatio gradientti funktiosta.

	SimpleIndex bestI=Lp->smallest_goodness();
	SimpleIndex add(bestI.dim(),-1,1);
	Data locations,gradients;
	int approx_nodes=0;
	
	while(add)
	{
		SimpleIndex i(bestI.dim(),-1,99999);
		i.get_values(bestI);
		i+=add;
		if((*Lp).in_lattice(i))
		{
			locations.enlarge((*Lp).point(i));
			gradients.enlarge((*Lp).node(i).gradient);
			approx_nodes++;
		}
		add++;
	}

	matrix A;
	Point b;

	linear_fit(A,b,locations,gradients);
	

	// Vastaus:
	p.set_location((-1.0) * A.inv() * b);
/*
	if(verbose)
	{
		Point unfixed=(*Lp).point(bestI);

		cout << "Total of " << n << " hyperplanes added to lattice" << endl;
		cout << "Final nodes in gradient approximation: " << approx_nodes << endl;
		cout << "Average size of lattice " << double(lattice_points)/lattice_number << endl;

		cout << "Difference without fix: " << (unfixed  - exact_median).length() << endl;
		cout << "Unfixed result: " << unfixed << endl;
		cout << "Improvement: " << ((unfixed  - exact_median).length()) - ((p.location() - exact_median).length()) << endl;
	}
*/	
	return p;
}

OjaPoint OjaData::medianLatticeApprox2()
{
	// Etsit��n ensimm�isen jaon tiheys
	double INITIAL_STEP;
	Point sz=max()-min();
	double sum=0.0;

	for(int i=0; i<dim(); i++)
		sum += sz[i]/LATTICE_GRIDSIZE;

	INITIAL_STEP=sum / dim();
	
	// Alustetaan muuttujat
	OjaPoint p(*this);
	Lattice L(min(),max(),INITIAL_STEP);
	Lattice *Lp=&L;
/*
	if(verbose)
	{
		cout << "Starting with " << L.points() << " points" << endl;
		cout << "First sample is " << set_size << " hyperplanes" << endl;
		cout << "Limit is " << chi2_limit << endl;
		cout << "Maximum lattice size at the end is " << epsilon << endl;
		cout << "Lattice measure " << (lattice_measure==1 ? "avg" : lattice_measure==2 ? "max" : "diam") << endl;
	}
*/
	int sets=set_size; // Arvottavan setin koko
	int set=1; // Arvottavan setin j�rjestysnumero
	int n=0; // Hilaan lis�ttyjen hypertasojen m��r�
	Index I(dim(),size()); // Arvotun hypertason indeksi
	Hyperplane H; // Arvottu hypertaso
	Point R(dim()); // V�liaikaismuuttuja gradientille
	matrix S(dim(),dim()); // Kovarianssimatriisi
	Point z; // V�liaikaismuuttja k�sitelt�v�lle hilapisteelle
	double D; // Arvotun hypertason "koko"
	double w; // Otospaino
	double D0=0.0; // Keskim��r�inen koko estimaatti
	double P0=0.0; // Keskim��r�inen samplaustn.
	int lattice_number=1; // Kuinka mones hila menossa.
	int lattice_points=L.points(); // Yhteenlaskettujen hila pisteiden m��r� kaikissa hiloissa.
	int total=1; // Arvottujen hypertasojen m��r�

	// Nollataan jokaisen solmun gradienttiestimaatti
	for(LatticeIterator j(L); j; j++)
		j.node().gradient=R;

#ifdef GRAPHICS
	set_lattice(&L);
#endif
	
	for(;;)
	{
		LOG("Sampling " << sets << " hyperplanes");
		LOG("Lattice have " << Lp->points() << " points");

#ifdef GRAPHICS
		message(sets);
#endif
		
		// Poimitaan satunnaisia hypertasoja 'SAMPLE_SIZE*4^set' kappaletta
		
		for(int sample=sets; sample; sample--)
		{
			double prob;
			
			do
			{
				I.random();

				H.get(*this,I);

				D=0;
				for(int j=0; j<dim(); j++)
					for(int k=0; k<dim(); k++)
					{
						double x_jk;
						x_jk=operator[](I[j])[k];
						
						D+=x_jk*x_jk;
					}
				
				D0=(double(total-1)/double(total)) * D0 + (1.0/double(total))*D;
				prob=1-exp(-D/D0);
				P0=(double(total-1)/double(total)) * P0 + (1.0/double(total))*prob;
				w=(P0/prob);

				total++;
			}
			while(Uniform(0,1) >= prob);

			n++;
	
			// Lis�t��n hypertaso kerroksen gradienttiestimaatteihin
			for(LatticeIterator j(*Lp); j; j++)
			{
				z=j.point();
				R=j.node().gradient;
				
				R = (double(n-1)/double(n)) * R
					+w * ((1.0/double(n)) * H.side(z)) * H.normal();
				
				j.node().gradient=R;
			}

			// Lis�t��n hypertaso kovarianssimatriisiin
 			S= S *(double(n-1)/double(n))
 				+ ::covariance(H.normal(),H.normal()) *
				(1.0/double(n)) * w;
		}

		// Lasketaan testisuure kaikille kerroksen hilapisteille
		matrix Sinv=S.inv();
		for(LatticeIterator j(*Lp); j; j++)
		{
  			R=j.node().gradient;	
  			j.node().goodness=double(n) * (R | (Sinv * R));
		}

#ifdef GRAPHICS	
		show_lattice();
		wait_if_pause();
#endif
		
		// Tarkistetaan, joko luottamusalue on riitt�v�n pieni
		
		double dist;

		dist=confidence_size(Lp);

		LOG("We are now closer than " << dist);

		if(dist < epsilon) break;

		errif(set==SAMPLE_GIVEUP_LIMIT,"Giving up after " << set << " sets");

		// Tarkennetaan hilaa
		if(!Lp->focus_under(chi2_limit,true))
		{
			// Ei sopivia. Otetaan paras ja sen naapurit.
			Lp->focus_on(Lp->smallest_goodness(),Lp->smallest_goodness(),true);
		}
		
		// Muodostetaan alkuarvot tihennetylle hilalle
		Lp=&Lp->sub();
		Lp->update_from_parent();
#ifdef GRAPHICS	
		show_lattice();
		wait_if_pause();
#endif

		lattice_number++;
		lattice_points+=Lp->points();
		
		// M��ritell��n seuraavan kierroksen otoskoko
		sets *= SAMPLE_SIZE_MULTIPLIER;

		set++;
	}
	
#ifdef GRAPHICS	
	set_lattice(0);
#endif
	p.set_location(Lp->point(Lp->smallest_goodness()));
/*
	if(verbose)
	{
		cout << "Total of " << total << " hyperplanes sampled" << endl;
		cout << "Total of " << n << " hyperplanes added to lattice" << endl;
		cout << "Sample ratio " << double(n)/total << endl;
		cout << "Average size of lattice " << double(lattice_points)/lattice_number << endl;
	}
*/
	return p;
}

/* Paramterit:

   set_size - kerralla arvottavien hypertasojen m��r�
   limit_percentage - luottamusprosentti
   epsilon - hilapistejaon maksimikoko
*/
OjaPoint OjaData::medianLatticeApprox3(list<Hyperplane>* store,list<Index>* idxstore,matrix* cov)
{
	// Etsit��n ensimm�isen jaon tiheys
	Point sz=(1.0/LATTICE_GRIDSIZE)*(max()-min());
	Point middlepoint=0.5*(max()+min());
	SimpleIndex splits(dim(),0,LATTICE_GRIDSIZE);
	
	double maxstep=-1.0;
	splits.fill(4);
	
	for(int i=0; i<dim(); i++)
		if(maxstep < 0 || maxstep < sz[i])
			maxstep=sz[i];
	for(int i=0; i<dim(); i++)
		sz[i]=(0.5*LATTICE_GRIDSIZE)*maxstep;

	FreeLattice L(middlepoint-sz,middlepoint+sz,splits); // Hila

	// Alustetaan muuttujat
	FreeLattice oldL(L); // Varakopio luottamusalueen palauttamiseen
	Index I(dim(),size()); // Apumuuttuja indeksin samplaukseen
	Hyperplane H; // Samplattu hypertaso
	Point z; // Apumuuttuja, kulloinenkin hilapiste
	Point R(dim()); // apumuuttuja gradientin laskentaan
	matrix S(dim(),dim()); // Kovarianssimatriisi
	matrix oldS(S); // Kovarianssimatriisin varakopio
	SimpleIndex bestI; // Kulloinkin parhaan pisteen indeksi
	int lattice_number=1; // Kuinka mones hila menossa.
	int lattice_points=L.points(); // Yhteenlaskettujen hila pisteiden m��r� kaikissa hiloissa.
	int restores=0; // Kuinka usein jouduttiin peruuttamaan.
	int restored_planes=0; // Kuinka monta planea s�mpl�ttiin turhaan
	int max_set_size=set_size*2; // Adaptiivisen moodin rajoitin, joka kasvaa hilakerroksittain

	/*
	if(verbose)
	{
		cout << "One sample is " << set_size << " hyperplanes" << endl;
		cout << "Limit is " << chi2_limit << endl;
		cout << "Maximum lattice size at the end is " << epsilon << endl;
		cout << "Starting with " << L.points() << " points" << endl;
	}
	*/
	
	// Nollataan jokaisen solmun gradienttiestimaatti
	for(LatticeIterator j(L); j; j++)
		j.node().gradient=R;
	
#ifdef GRAPHICS
	set_lattice(&L);
#endif
	
	int oldn,n=0; // Hypertasojen kokonaism��r�
	for(int cntr = 1, cntrmax = 5000;; cntr++)
	{
		oldL=L;
		oldS=S;
		oldn=n;
		
//		LOG("Lattice have " << L.points() << " points, sampling " << set_size << " planes. " << cntr);
	  
		// Poimitaan satunnaisia hypertasoja 'SAMPLE_SIZE3' kappaletta
		for(int sample=set_size; sample; sample--)
		{
			// Arvontatapa riippuu siit�, ovatko hypertasot generoitu valmiiksi
			if(planes)
			{       
				H=hyperplane(random(0,planes-1));
			}
			else
			{
				I.random();
				if(idxstore)
					(*idxstore).push_back(I);
				H.get(*this,I);
			}
			if(store)
				(*store).push_back(H);

			n++;
			
			// Lis�t��n hypertaso kerroksen gradienttiestimaatteihin, sellaisissa
			// pisteiss�, jotka eiv�t ole poissa luottamusalueelta
			for(FreeLatticeIterator j(L); j; j++)
			{
				z=j.point();
				R=j.node().gradient;
				
				R = (double(n-1)/double(n)) * R
					+((1.0/double(n)) * H.side(z)) * H.normal();
				
				j.node().gradient=R;
			}
			
			// Lis�t��n hypertaso kovarianssimatriisiin
 			S= S *(double(n-1)/double(n))
 				+ ::covariance(H.normal(),H.normal()) * (1.0/double(n));
			
		}

		// Lasketaan testisuure kaikille kerroksen hilapisteille
		matrix Sinv=S.inv();
		for(FreeLatticeIterator j(L); j; j++)
		{
  			R=j.node().gradient;	
  			j.node().goodness=double(n) * (R | (Sinv * R));
		}

		bestI=L.smallest_goodness();
		
		L.remove_nodes_over(chi2_limit);
//		LOG("We have now " << L.points() << " points left");

		if(L.points())
			bestI=L.smallest_goodness();
		
#ifdef GRAPHICS
		if(!store)
		{
			show_lattice();
			wait_if_pause();
		}
#endif

		if(L.points() <= 1 || cntr == cntrmax)  // if there is one point or the iterations limit is reached
		{
 			if(L.points()==0)
 			{
 				L=oldL;
				S=oldS;
				n=oldn;
 				LOG("Restoring old lattice");
				restores++;
				restored_planes+=set_size;
				
				if(adaptive)
					set_size=(set_size >= 2 ? set_size/2 : 1);
 			}
 			else
			{		
				if(L.box_average_edge_length() < epsilon)
					break;

				if (cntr == cntrmax)
				{
					LOG("Too many iterations. Generating a new grid around the best point.");
					cntrmax = 100;
				}
				cntr = 1;

				LOG("Lattice grid size " << L.box_average_edge_length());

				L.focus_on(bestI,bestI,true);
				lattice_number++;
				lattice_points+=L.points();

				if(adaptive)
				{
					max_set_size*=2;
				}

                // KOKEILU
				// Nollataan estimaatit
// 				Point R(dim());
// 				for(LatticeIterator j(L); j; j++)
// 					j.node().gradient=R;
// 				S=matrix(dim(),dim()); 
// 				n=1;
			}
			
#ifdef GRAPHICS
			if(!store)
			{
				show_lattice();
				wait_if_pause();
			}
#endif
		}
		else // L.points() > 1
		{
			if(adaptive && L.points() > 5 && set_size < max_set_size)
				set_size*=2;
		}
	}

#ifdef GRAPHICS	
	set_lattice(0);
#endif

	// Otetaan hilapiste ja sen kolme lupaavinta l�hinaapuria, joiden avulla lasketaan
	// lineaarinen aproksimaatio gradienttifunktiosta.

// 	SimpleIndex add(bestI.dim(),-1,1);
// 	Data locations(dim(),dim()+1),gradients(dim(),dim()+1);
// 	double* goodness;
// 	int found=0;
// 	goodness=new double[dim()+1];
// 	for(int i=0; i<dim()+1; i++)
// 		goodness[i]=-1.0;
	
// 	while(add)
// 	{
// 		SimpleIndex i(dim(),-1,99999);
// 		i.get_values(bestI);
// 		i+=add;
// 		if(L.in_lattice(i))
// 		{
// 			for(int k=0; k<dim()+1; k++)
// 				if(goodness[k]==-1.0 || goodness[k] > L.node(i).goodness)
// 				{
// 					for(int j=dim(); j>k; j--)
// 					{
// 						goodness[j]=goodness[j-1];
// 						locations[j]=locations[j-1];
// 						gradients[j]=gradients[j-1];
// 					}
					
// 					goodness[k]=L.node(i).goodness;
// 					locations[k]=L.point(i);
// 					gradients[k]=L.node(i).gradient;

// 					found++;
// 					break;
// 				}
// 		}
// 		add++;
// 	}
// 	delete goodness;

	
	// Vastaus:
	OjaPoint p(*this);
// 	if(found < dim()+1)
// 	{
// 		cerr << "warnign: unable to make final fix; not enough points" << endl;
		p.set_location(L.point(bestI));
// 	}
// 	else
// 	{
// 		matrix A;
// 		Point b;

// 		linear_fit(A,b,locations,gradients);
	
// 		p.set_location((-1.0) * A.inv() * b);
// 	}
	
//XXX	if(verbose)
//XXX	{
//XXX		Point unfixed=L.point(bestI);
		/*
		cout << "Total of " << n << " hyperplanes added to lattice" << endl;
		cout << "Average size of lattice " << double(lattice_points)/(lattice_number+restores) << endl;
		cout << restores << " times returned back to older lattice" << endl;
		cout << restored_planes << " useless hyperplanes sampled" << endl;

		if(!exact_median.is_nil())
			cout << "Difference without fix: " << (unfixed  - exact_median).length() << endl;
		cout << "Unfixed result: " << unfixed << endl;
		if(!exact_median.is_nil())
			cout << "Improvement: " << ((unfixed  - exact_median).length()) - ((p.location() - exact_median).length()) << endl;
		*/
//XXX	}
	
	if(cov)
		*cov=S;
	
	return p;
}

static int _weight(int *weight,const Index& I)
{
	int sum=weight[I[0]];
	for(int i=1; i<I.dim(); i++)
		sum *= weight[I[i]];

	return sum;
}

// Bootstrappaava versio, joka k�ytt�� valmiiksi laskettuja hypertasoja
OjaPoint OjaData::medianBootstrap(const list<Hyperplane>& store,const list<Index>& idxstore,const matrix& CovInv)
{
	int* w; // Bootstrap otoksen kertoimet kullekin datapisteelle
	int points=size(); // Havaintoaineiston koko

	// Muodostetaan bootstrap painokertoimet
	w=new int[points];
	errif(!w,"OjaData::medianBootstrap: out of memory");

	for(int i=0; i<points; i++)
		w[i]=0;
	for(int i=0; i<points; i++)
		w[random(0,points-1)]++;

	// Etsit��n ensimm�isen jaon tiheys
	double INITIAL_STEP;
	Point sz=max()-min();
	double sum=0.0;

	for(int i=0; i<dim(); i++)
		sum += sz[i]/LATTICE_GRIDSIZE;

	INITIAL_STEP=sum / dim();
	
	// Alustetaan muuttujat
	FreeLattice L(min(),max(),INITIAL_STEP); // Hila
	FreeLattice oldL(L); // Varakopio luottamusalueen palauttamiseen
	int oldN=0; // Varakopio hypertasojen m��r�st�
	Index I(dim(),size()); // Apumuuttuja indeksin samplaukseen
	list<Hyperplane>::const_iterator hyp=store.begin(); // Listan seuraavaksi k�ytett�v� hypertaso
	list<Index>::const_iterator idx=idxstore.begin(); // Listan seuraavaksi k�ytett�v� indeksi
	Hyperplane H; // Samplattu hypertaso
	Point z; // Apumuuttuja, kulloinenkin hilapiste
	Point R(dim()); // apumuuttuja gradientin laskentaan
	SimpleIndex bestI; // Kulloinkin parhaan pisteen indeksi
	
	// Nollataan jokaisen solmun gradienttiestimaatti
	for(LatticeIterator j(L); j; j++)
		j.node().gradient=R;

#ifdef GRAPHICS
	set_lattice(&L);
#endif
	
	int n=0;
	int totalw=0; // Painojen kumulatiivinen summa
	int wg;
	
	for(;;)
	{	
		oldL=L;
		oldN=n;
		
		// Poimitaan satunnaisia hypertasoja 'SAMPLE_SIZE3' kappaletta
		for(int sample=set_size; sample; sample--)
		{
			// Kaivetaan vanha taso tai lopetetaan, mik�li lopussa
			if(hyp==store.end())
			{
				LOG("Out of hyperplanes");
				goto all_done;
			}
			
			I=*idx;
			H=*hyp;
			idx++;
			hyp++;

			// Lis�t��n hypertaso kerroksen gradienttiestimaatteihin, sellaisissa
			// pisteiss�, jotka eiv�t ole poissa luottamusalueelta
			wg=_weight(w,I);

			if(wg)
			{				
				n++;
				
				for(FreeLatticeIterator j(L); j; j++)
				{
					z=j.point();
					R=j.node().gradient;
					
 					R = (1.0/double(totalw+wg)) * ((double(totalw) * R) + (H.side(z) * H.normal()));
//					R = (double(totalw)/double(totalw+wg)) * R + (double(wg)/double(totalw+wg)) * (H.side(z) * H.normal());
					
					j.node().gradient=R;
				}
				
				totalw+=wg;
			}
		} // for(int sample=set_size; sample; sample--)

		// Lasketaan testisuure kaikille kerroksen hilapisteille
		for(FreeLatticeIterator j(L); j; j++)
		{
  			R=j.node().gradient;	
  			j.node().goodness=double(n) * (R | (CovInv * R));
		}

		bestI=L.smallest_goodness();
		
		L.remove_nodes_over(chi2_limit);

		if(L.points())
			bestI=L.smallest_goodness();
		
#ifdef GRAPHICS	
		show_lattice();
		wait_if_pause();
#endif

		if(L.points() <= 1)
		{
 			if(L.points()==0)
 			{
 				L=oldL;
				n=oldN;
 			}
 			else
			{				
				if(L.box_diameter() < epsilon)
					break;
				
				L.focus_on(bestI,bestI,true);
			}
		}
		
	} // for(;;)

	// Siivotaan ja poistutaan

  all_done:

#ifdef GRAPHICS	
	set_lattice(0);
#endif
	
	OjaPoint med(*this); // Bootsrap otosta vastaava mediaani
	med.set_location(L.point(L.smallest_goodness()));

	delete[] w;
	
	return med;
}

OjaPoint OjaData::medianBootstrap(matrix& covariance,int how_many)
{
	matrix Cov;
	
	errif(planes,"OjaData::medianBootstrap: not implemented with -pH");
	
	list<Hyperplane> store;
	list<Index> idxstore;
	OjaPoint med,bootmed;

	//if(verbose)
	//	cout << "ESTIMATING OJA MEDIAN:" << endl;
	
	med=medianLatticeApprox3(&store,&idxstore,&Cov);
	//if(verbose)
	//	cout << "Median " << med.location() << endl;
	
	//if(verbose)
	//	cout << "ESTIMATING COVARIANCE:" << endl;

	Cov=Cov.inv();

	matrix ret(dim(),dim());

	for(int i=0; i<how_many; i++)
	{
		bootmed=medianBootstrap(store,idxstore,Cov);
		//if(verbose)
		//	cout << (i+1) << ". bootstrap " << bootmed.location() << endl;

		ret+=::covariance(bootmed.location()-med.location(),bootmed.location()-med.location());
	}

	ret=ret * (1.0 / double(how_many));
	
	covariance=ret;
	
	return med;
}

// Rakenne, joka s�ilytt�� yht� data pistett� vastaavan gradienttiestimaatin.
struct node
{
	int index; // datapisteen indeksi
	Point gradient; // gradienttiestimaatti
	double goodness; // testisuure

	node(int i,Point p,double g)
		{index=i; gradient=p; goodness=g;}
};

class chi2_limit_check : public unary_function<node,bool>
{
	double max;
	
public:
	chi2_limit_check(double m)
		{max=m;}
	bool operator()(const node& n) const
		{return n.goodness > max;}
};

bool operator<(const node& n1,const node& n2)
{
	return n1.goodness < n2.goodness;
}

OjaPoint OjaData::medianSimplexApprox()
{
	list<node> ok; // Lista j�ljell� olevista kelvollisista pisteist�.
	list<node>::iterator p; // Iteraattori edelliselle.
	Index I(dim(),size()); // Samplatun tason indeksi
	Hyperplane H; // Samplattu taso
	int n; // Kuinka monta tasoa on samplattu.
	int phase2=dim()*4; // Raja, jolloin siirryt��n yksitt�iseen samplaukseen.
	matrix S(dim(),dim()); // kovarianssimatriisin estimaatti.
	Point z,R; // Apumuuttujia
	chi2_limit_check test_fail(chi2_limit); // Testioperaattori
	Simplex Area(dim()); // d+1 parasta pistett�
/*	
	if(verbose)
	{
		cout << "Sample set is " << set_size << " hyperplanes" << endl;
		cout << "Limit is " << chi2_limit << endl;
		cout << "Phase 2 limit is " << phase2 << endl;
	}
*/
  start_again:

	n=0;
	ok.clear();
	
	// Alustetaan listaan kaikki pisteet
	for(int i=0; i<size(); i++)
		ok.push_back(node(i,Point(dim()),0.0));

	// P��iteraatio
	do
	{
		LOG("We have " << ok.size() << " points left, sampling " << set_size << " planes");
		
		// Poimitaan satunnaisia hypertasoja 'set_size' kappaletta
		for(int i=set_size; i; i--)
		{
			// Arvontatapa riippuu siit�, ovatko hypertasot generoitu valmiiksi
			if(planes)
			{
				H=hyperplane(random(0,planes-1));
			}
			else
			{
				I.random();
				H.get(*this,I);
			}
			n++;
			
			// Lis�t���n taso pisteiden gradienttiestimaatteihin
			for(p=ok.begin(); p!=ok.end(); p++)
			{
				z=(*this)[(*p).index];
				R=(*p).gradient;
				
				R = (double(n-1)/double(n)) * R
					+((1.0/double(n)) * H.side(z)) * H.normal();
				
				(*p).gradient=R;
			}
			
			// Lis�t��n hypertaso kovarianssimatriisiin
			S= S *(double(n-1)/double(n))
				+ ::covariance(H.normal(),H.normal()) * (1.0/double(n));
		}
		
		// Lasketaan testisuure kaikille j�ljell�oleville pisteille
		matrix Sinv=S.inv();
		for(p=ok.begin(); p!=ok.end(); p++)
		{
			R=(*p).gradient;
			(*p).goodness=double(n) * (R | (Sinv * R));
		}
		
		// Poistetaan rajan ylitt�neet.
		ok.remove_if(test_fail);
		
		// Tarkistetaan joko setti� on pienennetty tarpeeksi.
	}
	while(int(ok.size()) > phase2);

	if(int(ok.size()) < dim()+1)
	{
		phase2+=dim();
		
		//cerr << "warning: lost too many points, start again...(trying limit " << phase2 << ")" << endl;
		goto start_again;
	}

	int valid; // Kelvollisten pisteiden m��r�.
	int ptsin; // Pisteit� simpleksin sis�ll�.

	// J�lkimm�inen vaihe: ei poisteta en�� pisteit�.
	do
	{
		do
		{
			// Arvontatapa riippuu siit�, ovatko hypertasot generoitu valmiiksi
			if(planes)
			{
				H=hyperplane(random(0,planes-1));
			}
			else
			{
				I.random();
				H.get(*this,I);
			}
			n++;
			
			// Lis�t���n taso pisteiden gradienttiestimaatteihin
			for(p=ok.begin(); p!=ok.end(); p++)
			{
				z=(*this)[(*p).index];
				R=(*p).gradient;
				
				R = (double(n-1)/double(n)) * R
					+((1.0/double(n)) * H.side(z)) * H.normal();
				
				(*p).gradient=R;
			}
			
			// Lis�t��n hypertaso kovarianssimatriisiin
			S= S *(double(n-1)/double(n))
				+ ::covariance(H.normal(),H.normal()) * (1.0/double(n));
			
			// Lasketaan testisuure kaikille j�ljell�oleville pisteille
			matrix Sinv=S.inv();
			for(p=ok.begin(); p!=ok.end(); p++)
			{
				R=(*p).gradient;
				(*p).goodness=double(n) * (R | (Sinv * R));
			}
			
			// Lasketaan kelvolliset.
			valid=0;
			for(p=ok.begin(); p!=ok.end(); p++)
				if(!test_fail(*p))
					valid++;
		}
		while(valid > dim()+1);

		// J�rjestet��n paremmuusj�rjestykseen.
		ok.sort();
	
		// T�ytet��n simpleksi.
		p=ok.begin();
		for(int i=0; i<dim()+1; i++)
		{
//			cerr << (*this)[p->index] << endl;
			Area.set_column(i,(*this)[p->index]);
			p++;
		}
	
		// Tarkistetaan lopetusehto
		ptsin=0;
		for(p=ok.begin(); p!=ok.end(); p++)
			if(Area.contains((*this)[p->index]))
				ptsin++;

	}while(ptsin > dim()+1);

	//if(verbose && !exact_median.is_nil())
	//	cout << "Simplex contains median: " << Area.contains(exact_median) << endl;
	
	// Lasketaan keskipiste.
 	LOG("Result contains " << dim()+1 << " of " << ok.size() << " points");
	OjaPoint med(*this);

 	Point avg(dim());
 	p=ok.begin();
 	for(int i=dim()+1; i; i--)
 	{
 		avg+=(*this)[p->index];
 		p++;
 	}
 	avg=(1.0/(dim()+1)) * avg;
  	med.set_location(avg);
//	med.set_location((*this)[ok.begin()->index]);

	if(verbose)
	{
		//cout << "Final simplex" << endl;
		p=ok.begin();
		for(int j=dim()+1; j; j--,p++)
// 			if(trace)
// 				cerr << "  " << (*this)[p->index] << endl;
// 			else
// 				cout << "  " << (*this)[p->index] << "  gradient: " << p->gradient << "  test value: " << p->goodness << endl;
// 		
		// Etsit��n tarkkaa mediaania l�hin havainto.
		if(!exact_median.is_nil())
		{
			double dist=0.0;
			//int bestidx=0;
			for(int i=0; i<size(); i++)
				if(i==0 || exact_median.dist((*this)[i]) < dist)
				{
					dist=exact_median.dist((*this)[i]);
					//bestidx=i;
				}

		//	cout << "Distance to nearest data " << dist << endl;
			p=ok.begin();
		//	for(int j=dim()+1; j; j--,p++)
		//		if(p->index == bestidx)
		//			cout << "Nearest point belongs to simplex" << endl;
		}
	}
		
	return med;
}

