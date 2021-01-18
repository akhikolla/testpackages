#pragma once

#include "SubsetGenerator.h"
#include "Matrix2D.h"
#include <vector>


using std::vector;

class PipeElement
{
protected:
	PipeElement *next;
public:
	PipeElement(void)
	{
		next=0;
	};
	
	void isFollowedBy(PipeElement *p)
	{
		next=p;
	};
	virtual void reset(void)
	{
	};
};

class Matrix2DPipeElement:public PipeElement
{
public:
	virtual void compute(Matrix2D* m2,Matrix2D* mReduced)=0;
};

class doubleVecPipeElement:public PipeElement
{
public:
	virtual void compute(double *value, int anz)=0;
};

template <class T> class TemplatePipeElement:public PipeElement
{
public:
	virtual void compute(T &value)=0;
};

template <class T> class SumPipe:public TemplatePipeElement<T>
{
	T sum;
	using PipeElement::reset;	// Without this the clang throws a warning
public:
	void compute(T &a)	{	sum+=a;	};
	void reset(T a)	{  sum=a; };
	T getResult(void){	return sum;	};
};

template <class T> class AbsoluteValuePipe:public TemplatePipeElement<T>
{
	T f;
public:
	void compute(T &a)
	{
		f=fabs(a);
		((TemplatePipeElement<double>*)this->next)->compute(f);
	};
};

template <class T> class StoreItPipe:public TemplatePipeElement<T>
{
	T value;
public:
	void compute(T &a) {   value=a;	};
	T getResult(void)	{	return value;	};
};

class DetPipe:public doubleVecPipeElement
{
	int* validRows;
	Vector *mu;
	int DIM;
public:
	virtual inline ~DetPipe()
	{
	delete[] validRows;
	};

	DetPipe(int _DIM):doubleVecPipeElement()
	{
		DIM=_DIM;
		validRows = new int[DIM+1];
	};
	void compute(double *value, int anz)
	{
		double f=value[0];

		double sign=-1;
		for(int i=1;i<anz;i++)
		{
			f=f+sign*mu->getValue(i-1)*value[i];
			sign*=(-1);
		}

		((TemplatePipeElement<double>*)next)->compute(f);
	};
	void setMu(Vector *_mu)
	{
		mu=_mu;
	};
};


class GradientPipe:public doubleVecPipeElement
{
	int* validRows;
	Vector *mu;
	Vector *nabla;
	int DIM;

    DetPipe compDet;
	StoreItPipe<double> detValue;
public:
        virtual inline ~GradientPipe()
   	{
	delete[] validRows;
	delete nabla;
	};

	GradientPipe(int _DIM):doubleVecPipeElement(),compDet(_DIM)
	{
		DIM=_DIM;
		validRows = new int[DIM+1];
		nabla=new Vector(DIM);
		compDet.isFollowedBy(&detValue);
	};
	void compute(double *value, int anz)
	{
		compDet.compute(value,anz);
		double det=detValue.getResult();
 
		double sign = 0; 
		if (det < 0)
			sign = -1.0;
		else if (det > 0)
			sign = 1.0;

		for(int i=0;i<nabla->getSize();i++)
			nabla->setValue(i,0);

		double f;

		for(int i=1;i<anz;i++)
		{
			f=sign*value[i];
			nabla->setValue(i-1,f);
			sign*=(-1);
		}
 
		((SumPipe<Vector>*)next)->compute(*nabla);
	};
	void setMu(Vector *_mu)
	{
		mu=_mu;
		compDet.setMu(_mu);
	};
};

class SubDetPipe:public Matrix2DPipeElement
{
	double *value;
	int* validRows2;
	int DIM;
public:
	SubDetPipe(int _DIM)
	{
		DIM=_DIM;
		value=new double[DIM+1];
		validRows2 = new int[DIM+1];
	};
	~SubDetPipe(void)
	{
		delete[] validRows2;
		delete[] value;
	};
	void compute(Matrix2D* m2,Matrix2D* mReduced)
	{
//		FILE * debugDatei;  
	       /*
		Die Determinanten werden zu test zwecken in eine Datei geschrieben 
		um sie dann anschlie�end in maple zu vergleichen
			*/
  		
		
		for(int i = 0; i < mReduced->getNumberRows(); i++)
		{
			int k = 0;
			for(int j = 0; j < mReduced->getNumberRows(); j++)
			{
				if((i) != j)
				{
					validRows2[k] = j;
					k++;
				}
			}
			
  			value[i]=mReduced->reduced_determinant_pivoting(validRows2, mReduced->getNumberRows() - 1, mReduced->getNumberColumns());
  			
  			
		}
		((doubleVecPipeElement*)next)->compute(value, m2->getNumberRows());
		
		
		
	};
};


class RepeaterPipe:public doubleVecPipeElement
{
	vector<int> anzahl;
	vector<double*> tab;
public:
	virtual inline ~RepeaterPipe()
	{
	for(unsigned int i=0;i<tab.size();i++)
	delete[] tab[i];
 	};

	void compute(double *value, int anz)
	{
		anzahl.push_back(anz);
		double *vec=new double[anz];

		for(int i=0;i<anz;i++)
		   vec[i]=value[i];

		tab.push_back(vec);

		((doubleVecPipeElement*)next)->compute(value, anz);
	};
	void reCompute(void)
	{
		for(unsigned int i=0;i<tab.size();i++)
  		   ((doubleVecPipeElement*)next)->compute(tab[i], anzahl[i]);
	};
};



class SubsetPipe:public PipeElement
{
	SubsetGenerator* g;
	Matrix2D *points;
	Matrix2D* m2;
	Matrix2D* mReduced;
	int DIM;
public:
	SubsetPipe(SubsetGenerator* _g, Matrix2D *_points)
	{
		g=_g;
		points=_points;
		DIM=points->getNumberRows();
		m2 = new Matrix2D(DIM+1,DIM+1);
		mReduced = new Matrix2D(DIM+1,DIM+0);
	};
	~SubsetPipe(void)
	{
		delete m2;
		delete mReduced;
	};
	void compute(Vector *mu)
	{
		g->reset(m2, mReduced, *points, *mu);

		while(!g->isFinished())
		{
			g->getNextSubset(m2,mReduced, *points);
			((Matrix2DPipeElement*)next)->compute(m2,mReduced);
		}
	}
};


 
class ComputeObjectiveFunction
{
	Vector mu;
	Matrix2D* points;
	SubsetPipe subset;
	DetPipe det;	
	SubDetPipe subDet;
	RepeaterPipe repeat;
	AbsoluteValuePipe<double> absV;
	SumPipe<double> sum;

	bool firstRun;
	bool STORE_SUBDETERMINANTS;
public:
	ComputeObjectiveFunction(SubsetGenerator *g, Matrix2D* _points,bool _STORE_SUBDETERMINANTS)
		:subset(g,_points)
		,det(_points->getNumberRows())
		,subDet(_points->getNumberRows())
	{
		points=_points;
		subset.isFollowedBy(&subDet);
		
		STORE_SUBDETERMINANTS=_STORE_SUBDETERMINANTS;
		if (STORE_SUBDETERMINANTS)
		{
		  subDet.isFollowedBy(&repeat);
		  repeat.isFollowedBy(&det);
		}
		else
		{
		  subDet.isFollowedBy(&det);
		}
		det.isFollowedBy(&absV);
		absV.isFollowedBy(&sum);
		firstRun=true;
	};
    double compute(Vector *mu)
	{
		subset.reset();
		subDet.reset();
		det.reset();
		repeat.reset();
		sum.reset(0);

		det.setMu(mu);

		if (firstRun || ! STORE_SUBDETERMINANTS)
		{
			subset.compute(mu);
  			firstRun=false;
		}
		else
		{
			repeat.reCompute();
		}

		double f=sum.getResult();

		return f;
	};
};


class ComputeNabla
{
	Vector mu;
	Matrix2D* points;
	SubsetPipe subset;
	GradientPipe nabla;
	SubDetPipe subDet;
	RepeaterPipe repeat;

	SumPipe<Vector> sum;

	bool firstRun;
	bool STORE_SUBDETERMINANTS;
public:
	ComputeNabla(SubsetGenerator *g, Matrix2D* _points, bool _STORE_SUBDETERMINANTS)
		:subset(g,_points)
		,nabla(_points->getNumberRows())
		,subDet(_points->getNumberRows())
	{
		points=_points;
		subset.isFollowedBy(&subDet);
		STORE_SUBDETERMINANTS=_STORE_SUBDETERMINANTS;
		if (STORE_SUBDETERMINANTS)
		{
		  subDet.isFollowedBy(&repeat);
		  repeat.isFollowedBy(&nabla);
		}
		else
		{
		  subDet.isFollowedBy(&nabla);
		}
		nabla.isFollowedBy(&sum);
		firstRun=true;
	};
    Vector compute(Vector *mu)
	{
		subset.reset();
		subDet.reset();
		nabla.reset();
		repeat.reset();
		Vector null(mu->getSize());

		sum.reset(null);

		nabla.setMu(mu);

		if (firstRun || ! STORE_SUBDETERMINANTS)
		{
			subset.compute(mu);
  			firstRun=false;
		}
		else
		{
			repeat.reCompute();
		}

		Vector f=sum.getResult();

		return f;
	};
};


