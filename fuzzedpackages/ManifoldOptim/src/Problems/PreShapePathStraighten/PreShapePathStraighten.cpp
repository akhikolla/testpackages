
#include "PreShapePathStraighten.h"

/*Define the namespace*/
namespace ROPTLIB{

	PreShapePathStraighten::PreShapePathStraighten(integer innumP, integer indim, integer innumC)//variable???????
	{
		numP = innumP;
		dim = indim;
		numC = innumC;
	};

	PreShapePathStraighten::~PreShapePathStraighten(void)
	{
	};

	double PreShapePathStraighten::f(Variable *x) const    //cost function  ?????????????????????
	{
		const double *Path_x = x->ObtainReadData();
		Vector *Grad_temp = x->ConstructEmpty(); //根据x size获取相同size
		SharedSpace *Temp = new SharedSpace(Grad_temp); //把临时变量放到sharedspace里
		double *temp = Grad_temp->ObtainWriteEntireData();


		integer stp = numC - 1;
		//integer total_N = numP*dim*numC;
		//double *temptemp = new double[numP*dim*numC];
		//double *trapz_temp = new double[numC];
		//double intv;
		double result;

		//OUTSTREAM << x->Getlength() << std::endl;
		//OUTSTREAM << egf->Getlength() << std::endl;


		//Compute Dalpha
		for (integer t = 0; t < numC; t++)
		{
			if (t != 0)
			{
				for (integer j = 0; j < dim; j++)
				{
					for (integer i = 0; i < numP; i++)
					{
						//OUTSTREAM << t*numP*dim+j*numP+i << std::endl;
						//OUTSTREAM << stp <<std::endl;
						//OUTSTREAM << Path_x[t*numP*dim+j*numP+i]<< std::endl;
						//OUTSTREAM << (t-1)*numP*dim+j*numP+i << std::endl;
						//OUTSTREAM << Path_x[(t-1)*numP*dim+j*numP+i] << std::endl;
						temp[t*numP*dim + j*numP + i] = static_cast<double> (stp)*(Path_x[t*numP*dim + j*numP + i] - Path_x[(t - 1)*numP*dim + j*numP + i]);
					}
				}
			}
			//Project c(tau/numC) into T_alpha(M)
			if (t == 0)
			{
				for (integer j = 0; j < dim; j++)
				{
					for (integer i = 0; i < numP; i++)
					{
						temp[j*numP + i] = 0.0;
					}
				}
			}
			else
			{
				Item_2(Path_x + t*numP*dim, numP, dim, temp + t*numP*dim);
			}
		}
		//ForDebug::Print("Dalpha", temp, numP, dim,numC);

		x->AddToTempData("Dalpha", Temp);

		result = Domain->Metric(x, Grad_temp, Grad_temp)*0.5;
		//    OUTSTREAM << result<<std::endl;

		//    //Compute Cost function
		//    dcopy_(&total_N, const_cast<double *> (temp), &GLOBAL::IONE, temptemp, &GLOBAL::IONE);
		//    for (integer i = 0; i < numC; i++) {
		//        trapz_temp[i] = InnerProd_Q(temptemp+i*numP*dim, temptemp+i*numP*dim, numP, dim);
		//    }
		//    
		//    intv = 1.0/(numC-1);
		//    result = 0.5*ElasticCurvesRO::Trapz(trapz_temp, numC, intv);
		//    //OUTSTREAM << "PreShapePathStraighten::f, TODO" << std::endl;//---
		return result;
	};


	void PreShapePathStraighten::EucGrad(Variable *x, Vector *egf) const   //Dalpha需要提前assign初值么
	{
		const SharedSpace *Temp = x->ObtainReadTempData("Dalpha");
		Temp->GetSharedElement()->CopyTo(egf);
	};



	void PreShapePathStraighten::Item_1(const double *q, integer innumP, integer indim, double *q_c)
	{
		integer NXD = innumP*indim;
		integer iter;
		double coeff, coeff1, coeff2, dq_abs;
		const double TOL = 1e-10;
		double *qnorm = new double[innumP];
		double *tt = new double[indim];
		double *temp = new double[indim];
		double *b = new double[indim*innumP*indim];
		double *J = new double[indim*indim];
		double *beta = new double[indim];
		double *dq = new double[innumP*indim];
		integer *IPIV = new integer[indim];
		integer INFO;

		dcopy_(&NXD, const_cast<double *> (q), &GLOBAL::IONE, q_c, &GLOBAL::IONE);
		coeff = 1.0 / std::sqrt(InnerProd_Q(q_c, q_c, innumP, indim));
		dscal_(&NXD, &coeff, q_c, &GLOBAL::IONE);

		for (integer i = 0; i < innumP; i++) {
			qnorm[i] = dnrm2_(&indim, q_c + i, &innumP);     //check!!!!!!!!!!
		}

		for (integer t = 0; t < indim; t++) {
			tt[t] = InnerProd_Q(q_c + t*innumP, qnorm, innumP, 1);
		}

		iter = 0;
		while (dnrm2_(&indim, tt, &GLOBAL::IONE) > TOL && iter < 100) {
			//*******
			for (integer t = 0; t < innumP; t++) {
				qnorm[t] = dnrm2_(&indim, q_c + t, &innumP);
				for (integer d = 0; d < indim; d++) {
					//find the d-th temp
					for (integer i = 0; i < indim; i++) {
						temp[i] = 2.0*q_c[i*innumP + t] * tt[d];
					}
					//find b
					for (integer i = 0; i < indim; i++) {
						b[d*indim*innumP + i*innumP + t] = q_c[d*innumP + t] * q_c[i*innumP + t] / qnorm[t] + qnorm[t] * (i == d) - temp[i];
					}
				}
			}
			//********
			for (integer i = 0; i < indim; i++) {
				for (integer j = 0; j < indim; j++) {
					J[i*indim + j] = InnerProd_Q(b + i*innumP*indim, b + j*innumP*indim, innumP, indim);
				}
			}
			//beta = inv(J)*(-residual');
			//***************************************Check!!!**************************************
			//***************************************Check!!!**************************************
			//***************************************Check!!!**************************************
			//***************************************Check!!!**************************************
			dcopy_(&indim, tt, &GLOBAL::IONE, beta, &GLOBAL::IONE);
			coeff = -1.0;
			dscal_(&indim, &coeff, beta, &GLOBAL::IONE);
			dgesv_(&indim, &GLOBAL::IONE, J, &indim, IPIV, beta, &indim, &INFO);
			for (integer i = 0; i < innumP*indim; i++) {
				dq[i] = 0.0;
			}
			for (integer d = 0; d < indim; d++) {
				for (integer i = 0; i < indim; i++) {
					for (integer j = 0; j < innumP; j++) {
						dq[i*innumP + j] = dq[i*innumP + j] + beta[d] * b[d*innumP*indim + i*innumP + j];
					}
				}
			}
			dq_abs = std::sqrt(InnerProd_Q(dq, dq, innumP, indim));
			coeff1 = std::cos(dq_abs);
			coeff2 = std::sin(dq_abs) / dq_abs;
			dscal_(&NXD, &coeff1, q_c, &GLOBAL::IONE);
			daxpy_(&NXD, &coeff2, dq, &GLOBAL::IONE, q_c, &GLOBAL::IONE);

			for (integer i = 0; i < innumP; i++) {
				qnorm[i] = dnrm2_(&indim, q_c + i, &innumP);          ////check
			}
			for (integer i = 0; i < indim; i++) {
				tt[i] = InnerProd_Q(q_c + i*innumP, qnorm, innumP, 1);
			}

			iter++;
		}

		if (iter >= 100) {
			OUTSTREAM << "Item_1: Iterations exceeded 100" << std::endl;
		}

		delete[] qnorm;
		delete[] tt;
		delete[] temp;
		delete[] b;
		delete[] J;
		delete[] beta;
		delete[] dq;
		delete[] IPIV;
	}



	void PreShapePathStraighten::Item_2(const double *q, integer innumP, integer indim, double *w) //need duplicate w??
	{
		double temp0, coeff;  //integral
		double *temp = new double[indim];
		double *qnorm = new double[innumP];
		integer NXD = innumP*indim;
		//OUTSTREAM << "t1" << std::endl;//---
		double *w_c = new double[NXD];
		//OUTSTREAM << "t2" << std::endl;//---
		double *tt = new double[indim];
		//OUTSTREAM << "t3" << std::endl;//---
		double *b = new double[indim*innumP*indim];
		//OUTSTREAM << "t4" << std::endl;//---



		temp0 = InnerProd_Q(w, q, innumP, indim);
		//OUTSTREAM << "t5" << std::endl;//---
		temp0 = -temp0;
		//OUTSTREAM << "t6" << std::endl;//---
		dcopy_(&NXD, w, &GLOBAL::IONE, w_c, &GLOBAL::IONE);
		//OUTSTREAM << "t7" << std::endl;//---
		daxpy_(&NXD, &temp0, const_cast<double *> (q), &GLOBAL::IONE, w_c, &GLOBAL::IONE);
		//OUTSTREAM << "t8" << std::endl;//---

		for (integer i = 0; i < innumP; i++) {
			qnorm[i] = dnrm2_(&indim, const_cast<double *> (q)+i, &innumP);
		}
		//OUTSTREAM << "t9" << std::endl;//---

		//ForDebug::Print("qnorm:", qnorm, innumP);//---
		//ForDebug::Print("q:", q, innumP);//---
		//ForDebug::Print("q + innumP:", q + innumP, innumP);//---

		for (integer t = 0; t < indim; t++)
		{
			tt[t] = InnerProd_Q(q + t*innumP, qnorm, innumP, 1);
			//OUTSTREAM << "t:" << t << ", " << tt[t] << std::endl;//----

			//        double *PInnerProd = new double[innumP];
			//        ElasticCurvesRO::PointwiseInnerProd(q, qnorm, 1, innumP, PInnerProd);
			//        ForDebug::Print("Pinner", PInnerProd, innumP);
			//        double intv = 1.0/static_cast<double>((innumP-1));
			//        OUTSTREAM << "intv" <<intv<<std::endl;
			//        double check = ElasticCurvesRO::Trapz(PInnerProd, innumP, intv);
			//        OUTSTREAM <<"innumP" << innumP <<std::endl;
			//        OUTSTREAM <<"check"<< check <<std::endl;
		}
		//OUTSTREAM << "t10" << std::endl;//---

		//OUTSTREAM << tt[0]<<std::endl;
		//OUTSTREAM << tt[1] <<std::endl;

		for (integer t = 0; t < innumP; t++)
		{
			for (integer d = 0; d < indim; d++)
			{
				//find the d-th temp
				for (integer i = 0; i < indim; i++)
				{
					temp[i] = 2.0*q[i*innumP + t] * tt[d];
				}
				//find b
				for (integer i = 0; i < indim; i++)
				{
					b[i*indim*innumP + d*innumP + t] = q[d*innumP + t] * q[i*innumP + t] / qnorm[t] + qnorm[t] * (i == d) - temp[i];
				}
			}
		}
		//OUTSTREAM << "t11" << std::endl;//---


		//orthonormal basis of the normal space
		coeff = 1.0 / std::sqrt(InnerProd_Q(b, b, innumP, indim));
		//OUTSTREAM << "t12" << std::endl;//---

		dscal_(&NXD, &coeff, b, &GLOBAL::IONE);
		//OUTSTREAM << "t13" << std::endl;//---

		for (integer i = 1; i < indim; i++) {
			for (integer j = 0; j < i; j++) {
				coeff = -InnerProd_Q(b + i*innumP*indim, b + j*innumP*indim, innumP, indim);
				daxpy_(&NXD, &coeff, b + j*innumP*indim, &GLOBAL::IONE, b + i*innumP*indim, &GLOBAL::IONE);
			}
			coeff = 1.0 / std::sqrt(InnerProd_Q(b + i*innumP*indim, b + i*innumP*indim, innumP, indim));
			dscal_(&NXD, &coeff, b + i*innumP*indim, &GLOBAL::IONE);
		}
		//OUTSTREAM << "t14" << std::endl;//---


		//Project w into T_q(Cc)
		for (integer i = 0; i < indim; i++) {
			coeff = -InnerProd_Q(w_c, b + i*innumP*indim, innumP, indim);
			daxpy_(&NXD, &coeff, b + i*innumP*indim, &GLOBAL::IONE, w_c, &GLOBAL::IONE);
		}
		//OUTSTREAM << "t15" << std::endl;//---


		dcopy_(&NXD, w_c, &GLOBAL::IONE, w, &GLOBAL::IONE);
		//OUTSTREAM << "t16" << std::endl;//---

		delete[] temp;
		delete[] qnorm;
		delete[] w_c;
		delete[] tt;
		delete[] b;
	}


	//q1 q2 do not change
	void PreShapePathStraighten::Item_3(const double *w, const double *q1, const double *q2, integer innumP, integer indim, double *wbar)
	{
		double temp0, coeff, l;
		//OUTSTREAM << innumP << "," << indim << std::endl;//---
		double *q = new double[innumP*indim], *wtilde = new double[innumP*indim];
		//OUTSTREAM << "h1" << std::endl;//---
		double *qnorm = new double[innumP];
		//OUTSTREAM << "h2" << std::endl;//---
		double *tt = new double[indim];
		//OUTSTREAM << "h3" << std::endl;//---
		double *temp = new double[indim];
		//OUTSTREAM << "h4" << std::endl;//---
		double *b = new double[indim*innumP*indim];
		//OUTSTREAM << "h5" << std::endl;//---
		integer NXD = innumP*indim;
		//OUTSTREAM << "h6" << std::endl;//---

		temp0 = InnerProd_Q(w, q2, innumP, indim);
		for (integer i = 0; i < innumP*indim; i++) {
			q[i] = q1[i] + q2[i];
		}
		//OUTSTREAM << "h7" << std::endl;//---
		coeff = -2.0*temp0 / InnerProd_Q(q, q, innumP, indim);
		//OUTSTREAM << "h8" << std::endl;//---
		dcopy_(&NXD, const_cast<double *> (w), &GLOBAL::IONE, wtilde, &GLOBAL::IONE);
		//OUTSTREAM << "h9" << std::endl;//---
		daxpy_(&NXD, &coeff, q, &GLOBAL::IONE, wtilde, &GLOBAL::IONE);
		//OUTSTREAM << "h10" << std::endl;//---
		l = std::sqrt(InnerProd_Q(wtilde, wtilde, innumP, indim));   //need std??????
		//OUTSTREAM << "h11" << std::endl;//---

		for (integer i = 0; i < innumP; i++) {
			qnorm[i] = dnrm2_(&indim, const_cast<double *> (q2)+i, &innumP);     //check!!!!!!!!!!
		}
		for (integer t = 0; t < indim; t++) {
			tt[t] = InnerProd_Q(q2 + t*innumP, qnorm, innumP, 1);
		}

		for (integer t = 0; t < innumP; t++) {
			for (integer d = 0; d < indim; d++) {
				//find the d-th temp
				for (integer i = 0; i < indim; i++) {
					temp[i] = 2.0*q2[i*innumP + t] * tt[d];
				}
				//find b
				for (integer i = 0; i < indim; i++) {
					b[d*indim*innumP + i*innumP + t] = q2[d*innumP + t] * q2[i*innumP + t] / qnorm[t] + qnorm[t] * (i == d) - temp[i];
				}
			}
		}
		//Orthonormal basis of the normal space
		coeff = 1.0 / std::sqrt(InnerProd_Q(b, b, innumP, indim));
		dscal_(&NXD, &coeff, b, &GLOBAL::IONE);
		for (integer i = 1; i < indim; i++) {
			for (integer j = 0; j < i; j++) {
				coeff = -InnerProd_Q(b + i*innumP*indim, b + j*innumP*indim, innumP, indim);
				daxpy_(&NXD, &coeff, b + j*innumP*indim, &GLOBAL::IONE, b + i*innumP*indim, &GLOBAL::IONE);
			}
			coeff = 1.0 / std::sqrt(InnerProd_Q(b + i*innumP*indim, b + i*innumP*indim, innumP, indim));
			dscal_(&NXD, &coeff, b + i*innumP*indim, &GLOBAL::IONE);
		}

		//Project w into T_q(Cc)
		for (integer i = 0; i < indim; i++) {
			coeff = -InnerProd_Q(wtilde, b + i*innumP*indim, innumP, indim);
			daxpy_(&NXD, &coeff, b + i*innumP*indim, &GLOBAL::IONE, wtilde, &GLOBAL::IONE);
		}

		//Rescale
		if (std::sqrt(InnerProd_Q(wtilde, wtilde, innumP, indim)) > 1e-12) {
			coeff = l / std::sqrt(InnerProd_Q(wtilde, wtilde, innumP, indim));
			dscal_(&NXD, &coeff, wtilde, &GLOBAL::IONE);
		}

		dcopy_(&NXD, wtilde, &GLOBAL::IONE, wbar, &GLOBAL::IONE);    //if necessary???????

		delete[] q;
		delete[] wtilde;
		delete[] qnorm;
		delete[] tt;
		delete[] temp;
		delete[] b;
		//OUTSTREAM << "h6" << std::endl;//---
	}


	//Calculating inner prod of two points q1, q2 in shape space.
	double PreShapePathStraighten::InnerProd_Q(const double *q1, const double *q2, integer innumP, integer indim)
	{
		double intv, result;
		double *PInnerProd = new double[innumP];
		ElasticCurvesRO::PointwiseInnerProd(q1, q2, indim, innumP, PInnerProd);
		intv = 1.0 / static_cast<double>((innumP - 1));
		result = ElasticCurvesRO::Trapz(PInnerProd, innumP, intv);

		delete[] PInnerProd;
		return result;
	}
} /*end of ROPTLIB namespace*/
