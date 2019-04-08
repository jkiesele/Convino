/*
 * fitfunction.cpp
 *
 *  Created on: 29 Nov 2016
 *      Author: jkiesele
 */




#include "fitfunction.h"
#include "combiner.h"
#include <stdexcept>

//MP support here is still experimental
#ifdef USE_MP
#include <omp.h>
#endif



void fitfunctionBase::eval(const double*x, double& f, double* df)const{
	if(!c) throw std::logic_error("fitfunctionBase::eval no combiner associated");

	f=0;
	long double f_int=0;

	const std::vector<measurement>& meass=c->measurements_;
	size_t nmeas=meass.size();

//#ifdef USE_MP
//	size_t maxthreads=omp_get_max_threads();
//	omp_set_num_threads(std::min(nmeas,maxthreads));
//#pragma omp parallel for reduction(+:f)
//#endif
	for(size_t i=0;i<nmeas;i++){
		double f_local=meass.at(i).evaluate(x,df,c->lh_mod_==combiner::lh_mod_pearson,nDim());
		f_int+=(long double)f_local;
	}
	if(meass.size()<1)
	    throw std::out_of_range("fitfunctionBase::eval: No measurements associated");


	//this part needs more care, make it  iterative  use augmented Lagrangian method
	//mu and lambda parameters will be part of combiner.
	//if(meass.at(0).getExcludeBin()>=0){
	//    double combsumdiff = 1. - meass.at(0).getCombSum(x);
	//    f_int+= c->auglagrangemu_/2. * combsumdiff * combsumdiff  - c->auglagrangelambda_ * combsumdiff ;
	//}
	const size_t nsys=c->npars_-c->nest_;

	const double sqrt2 = std::sqrt(2);

	if(f_int!=f_int){
		throw std::runtime_error("fitfunctionBase::eval: nan in measurement chi2");
	}


	for(size_t i=0;i<nsys;i++){
		double pi=x[i];

		if(c->allsysparas_.at(i).getType() == parameter::para_unc_lognormal){//this is switched off for now
			if(x[i]>-1)
				pi = std::log(1 + x[i]); //this needs some re-work and transofrmation before, TBI
			else
				pi = 1e4*x[i];
		}
		for(size_t j=i;j<nsys;j++){

			double pj=x[j];
			if(c->allsysparas_.at(j).getType() == parameter::para_unc_lognormal){
				if(x[j]>-1)
					pj = sqrt2* std::log(1 + x[j]);
				else
					pj = 1e4*x[j];
			}
			if(i==j)
			    f_int += (long double) (pi * c->inv_priors_[i][j] * pj);
			else
			    f_int += (long double) (2.* pi * c->inv_priors_[i][j] * pj);
		}
	}

	if(f_int!=f_int){
		for(int i=0;i<nDim();i++){
			std::cout << x[i] << std::endl;
		}
		throw std::runtime_error("fitfunctionBase::eval: nan in external correlations");
	}
	f=f_int;


}


