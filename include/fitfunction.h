/*
 * futfunction.h
 *
 *  Created on: 29 Nov 2016
 *      Author: jkiesele
 */

#ifndef INCLUDE_FUTFUNCTION_H_
#define INCLUDE_FUTFUNCTION_H_


#include "Math/IFunction.h"
#include <stdexcept>
#include "combiner.h"
#include <iostream>

/*
 * Interface for function to be minimised, gradients etc. detached from combiner for better overview
 */

class fitfunctionBase{
public:
	fitfunctionBase():c(0){}
	fitfunctionBase(const combiner* c_in):c(c_in){}
	fitfunctionBase(const fitfunctionBase& r):c(r.c){}
	fitfunctionBase& operator = (const fitfunctionBase& r){	c=r.c;return *this;}
protected:
	void eval(const double*x, double& f, double* df)const;
	inline int nDim()const{
		if(!c) throw std::logic_error("fitfunctionBase::nDim no combiner associated");
		return c->npars_;
	}
private:
	const combiner * c;

};

#define INHERITFROM ROOT::Math::IBaseFunctionMultiDim //
//#define INHERITFROM ROOT::Math::IGradientFunctionMultiDim
//#define INHERITFROM  ROOT::Minuit2::FCNGradientBase

class fitfunctionGradient:  public fitfunctionBase,
public INHERITFROM {

public:
	fitfunctionGradient():fitfunctionBase(),INHERITFROM(){}
	fitfunctionGradient(const combiner* c_in):fitfunctionBase(c_in),INHERITFROM(){}
	fitfunctionGradient(const fitfunctionGradient& r):fitfunctionBase(r){}
	fitfunctionGradient& operator = (const fitfunctionGradient& r) {fitfunctionBase::operator =(r); return *this;}
	unsigned int	NDim() const{
		return fitfunctionBase::nDim();
	}
	INHERITFROM*	Clone() const{
		return new fitfunctionGradient(*this);
	}
	void Gradient(const double* x, double* grad) const{
		double f=0;
		FdF(x,f,grad);
	}
	void FdF(const double* x, double& f, double* df)const{
		fitfunctionBase::eval(x,f,df);
	}
/*
	std::vector<double> Gradient(const std::vector<double>& x) const{
		std::vector<double> df(NDim(),0);
		double f;
		eval(&x.at(0),f,&df.at(0));
		return df;
	}

	double operator()(const std::vector<double>& x) const{
		std::vector<double> df(NDim(),0);
		double f;
		eval(&x.at(0),f,&df.at(0));
		return f;
	}

	double Up() const{return 1;}
*/
private:
	double	DoDerivative(const double* x, unsigned int icoord) const{return 0;}
	double	DoEval(const double* x) const{
		double f=0;
		static std::vector<double> df(NDim());
		eval(x,f,&df.at(0));
		return f;
	}

};

#endif /* INCLUDE_FUTFUNCTION_H_ */
