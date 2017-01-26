/*
 * simpleFitter.h
 *
 *  Created on: May 22, 2014
 *      Author: kiesej
 */

#ifndef SIMPLEFITTER_H_
#define SIMPLEFITTER_H_

#include <TROOT.h>
#include <vector>
#include "Math/Functor.h"
#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include <algorithm>

#include "triangularMatrix.h"




/**
 * Wrapper around TMinuit and other parallelization tools
 * In principle allows for some multithreading
 */
class simpleFitter{
public:
	enum printlevels{pl_silent,pl_normal,pl_verb}; //TBI
	enum mintypes{type_minuit, type_gsl_nonlinear, type_gsl,type_gsl_simanealing};

	void setType(mintypes t) {type_=t;}

	simpleFitter();
	~simpleFitter();
	simpleFitter(const simpleFitter&);
	simpleFitter& operator=(const simpleFitter&);
	/**
	 * Available fit modes:
	 * <some docu>
	 */
	void setMinFunction( ROOT::Math::Functor* f);
	void setMinFunction(const ROOT::Math::IGradientFunctionMultiDim &f);
	void setMinFunction(const ROOT::Math::IBaseFunctionMultiDim& f);

	/**
	 * these are then considered as start values
	 * if not defined, will assume 0 for all parameters
	 * and step sizes of 0.01 for all parameters
	 */
	void setParameters(const std::vector<double>& inpars,const std::vector<double>& steps=std::vector<double>());
	void setParameterNames(const std::vector<TString>& nms){paranames_=nms;}

	void setParameter(size_t idx,double value);
	void setParameterFixed(size_t idx,bool fixed=true);
	void setParameterLowerLimit(size_t idx,double value);
	void setParameterUpperLimit(size_t idx,double value);
	void removeParameterLowerLimit(size_t idx);
	void removeParameterUpperLimit(size_t idx);


	void setMaxCalls(unsigned  int calls){maxcalls_=calls;}

	void setAsMinosParameter(size_t parnumber,bool set=true);
	void setTolerance(double tol){tolerance_=tol;}


	static int printlevel;

	void setStrategy(int str){strategy_=str;}
	void feedErrorsToSteps();



	void fit();

	const std::vector<double> *getParameters()const{return &paras_;}
	/**
	 * these are signed!
	 */
	const std::vector<double> *getParameterErrUp()const{return &paraerrsup_;}
	/**
	 * these are signed!
	 */
	const std::vector<double> *getParameterErrDown()const{return &paraerrsdown_;}
	/**
	 * always > 0
	 */
	double getParameterErr(size_t idx)const;
	const double&  getParameter(size_t idx)const;
	const std::vector<TString> *getParameterNames()const {return &paranames_;}


	const double& getCorrelationCoefficient(size_t i, size_t j)const;
	const double& getHessianCoefficient(size_t i, size_t j)const;

	correlationMatrix getCorrelationMatrix()const;
	triangularMatrix getHessianMatrix()const;

	void getParameterErrorContribution(size_t a, size_t b,double & errup, double& errdown);
	/**
	 * gets the contribution of a to b by fixing the parameters a, repeating minos and returning the changes in errup and errdown
	 * failures are indicated by negative return values!
	 * requires parameter b to be a minos parameter
	 */
	void getParameterErrorContributions(std::vector<size_t> a, size_t b,double & errup, double& errdown);

	/**
	 * Fixes all parameters and thereby determines the statistical component of the
	 * error to parameter b only
	 */
	void getStatErrorContribution( size_t b,double & errup, double& errdown);

	size_t findParameterIdx(const TString& paraname)const;


	bool wasSuccess(){return minsuccessful_;}
	bool minosWasSuccess(){return minossuccessful_;}


	void setOnlyRunDummy(bool dummyrun){dummyrun_=dummyrun;}


	static bool debug;

	const double& getChi2Min()const{return chi2min_;}


private: //set some to protected if inheritance is needed
	double chi2min_;

	std::vector<double> paras_;
	std::vector<double> stepsizes_;
	std::vector<double> paraerrsup_,paraerrsdown_;
	std::vector<TString> paranames_;
	std::vector<bool> parafixed_;


	std::vector< std::vector<double> > paracorrs_;
	std::vector< std::vector<double> > hessian_;

	std::vector<int> minospars_;
	std::vector<std::pair<bool, double> > lowerlimits_;
	std::vector<std::pair<bool, double> > upperlimits_;

	unsigned int  maxcalls_;

	//definitely privateL

	//1D section
	//this is moved to global due to tminuit reasons
	// void chisq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


	bool checkSizes()const;

	bool minsuccessful_,minossuccessful_;

	double tolerance_;
	const ROOT::Math::Functor* functobemin_;
	const ROOT::Math::IGradientFunctionMultiDim* gradfunctobemin_;
	const ROOT::Math::IBaseFunctionMultiDim * basefunctobemin_;

	ROOT::Math::Minimizer * pminimizer_;

	void copyFrom(const simpleFitter&);

	ROOT::Math::Minimizer * invokeMinimizer()const;


	int strategy_;
	bool dummyrun_;

	mintypes type_;
};




#endif /* SIMPLEFITTER_H_ */
