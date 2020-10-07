/*
 * simpleFitter.cc
 *
 *  Created on: May 22, 2014
 *      Author: kiesej
 */

#include "simpleFitter.h"
#include <TROOT.h>
#include <TMinuit.h>
#include "Minuit2/MnStrategy.h"
#include <iostream>
#include "Math/Functor.h"
#include "Math/Factory.h"
#include <Minuit2/MnMachinePrecision.h>
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include <stdexcept>
#include "Minuit2/Minuit2Minimizer.h"
#include "TMinuitMinimizer.h"
//#include "Math/GSLMinimizer.h"
//#include "Math/GSLNLSMinimizer.h"
//#include "Math/GSLSimAnMinimizer.h"
#include "TMinuitMinimizer.h"
#include <limits>
#include <cfloat>
//use TMinuit
//define all the functions - nee dto be c-stylish








int simpleFitter::printlevel=0;
bool simpleFitter::debug=false;

simpleFitter::simpleFitter():chi2min_(1e6),maxcalls_(4e8),
		minsuccessful_(false),minossuccessful_(false),tolerance_(0.1),functobemin_(0),
		gradfunctobemin_(0),basefunctobemin_(0),pminimizer_(0),
		strategy_(2),dummyrun_(false),fastmode_(false),type_(type_minuit){
    ncontourpoints_=100;
}

simpleFitter::~simpleFitter(){
	if(pminimizer_){
		delete pminimizer_;
		pminimizer_=0;
	}
}

simpleFitter::simpleFitter(const simpleFitter& rhs){
	copyFrom(rhs);
}
simpleFitter& simpleFitter::operator=(const simpleFitter& rhs){
	copyFrom(rhs);
	return *this;
}

void simpleFitter::setMinFunction( ROOT::Math::Functor* f){
	gradfunctobemin_=0;
	basefunctobemin_=0;
	functobemin_= f;
}

void simpleFitter::setMinFunction(const ROOT::Math::IGradientFunctionMultiDim &f){
	functobemin_=0;
	basefunctobemin_=0;
	gradfunctobemin_=&f;
}
void simpleFitter::setMinFunction(const ROOT::Math::IBaseFunctionMultiDim& f){
	gradfunctobemin_=0;
	functobemin_=0;
	basefunctobemin_=&f;
}


void simpleFitter::setParameters(const std::vector<double>& inpars,const std::vector<double>& steps){
	paras_=inpars;
	if(steps.size()>0)
		stepsizes_=steps;
	else
		stepsizes_.resize(paras_.size(),1);
	paraerrsup_.resize(paras_.size(),0);
	paraerrsdown_.resize(paras_.size(),0);
	parafixed_.resize(paras_.size(),false);
	std::vector<double> dummy;
	dummy.resize(paras_.size(),0);
	paracorrs_.clear();
	paracorrs_.resize(paras_.size(),dummy);
	std::pair<bool, double> limits(false,0);
	lowerlimits_.resize(paras_.size(),limits);
	upperlimits_.resize(paras_.size(),limits);

	makecontours_.resize(paras_.size());

}

void simpleFitter::setParameter(size_t idx,double value){
	if(idx >= paras_.size())
		throw std::out_of_range("simpleFitter::setParameter: index out of range");
	paras_.at(idx)=value;
}
void simpleFitter::setParameterFixed(size_t idx,bool fixed){
	if(idx >= paras_.size())
		throw std::out_of_range("simpleFitter::setParameterFixed: index out of range");
	parafixed_.at(idx)=fixed;
}


void simpleFitter::setParameterFixed(const TString& paraname,bool fixed){
    setParameterFixed(findParameterIdx(paraname),fixed);
}

void simpleFitter::setParameterLowerLimit(size_t idx,double value){
	if(idx >= paras_.size())
		throw std::out_of_range("simpleFitter::setParameterLowerLimit: index out of range");
	std::pair<bool, double> limit(true,value);
	lowerlimits_.at(idx)=limit;
}
void simpleFitter::setParameterUpperLimit(size_t idx,double value){
	if(idx >= paras_.size())
		throw std::out_of_range("simpleFitter::setParameterUpperLimit: index out of range");
	std::pair<bool, double> limit(true,value);
	upperlimits_.at(idx)=limit;
}
void simpleFitter::removeParameterLowerLimit(size_t idx){
	if(idx >= paras_.size())
		throw std::out_of_range("simpleFitter::removeParameterLowerLimit: index out of range");
	std::pair<bool, double> limit(false,0);
	lowerlimits_.at(idx)=limit;
}
void simpleFitter::removeParameterUpperLimit(size_t idx){
	if(idx >= paras_.size())
		throw std::out_of_range("simpleFitter::removeParameterLowerLimit: index out of range");
	std::pair<bool, double> limit(false,0);
	upperlimits_.at(idx)=limit;
}

void simpleFitter::setAsMinosParameter(size_t parnumber,bool set){
	std::vector<int>::iterator it=std::find(minospars_.begin(),minospars_.end(),parnumber);
	size_t pos=it-minospars_.begin();
	if(set){
		if(pos >= minospars_.size())
			minospars_.push_back(parnumber);
	}
	else{
		if(minospars_.size()>0 && it!=minospars_.end()){
			minospars_.erase(it);
		}
	}
}


void simpleFitter::feedErrorsToSteps(){
	for(size_t i=0;i<stepsizes_.size();i++)
		stepsizes_.at(i)=paraerrsup_.at(i)/100;//*0.5; //good enough
}

void simpleFitter::fit(){

	if(!checkSizes()){
		std::cout << "ERROR IN simpleFitter::fit" <<std::endl;

		std::cout << "parameter:             " << paras_.size() << std::endl;
		std::cout << "parameter errors up:   " << paraerrsup_.size() << std::endl;
		std::cout << "parameter errors down: " << paraerrsdown_.size() << std::endl;
		std::cout << "parameter steps size:  " << stepsizes_.size() << std::endl;

		throw std::runtime_error("simpleFitter::fit: number of parameters<->stepsizes or parameters<->fitfunction ");
		//  return;
	}

	if(printlevel>0)
		std::cout << "simpleFitter::fit(): fitting "  << paras_.size() <<  " parameters." <<std::endl;

	minsuccessful_=false;

	// invoke min, set func, fill errors, chi2min, correlations
	if(dummyrun_){
		std::cout << "simpleFitter::fit: warning. Dummyrun mode -> will just fill dummies."<<std::endl;
		chi2min_=0;
		if(pminimizer_){
			delete pminimizer_;
			pminimizer_=0;
		}
		pminimizer_=invokeMinimizer();
		std::vector<double> dummy;
		dummy.resize(paras_.size(),0);
		paracorrs_.clear();
		paracorrs_.resize(paras_.size(),dummy);
		for(size_t i=0;i<paras_.size();i++){
			paraerrsdown_.at(i)=-1;
			paraerrsup_.at(i)=1;
			for(size_t j=0;j<paras_.size();j++){
				if(i!=j)
					paracorrs_.at(i).at(j)=0;
				else
					paracorrs_.at(i).at(j)=1;
			}
		}
		minossuccessful_=true;
		minsuccessful_=true;
		return;
	}


	//  if(functobemin_)
	if(pminimizer_){
		delete pminimizer_;
		pminimizer_=0;
	}

	pminimizer_=invokeMinimizer();
	if(printlevel>0)
		std::cout << "simpleFitter::fit(): minimizer prepared "<<std::endl;



	minsuccessful_=pminimizer_->Minimize();
	if(printlevel>0)
		std::cout << "simpleFitter::fit(): minimization done" <<std::endl;
	const double *xs = pminimizer_->X();
	const double * errs=pminimizer_->Errors();
	//feed back
	if(minospars_.size()>0)
		minossuccessful_=true;
	for(size_t i=0;i<paras_.size();i++){
		paras_.at(i)=xs[i];
		TString paraname="";
		if(paranames_.size()>i)
			paraname=paranames_.at(i);
		else
			paraname=(TString)"par"+(i);
		if(std::find(minospars_.begin(),minospars_.end(),i) != minospars_.end()){

			if(errs && !pminimizer_->GetMinosError(i, paraerrsdown_.at(i), paraerrsup_.at(i))){
				paraerrsdown_.at(i)=-errs[i];
				paraerrsup_.at(i)=errs[i];
				minossuccessful_=false;
			}
		}
		else if(errs){
			paraerrsdown_.at(i)=-errs[i];
			paraerrsup_.at(i)=errs[i];
		}
		if(printlevel>0){
			TString someblanks;
			if(paraname.Length()<20){
				size_t add=20-paraname.Length();
				for(size_t k=0;k<add;k++)
					someblanks+=" ";
			}
			if(paras_.at(i)>0)
				someblanks+=" ";
			if(printlevel>-1)
				std::cout << paraname+" "  << someblanks << paras_.at(i)<< " +"<< paraerrsup_.at(i) << " " << paraerrsdown_.at(i)<< std::endl;
		}
		try{
		    for(const auto& pb: makecontours_.at(i)){

		        TString paraname_b="";
		        if(paranames_.size()>pb)
		            paraname_b=paranames_.at(pb);
		        else
		            paraname_b=(TString)"par"+(pb);
		        if(printlevel>0 || true)
		            std::cout << "calculating contour for "<< paraname << " vs " <<paraname_b<<std::endl;

		        std::vector<double> x(ncontourpoints_);
		        std::vector<double> y(ncontourpoints_);
		        bool csucc = pminimizer_->Contour(i,pb,ncontourpoints_,&x.at(0),&y.at(0));

		        if(!csucc){
		            std::cout << "calculating contour for "<< paraname << " vs " <<paraname_b<< "failed" << std::endl;
		        }
		        else
		            contours_.push_back({{paraname,paraname_b}, {x,y}});


		    }
		}catch(std::exception& ex){
		    std::cout <<"exception in contour: "<<  ex.what() << ", makecontours_.size() "
		            << makecontours_.size() << ", i: "<< i << std::endl;
		    throw ex;
		}

		if(printlevel>0)
		    std::cout << paraname << std::endl;
	}


    if(printlevel>0){
        std::cout << "simpleFitter::fit(): minimisation done, errors done, filling info" << std::endl;}

	chi2min_=pminimizer_->MinValue();
	std::vector<double> dummy;
	dummy.resize(paras_.size(),0);
	paracorrs_.clear();
	paracorrs_.resize(paras_.size(),dummy);
	for(size_t i=0;i<paras_.size();i++){
		for(size_t j=0;j<paras_.size();j++){
			paracorrs_.at(i).at(j)=pminimizer_->Correlation(i,j);
		}
	}

	hessian_.resize(paras_.size(),std::vector<double> (paras_.size(),0));
	if(!fastmode_){
	    if(printlevel>0)
	        std::cout << "simpleFitter::fit(): getting Hessian" << std::endl;
	    double hessian_s[paras_.size()*paras_.size()];
	    pminimizer_->GetHessianMatrix(hessian_s);

	    for(size_t i=0;i<paras_.size();i++){
	        for(size_t j=0;j<paras_.size();j++){
	            hessian_.at(i).at(j)=hessian_s[i *paras_.size() + j];
	        }
	    }
	}

	if(printlevel>1){
		for(size_t i=0;i<paras_.size();i++){
			std::cout << paranames_.at(i) << "\t";
			for(size_t j=0;j<paras_.size();j++){
				std::cout << getCorrelationCoefficient(i,j) << " ";
			}
			std::cout << std::endl;
		}
	}

    if(printlevel>0)
        std::cout << "simpleFitter::fit(): done" << std::endl;

}

const double& simpleFitter::getParameter(size_t idx)const{
	if(idx>=paras_.size())
		throw std::out_of_range("simpleFitter::getParameter: index out of range");

	return paras_.at(idx);
}

double simpleFitter::getParameterErr(size_t idx)const{
	if(idx>=paraerrsup_.size())
		throw std::out_of_range("simpleFitter::getParameterErr: index out of range");
	if(paraerrsdown_.size()!=paraerrsup_.size())
		throw std::out_of_range("simpleFitter::getParameterErr: serious problem errors up/down not of same size");

	double err=fabs(paraerrsup_.at(idx));
	if(err<fabs(paraerrsdown_.at(idx)))err=fabs(paraerrsdown_.at(idx));
	return err;
}


const double& simpleFitter::getCorrelationCoefficient(size_t i, size_t j)const{
	if(i>=paracorrs_.size())
		throw std::out_of_range("simpleFitter::getCorrelationCoefficient: first index out of range");
	if(j>=paracorrs_.at(i).size())
		throw std::out_of_range("simpleFitter::getCorrelationCoefficient: second index out of range");

	return paracorrs_.at(i).at(j);

}

const double& simpleFitter::getHessianCoefficient(size_t i, size_t j)const{
	if(i>=hessian_.size())
		throw std::out_of_range("simpleFitter::getHessianCoefficient: first index out of range");
	if(j>=hessian_.at(i).size())
		throw std::out_of_range("simpleFitter::getHessianCoefficient: second index out of range");
	return hessian_.at(i).at(j);
}

correlationMatrix simpleFitter::getCorrelationMatrix()const{
	correlationMatrix m;
	if(paranames_.size())
		m=correlationMatrix(paranames_);
	else
		m=correlationMatrix(paras_.size());
	for(size_t i=0;i<m.size();i++){
		for(size_t j=i;j<m.size();j++){
			m.setEntry(i,j,paracorrs_.at(i).at(j));
		}
	}
	return m;
}

triangularMatrix simpleFitter::getCovarianceMatrix()const{
    auto hess = getHessianMatrix();
    if(!fastmode_)
        hess.invert();
    return hess;
}

triangularMatrix simpleFitter::getHessianMatrix()const{
	triangularMatrix m;
	if(paranames_.size())
		m=triangularMatrix(paranames_);
	else
		m=triangularMatrix(paras_.size());
	for(size_t i=0;i<m.size();i++){
		for(size_t j=i;j<m.size();j++){
			m.setEntry(i,j,hessian_.at(i).at(j));
		}
	}
	return m;
}



/**
 * always returns positive values if no fault
 */
void simpleFitter::getParameterErrorContribution(size_t a, size_t b,double & errup, double& errdown){
	if(!wasSuccess())
		throw std::logic_error("simpleFitter::getParameterErrorContribution: can only be called after fit");
	if(!pminimizer_)
		throw std::logic_error("simpleFitter::getParameterErrorContribution: can only be called after fit (minimizer pointer 0) IF YOU SEE THIS< SOMETHING IS SERIOUSLY WRONG");
	if(a>=paras_.size() || b>=paras_.size())
		throw std::out_of_range("simpleFitter::getParameterErrorContribution: one index out of range");
	if(parafixed_.at(a)){
		errdown=0;
		errup=0;
		return;
	}

	std::vector<size_t> as;
	as.push_back(a);
	if(dummyrun_){
		errup=1;
		errdown=-1;
	}
	else
		getParameterErrorContributions(as,b,errup,errdown);

}

void simpleFitter::getParameterErrorContributions(std::vector<size_t> a, size_t b,double & errup, double& errdown){

	if(!wasSuccess())
		throw std::logic_error("simpleFitter::getParameterErrorContributions: can only be called after fit");
	if(!pminimizer_)
		throw std::logic_error("simpleFitter::getParameterErrorContributions: can only be called after fit (minimizer pointer 0) IF YOU SEE THIS< SOMETHING IS SERIOUSLY WRONG");
	if(a.size()>=paras_.size() || b>=paras_.size())
		throw std::out_of_range("simpleFitter::getParameterErrorContributions: one index out of range");
	if(a.size()<1)
		throw std::out_of_range("simpleFitter::getParameterErrorContributions: input indices empty");
	//work on hesse errors
	if(dummyrun_){
		errup=100000;
		errdown=-100000;
		return;
	}


	double olderr=pminimizer_->Errors()[b];
	//std::cout << "olderr: " << olderr<<std::endl;
	ROOT::Math::Minimizer *  min=invokeMinimizer();
	for(size_t i=0;i<a.size();i++)
		min->SetFixedVariable(a.at(i),paranames_.at(a.at(i)).Data(),paras_.at(a.at(i)));


	bool succ=min->Minimize();
	if(succ){
		if(std::find(minospars_.begin(),minospars_.end(),b) != minospars_.end()){// minos para
			double merrd,merrup;
			double oldup=paraerrsup_.at(b);
			double olddown=paraerrsdown_.at(b);
			bool succ = min->GetMinosError(b, merrd,merrup);
			if(succ){
				if(olddown*olddown-merrd*merrd > 0 )
					errdown=sqrt(olddown*olddown-merrd*merrd);
				else if(fabs(olddown / merrd -1 ) < 0.01) //prob only numerical
					errdown=0;
				else
					errdown=2e12;//produce a nan

				if(oldup*oldup - merrup*merrup >0)
					errup=sqrt(oldup*oldup - merrup*merrup);
				else if(fabs(oldup / merrup -1) < 0.01) //prob only numerical
					errup=0;
				else
					errup=2e12;//produce a nan

				delete min;
				return;
			} //if not make symm
			std::cout << "simpleFitter::getParameterErrorContributions: warning. evaluating minos error failed. will revert to symmetric error"
					<<std::endl;
		}

		errdown=min->Errors()[b];
		//std::cout << "newerr: " << errdown<<std::endl;
		double reldiff= (olderr-errdown)/olderr;
		if(reldiff>=0){
			errdown=sqrt(olderr*olderr-errdown*errdown);}
		else{
			if(reldiff>-0.001) //purely numerical difference
				errdown=0;
			else
				errdown=-1000;
		}
		errup=errdown;
		if(a.size()<2){
			if(getCorrelationCoefficient(a.at(0),b)<0){
				errup=-errdown;
			}
			else{
				errdown=-errdown;
			}
		}
	}
	else{
		errdown=-200000;
		errup=errdown;
	}
	delete min;
}
void simpleFitter::getStatErrorContribution(size_t b,double & errup, double& errdown){

	if(!wasSuccess())
		throw std::logic_error("simpleFitter::getStatErrorContribution: can only be called after fit");
	if(!pminimizer_)
		throw std::logic_error("simpleFitter::getStatErrorContribution: can only be called after fit (minimizer pointer 0) IF YOU SEE THIS< SOMETHING IS SERIOUSLY WRONG");
	if(b>=paras_.size())
		throw std::out_of_range("simpleFitter::getStatErrorContribution:  index out of range");

	if(dummyrun_){
		errup=1000000;
		errdown=1000000;
		return;
	}

	//work on hesse errors
	double olderr=pminimizer_->Errors()[b];
	//std::cout << "olderr: " << olderr<<std::endl;
	ROOT::Math::Minimizer *  min=invokeMinimizer();
	for(size_t i=0;i<paranames_.size();i++){
		if(i==b) continue;
		min->SetFixedVariable(i,paranames_.at(i).Data(),paras_.at(i));
	}


	bool succ=min->Minimize();
	if(succ){
		errdown=min->Errors()[b];
		double reldiff= (olderr-errdown)/olderr;
		if(reldiff<0){
			if(reldiff>-0.001) //purely numerical difference
				errdown=0;
			else
				errdown=-100000;
		}

	}
	else{
		errdown=-200000;
	}
	errup=errdown;

	delete min;
}



size_t simpleFitter::findParameterIdx(const TString& paraname)const{
	size_t idx=std::find(paranames_.begin(), paranames_.end(),paraname) - paranames_.begin();
	if(idx<paranames_.size()) return idx;
	throw std::runtime_error("simpleFitter::findParameterIdx: name not found");
}




bool simpleFitter::checkSizes()const{


	if(paras_.size() != stepsizes_.size())
		return false;
	if(paras_.size() != paraerrsup_.size())
		return false;
	if(paras_.size() != paraerrsdown_.size())
		return false;
	if(paras_.size() != lowerlimits_.size())
		return false;
	if(paras_.size() != upperlimits_.size())
		return false;

	return true;
}



void simpleFitter::copyFrom(const simpleFitter& rhs){
	if(this==&rhs) return;

	chi2min_=rhs.chi2min_;

	paras_=rhs.paras_;
	stepsizes_=rhs.stepsizes_;
	paraerrsup_=rhs.paraerrsup_;
	paraerrsdown_=rhs.paraerrsdown_;
	paranames_=rhs.paranames_;
	parafixed_=rhs.parafixed_;


	paracorrs_=rhs.paracorrs_;
	hessian_=rhs.hessian_;

	minospars_=rhs.minospars_;
	lowerlimits_=rhs.lowerlimits_;
	upperlimits_=rhs.upperlimits_;

	maxcalls_=rhs.maxcalls_;


	minsuccessful_=rhs.minsuccessful_;
	minossuccessful_=rhs.minossuccessful_;

	tolerance_=rhs.tolerance_;
	functobemin_=rhs.functobemin_;
	gradfunctobemin_=rhs.gradfunctobemin_;
	basefunctobemin_=rhs.basefunctobemin_;

	pminimizer_=0;//gets recreated

	strategy_=rhs.strategy_;
	dummyrun_=rhs.dummyrun_;
	fastmode_=rhs.fastmode_;

	type_=rhs.type_;

	makecontours_=rhs.makecontours_;
	contours_ = rhs.contours_;

}

ROOT::Math::Minimizer * simpleFitter::invokeMinimizer()const{
	if(printlevel>0)
		std::cout << "simpleFitter::invokeMinimizer" <<std::endl;
	ROOT::Math::Minimizer* min=0;
	if(type_==type_minuit){
		//min=new TMinuitMinimizer();
		min=new ROOT::Minuit2::Minuit2Minimizer("");
		if(printlevel>0)
			std::cout << "simpleFitter::invokeMinimizer: got Minuit2Minimizer" <<std::endl;
	}
	//else if(type_==type_gsl){
	//	min=new ROOT::Math::GSLMinimizer();
	//	if(printlevel>0)
	//		std::cout << "simpleFitter::invokeMinimizer: got GSLMinimizer" <<std::endl;
	//}
	//else if(type_==type_gsl_nonlinear){
	//	min=new ROOT::Math::GSLNLSMinimizer();
	//	if(printlevel>0)
	//		std::cout << "simpleFitter::invokeMinimizer: got GSLNLSMinimizer" <<std::endl;
	//}
	//else if(type_==type_gsl_simanealing){
	//	min=new ROOT::Math::GSLSimAnMinimizer();
	//	if(printlevel>0)
	//		std::cout << "simpleFitter::invokeMinimizer: got GSLNLSMinimizer" <<std::endl;
	//}


	if(!min)
		throw std::runtime_error("simpleFitter::invokeMinimizer");
	if(printlevel>0){
		min->SetPrintLevel(printlevel);
	}
	else
		min->SetPrintLevel(0);

	if(functobemin_){
		min->SetFunction(*functobemin_);
	}
	else if(basefunctobemin_){
		min->SetFunction(*basefunctobemin_);
	}
	else if(gradfunctobemin_){
		min->SetFunction(*gradfunctobemin_);
	}
	else{
		throw std::logic_error("simpleFitter::invokeMinimizer: no function to be minimized");
	}

	min->SetMaxFunctionCalls(maxcalls_); // for Minuit/Minuit2
	min->SetMaxIterations(maxcalls_);  // for GSL
	min->SetTolerance(tolerance_);



	for(size_t i=0;i<paras_.size();i++){
		TString paraname="";
		if(paranames_.size()>i)
			paraname=paranames_.at(i);
		else
			paraname=(TString)"par"+(i);


		if(printlevel>1)
			std::cout << "creating parameter " << i<<"/"<<paras_.size() << ": " << paraname <<std::endl;
		if(parafixed_.at(i))
			min->SetFixedVariable((unsigned int)i,paraname.Data(),paras_.at(i));
		else
			min->SetVariable((unsigned int)i,paraname.Data(),paras_.at(i), stepsizes_.at(i));

		if(lowerlimits_.at(i).first && upperlimits_.at(i).first){
			if(lowerlimits_.at(i).second>=upperlimits_.at(i).second){
				std::string errstr="simpleFitter::invokeMinimizer: parameter ";
				errstr+= paraname.Data();
				errstr+=" has upper limit < lower limit";
				throw std::out_of_range(errstr);
			}
			min->SetLimitedVariable((unsigned int)i,paraname.Data(),paras_.at(i), stepsizes_.at(i),
					lowerlimits_.at(i).second,upperlimits_.at(i).second);
		}
		else if(lowerlimits_.at(i).first){
			min->SetLowerLimitedVariable((unsigned int)i,paraname.Data(),paras_.at(i), stepsizes_.at(i),
					lowerlimits_.at(i).second);
		}
		else if(upperlimits_.at(i).first){
			min->SetUpperLimitedVariable((unsigned int)i,paraname.Data(),paras_.at(i), stepsizes_.at(i),
					upperlimits_.at(i).second);
		}
	}

	min->SetPrecision(-1);//DBL_EPSILON); //16
	min->SetStrategy(strategy_);




	return min;
}




