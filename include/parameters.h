/*
 * parameters.h
 *
 *  Created on: 4 Jul 2016
 *      Author: jkiesele
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "TString.h"
#include <stdexcept>


/*
 * This is a mix of uncertainty and parameter here
 * the uncertainty class has asymmetric info, the parameter class
 * has the way to add uncertainties on those uncertainties... maybe that should be changed?
 *
 * Maybe add this to the uncertainty parameter class? ... but then the uncertainty parameter class
 * needs parameter associations?... hmm...
 *
 * The uncertainty class is already per bin, per systematic.
 * maybe this could just point to a free parameter (or a shared one).. that's probably better than blowing up
 * the parameter class
 *
 * uncertainty now inherits from parameter.
 * actually, Lk_ then describes a full uncertainty covariance in (est,unc) with
 * parameters loose coupled to the underlying source uncertainty - or fully coupled, that's the
 * technical challenge (parameter associations).
 * This part cannot cover contributions across measurements, which are only incorporated through the
 * sources of uncertainties.
 * One COULD add additional couplings from an off-diagonal block... could be part of
 * combiner and would NOT require additional parameters, only penalty terms and an interface..
 * the interface could use the naming schemes
 *
 */


class parameter{
public:

	enum para_type{para_unc_absolute,para_unc_relative,para_unc_lognormal,para_estimate};

	parameter():startval_(0),type_(para_unc_absolute),associatedto_(0),stat_(1){}
	virtual ~parameter(){}

	void setType(para_type t){type_=t;}
	const para_type& getType()const{return type_;}

	bool isObservable()const{return type_==para_estimate;}
	bool isSystematic()const{return type_!=para_estimate;}

	const TString& name()const{return name_;}
	void setName(const TString& n){name_=n;}

	const double& getNominalVal()const{return startval_;}
	virtual void setNominalVal(const double & val){startval_=val;}

	const size_t& getAsso()const{return associatedto_;}
	void setAsso(const size_t & s){associatedto_=s;}


	void setStat(const double& s){stat_=s;}
	const double& stat()const{return stat_;}

protected:
	double startval_;
private:
	TString name_;
	para_type type_;

protected:
	size_t associatedto_;

	//these are up to size of global parameters - increases memory usage but makes evaluation faster
    std::vector<size_t> delta_associations_to_pars_;
	std::vector<double> k_sigmasqs_;

	double stat_;
};




#endif /* CROSSSECTION_SUMMER16_PARAMETERS_H_ */
