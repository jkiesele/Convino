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





class parameter{
public:

	enum para_type{para_unc_absolute,para_unc_relative,para_unc_lognormal,para_estimate};

	parameter():startval_(0),type_(para_unc_absolute),associatedto_(0){}
	virtual ~parameter(){}

	void setType(para_type t){type_=t;}
	const para_type& getType()const{return type_;}

	bool isObservable()const{return type_==para_estimate;}
	bool isSystematic()const{return type_!=para_estimate;}

	bool isRelative()const{return type_==para_unc_relative || type_==para_unc_lognormal;}

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
	size_t associatedto_;
	double stat_;
};




#endif /* CROSSSECTION_SUMMER16_PARAMETERS_H_ */
