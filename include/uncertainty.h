/*
 * uncertainty.h
 *
 *  Created on: 25 Apr 2018
 *      Author: jkiesele
 */

#ifndef INCLUDE_UNCERTAINTY_H_
#define INCLUDE_UNCERTAINTY_H_

#include <iostream>
#include <string>
#include <math.h>
#include "parameters.h"

/*
 *
 * Make this inherit from parameter and then cast to it in the loop?
 *
 */


class uncertainty: public parameter{
public:
    uncertainty():assoestname_(""),upvar_(0),downvar_(0),uncofunc2_(-100),uncofunc_paridx_(0){}

    /**
     * up variation and down variation are signed
     */
    uncertainty(const double& up, const double& down):assoestname_(""),
            upvar_(up),downvar_(down),uncofunc2_(-100),uncofunc_paridx_(0){}

    const double& upVar()const{return upvar_;}
    const double& downVar()const{return downvar_;}

    bool isRelative()const{return type_==para_unc_relative || type_==para_unc_lognormal;}

    void setAssoEstimateName(const TString& name){
        assoestname_=name;
    }
    const TString& getAssoEstimateName()const{
        return assoestname_;
    }

    const size_t& getSubParaAsso()const{
        return uncofunc_paridx_;
    }

    void setSubParaAsso(size_t asso){
        uncofunc_paridx_=asso;
    }

    inline double symm()const{
        double sign=1;
        if(upvar_<downvar_)sign=-1;
        if(fabs(upvar_)>fabs(downvar_)) return sign*(double)fabs(upvar_);
        else return sign*fabs(downvar_);
    }

    inline bool hasUnvOfUnc()const{
        return uncofunc2_>0;
    }


    void setUncSigmaSq(double sigmasq){
        uncofunc2_=sigmasq;
    }

    const double& getUncSigmaSq()const{
        return uncofunc2_;
    }


    inline double eval(const double* pars)const{
        if(hasUnvOfUnc()){
            return eval_priv(pars[uncofunc_paridx_]);
        }
        else
            return eval_priv(pars[associatedto_]);
    }

    inline const double& getsqrtpenalty(const double* pars)const{
        return pars[associatedto_];
    }

    double getUncOfUncPenalty(const double* pars)const{
        if(hasUnvOfUnc()){
            double diff = pars[uncofunc_paridx_] - pars[associatedto_];
            return diff*diff/uncofunc2_;
        }
        else{ return 0;}
    }


    /**
     * either just a number or a string of the format
     * (+X-Y), where X is the upward variation and Y the downward variation.
     * (-X+Y) or (+X+Y) or (-X-Y) are also possible and the sign information
     * is propagated accordingly
     */
    void readFromString(std::string);

private:
    inline double eval_priv(const double& lambda)const{
        if(lambda>=0)return upvar_*lambda;
        else return downvar_*(double)fabs(lambda);
    }
    TString assoestname_;

    double upvar_;
    double downvar_;

    double uncofunc2_;
    size_t uncofunc_paridx_;

};

std::ostream& operator<<(std::ostream& os,  uncertainty& u);

#endif /* INCLUDE_UNCERTAINTY_H_ */
