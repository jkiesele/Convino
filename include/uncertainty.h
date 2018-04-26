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

class uncertainty{
public:
    uncertainty():upvar_(0),downvar_(0){}

    /**
     * up variation and down variation are signed
     */
    uncertainty(const double& up, const double& down):
        upvar_(up),downvar_(down){}

    const double& upVar()const{return upvar_;}
    const double& downVar()const{return downvar_;}

    inline double symm()const{
        double sign=1;
        if(upvar_<downvar_)sign=-1;
        if(fabs(upvar_)>fabs(downvar_)) return sign*fabs(upvar_);
        else return sign*fabs(downvar_);
    }

    double eval(const double& lambda)const{
        if(lambda>=0)return upvar_*lambda;
        else return downvar_*fabs(lambda);
    }

    /**
     * either just a number or a string of the format
     * (+X-Y), where X is the upward variation and Y the downward variation.
     * (-X+Y) or (+X+Y) or (-X-Y) are also possible and the sign information
     * is propagated accordingly
     */
    void readFromString(std::string);

private:
    double upvar_;
    double downvar_;

};

std::ostream& operator<<(std::ostream& os,  uncertainty& u);

#endif /* INCLUDE_UNCERTAINTY_H_ */
