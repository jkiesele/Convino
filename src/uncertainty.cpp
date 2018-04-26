/*
 * uncertainty.cpp
 *
 *  Created on: 25 Apr 2018
 *      Author: jkiesele
 */


#include "../include/uncertainty.h"
#include "TString.h"


void uncertainty::readFromString(std::string s){
    /*
     * either just a number or a string of the format
     * (+X-Y), where X is the upward variation and Y the downward variation.
     * (-X+Y) or (+X+Y) or (-X-Y) are also possible and the sign information
     * is propagated accordingly
     */
    if(s.length()<1)
        throw std::runtime_error("uncertainty::readFromString: string empty");
    if(s.find("(")==std::string::npos && s.find(")")==std::string::npos){ //simple format
        upvar_=atof(s.data());
        downvar_=-upvar_;
        return;
    }

    if(s.at(s.length()-1)!=')')
        throw std::runtime_error("uncertainty::readFromString: string format wrong: "+s);
    if(s.at(0)!='(')
        throw std::runtime_error("uncertainty::readFromString: string format wrong: "+s);
    if(s.length()<6)
        throw std::runtime_error("uncertainty::readFromString: string format wrong: "+s);

    TString tst=s;
    tst.ReplaceAll(" ","");
    s=tst.Data();
    //remove parentheses
    s=s.substr(1,s.length()-2);

    double secondsign=-1;
    char separator='-';
    size_t seppos=s.find_last_of('-');
    if(seppos==std::string::npos || seppos==0){
        separator='+';
        seppos=s.find_last_of('+');
        secondsign=1;
    }
    std::string firststr=s.substr(0,seppos);
    std::string secstr=s.substr(seppos,s.length());

    if(firststr.at(0)=='+')
        firststr=firststr.substr(1,firststr.length());
    if(secstr.at(0)=='+')
        secstr=secstr.substr(1,secstr.length());


    upvar_ = atof(firststr.data());
    downvar_ = atof(secstr.data());

}


std::ostream& operator<<(std::ostream& os,  uncertainty& u){
    if(u.upVar()>=0)
        os << "(+" << u.upVar();
    else
        os << "(" << u.upVar();
    if(u.downVar()>=0)
        os << "+" << u.downVar();
    else
        os << u.downVar();
    os << ")";

    return os;
}

