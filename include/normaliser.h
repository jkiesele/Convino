/*
 * normaliser.h
 *
 *  Created on: 5 Apr 2019
 *      Author: jkiesele
 */

#ifndef INCLUDE_NORMALISER_H_
#define INCLUDE_NORMALISER_H_

#include "combinationResult.h"
#include "triangularMatrix.h"

class TH1D;
class TH2D;
class normaliser{
public:

    normaliser():iterations_(1e6),seed_(0),addbin_(-1){}

    void setInput(const combinationResult& res);
    void setInput(const TH1D* , const TH2D* );
    void setInput(const std::vector<double>* , const triangularMatrix* );


    void setIterations(size_t iterations){
        iterations_=iterations;
    }
    void setSeed(int seed){
        seed_=seed;
    }

    void clear(){
        out_res = combinationResult();
    }

    void addFloatingBin(int number, const TString& name=""){
        addbin_ = number;
        addbin_name = name;
        if(name.Length()<1)
            addbin_name="add_bin";
    }

    combinationResult getNormalised(double addUncertFraction=0., int bin=-1)const;

    TH1D * getNormalisedTH1D()const;
    TH2D * getNormalisedCovarianceTH2D(double addUncertFraction=0., int bin=-1)const;

private:
    //helper
    double normaliseVector(std::vector<double> & v)const;

    size_t iterations_;
    int seed_;
    combinationResult in_res_;
    mutable combinationResult out_res;
    int addbin_;
    TString addbin_name;
};


#endif /* INCLUDE_NORMALISER_H_ */
