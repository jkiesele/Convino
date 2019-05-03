/*
 * multiGausGenerator.h
 *
 *  Created on: 5 Apr 2019
 *      Author: jkiesele
 */

#ifndef INCLUDE_MULTIGAUSGENERATOR_H_
#define INCLUDE_MULTIGAUSGENERATOR_H_


#include "triangularMatrix.h"
#include <vector>
#include "TMatrixD.h"
#include "TVectorD.h"

class TRandom3;
class multiGausGenerator{
public:
    multiGausGenerator(int seed=123);
    ~multiGausGenerator();

    void setCovariance(const triangularMatrix& m);

    const std::vector<double>& generate()const;

    void reset(){counter_=0;}
private:
    int seed_;
    triangularMatrix m_;
    TMatrixD transform_;
    mutable std::vector<double> tmp_out_;
    mutable TVectorD generated_;
    mutable size_t counter_;
    TRandom3 * rand_;
};



#endif /* INCLUDE_MULTIGAUSGENERATOR_H_ */
