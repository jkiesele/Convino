/*
 * multiGausGenerator.cpp
 *
 *  Created on: 5 Apr 2019
 *      Author: jkiesele
 */


#include "multiGausGenerator.h"
#include "TMatrixD.h"
#include "matrixHelper.h"
#include <iostream>
#include "TMatrixDEigen.h"
#include "TRandom3.h"


multiGausGenerator::multiGausGenerator(int seed):seed_(seed),counter_(0){

    rand_ = new TRandom3(seed_);

}
multiGausGenerator::~multiGausGenerator(){
    delete rand_;
}

const std::vector<double>& multiGausGenerator::generate()const{

    counter_++;
    for(size_t i=0;i<m_.size();i++)
        generated_[i] = rand_->Gaus(0,1);

    TVectorD out = transform_ * generated_;
    for(size_t i=0;i<m_.size();i++)
        tmp_out_[i]=out[i];
    return tmp_out_;
}

void multiGausGenerator::setCovariance(const triangularMatrix& m_in){

    TMatrixD tm;
    m_in.toTMatrix(tm);
    matrixHelper helper(tm);
    if (! helper.checkMatrixPosDefinite(tm)){
        std::cout << m_in << std::endl;
        throw std::runtime_error("multiGausGenerator::setCovariance: input covariance must be positive definite. Check input.");
    }

    m_=m_in;
    tmp_out_.clear();
    tmp_out_.resize(m_.size());
    generated_.ResizeTo(m_.size());



    const TMatrixDEigen e_eigenvecs(tm);
    const TMatrixD eigenvals(e_eigenvecs.GetEigenValues());
    const TMatrixD eigenvecs(e_eigenvecs.GetEigenVectors());


    TMatrixD sqrt_eigenvals_diag(eigenvals);
    for(int i=0;i<eigenvals.GetNrows(); i++)
        for(int j=0;j<eigenvals.GetNrows(); j++){
            sqrt_eigenvals_diag[i][j]=std::sqrt(eigenvals[i][j]);
        }

    transform_.ResizeTo(eigenvecs);
    transform_ = eigenvecs * sqrt_eigenvals_diag;



}
