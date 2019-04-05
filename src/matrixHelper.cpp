/*
 * matrixHelper.cpp
 *
 *  Created on: 14 Oct 2016
 *      Author: jkiesele
 */


#include "matrixHelper.h"
#include "TVectorD.h"
#include "helpers.h"
#include "TRandom3.h"
#include "TMatrixDEigen.h"
#include "textFormatter.h"
//#include "averager.h"
#ifdef USE_MP
#include <omp.h>
#endif

#ifdef ARMADILLO
#include <armadillo>
#endif


bool matrixHelper::debug=false;


std::ostream& operator<<(std::ostream& os, const TMatrixD& m){
    std::streamsize save=os.width();
    os<<'\n';
    os.width(4);
    size_t maxnamewidth=5;

    for(int i=0;i<m.GetNrows();i++){
        os.width(maxnamewidth+1);
        os << std::left << i;
        for(int j=0;j<m.GetNcols();j++){
            double entrd=round(m[i][j] ,0.0001);
            std::string entr=toString(entrd);
            os << textFormatter::fixLength(entr,6);
            os <<" ";
        }
        os<<'\n';
    }

    os.width(save);
    return os;

}

matrixHelper::matrixHelper(const TMatrixD &m):strategy_(mhs_weighted),origm_(m),threshold_(0.001),rand_(0),ncalls_(100),
		alwaysset_(-1){

}

bool matrixHelper::checkMatrixPosDefinite(const TMatrixD &m)const{
	if(m.Determinant()>0){
		const TMatrixDEigen eigenvecs(m);
		const TMatrixD eigenvals = eigenvecs.GetEigenValues();
		bool allgood=true;
		for(int i=0;i<eigenvals.GetNrows();i++){
			if(eigenvals[i][i]<=0){
				allgood=false;
			}
		}
		if(allgood)return true;
	}
	return false;
}

void matrixHelper::setKeepEntries(const TMatrixD& in){
	keepentries_.ResizeTo(in.GetNrows(),in.GetNcols());
	keepentries_=in;
}

bool matrixHelper::getPositiveDefinite( TMatrixD& m)const{
	size_t calls=0;
	return getPositiveDefinite(m,calls);
}
bool matrixHelper::getPositiveDefinite( TMatrixD& m, size_t & neededcalls)const{
	//if(!origm_)
	//	throw std::runtime_error("matrixHelper::getPositiveDefinite: original matrix not set");
	m.ResizeTo(origm_);
	m=origm_;

	return getPositiveDefinitePriv(m,&neededcalls);

}





bool matrixHelper::getPositiveDefinitePriv( TMatrixD& morig, size_t * neededcalls)const{


	const int nrows=morig.GetNrows();
	if(invweights_.size()<(size_t)nrows)
		throw std::runtime_error("matrixHelper::getPositiveDefiniteExp: give wiehgt vector with correct size");

	//TMatrixD m2=morig;
	if(strategy_ == mhs_keepentries || strategy_==mhs_keepandweight || strategy_==mhs_keepandadapt)
		if(keepentries_.GetNcols() != keepentries_.GetNrows() || keepentries_.GetNrows()!= nrows)
			throw std::out_of_range("matrixHelper::getPositiveDefinitePriv: keepentries matrix size wrong");

	//produce the weight matrix
	//change for the keepentries_ thing

	TMatrixD weightmatrix=morig;
	double weightproduct=1;
	if(strategy_==mhs_weighted || strategy_==mhs_keepandweight){
		double maxweight=-1;
		double averageweight=0;
		for(int i=0;i<nrows;i++){
			for(int j=0;j<nrows;j++){
				double thisweight=fabs(invweights_.at(i)*invweights_.at(j));
				weightproduct*=thisweight;
				averageweight+=thisweight;
				if(thisweight>maxweight)maxweight=thisweight;
			}
		}
		averageweight/= (float)(nrows*nrows);
		for(int i=0;i<nrows;i++){
			for(int j=0;j<nrows;j++){
				weightmatrix[i][j]=fabs(invweights_.at(i)*invweights_.at(j))/averageweight;//maxweight;
				if(strategy_==mhs_keepandweight && keepentries_[i][j])
					weightmatrix[i][j]=1;
				if(i==j)
					weightmatrix[i][j]=invweights_.at(i)/sqrt(maxweight); //not used but nice for print
			}
		}
		if(debug)
			std::cout << weightmatrix << std::endl;
	}
	else if(strategy_==mhs_keepentries){
		weightmatrix=keepentries_;
	}

	//build a weight matrix
	//	TMatrixD previous=morig;


#ifdef ARMADILLO
	///go to armadillo implementation
	arma::mat m2(nrows,nrows,arma::fill::zeros);
	for(int i=0;i<nrows;i++){
		for(int j=0;j<nrows;j++){
			m2.at(i,j)=morig[i][j];
		}
	}
	arma::mat armaorig=m2;
#endif
	TMatrixD m2=morig;
	TMatrixD armaorig=m2;

	bool allgood=true;
	//	double thisdiff=0;
	//	double lastdiff=0;

	double absfullcorrection=0; //for estimating the convergence behaviour
//	double previousabsfullcorrection=0;
//	double frobdistance=0;

	for(size_t call=0;call<ncalls_+1;call++){
		if(neededcalls)
			*neededcalls=call+1;
		/*
		 * ** reset important values to originals
		 */
		absfullcorrection=0;
		allgood=true;
		if(debug) std::cout << std::endl;

#ifdef ARMADILLO
		arma::mat nosetbackm=m2;
#endif

		for(int i=0;i<nrows;i++){
			for(int j=i;j<nrows;j++){
				if(i==j) continue;
#ifdef ARMADILLO
				double & m_element=m2.at(i,j);
				double & m_trans_element=m2.at(j,i);
				const double& orig_element=armaorig.at(i,j);
#else
				double & m_element=m2[i][j];
				double & m_trans_element=m2[j][i];
				const double& orig_element=armaorig[i][j];
#endif
				if(m_element!=m_element){
					throw std::runtime_error("matrixHelper::getPositiveDefinitePriv: threshold probably set too low, numerically instable");
				}
				if(m_element>1)
					throw std::runtime_error("matrixHelper::getPositiveDefinitePriv: entries>1...");



				bool setback=false;
				//new
				if(strategy_==mhs_keepentries){
					if(keepentries_[i][j]){
						if(fabs(orig_element-m_element)/1e-4 > threshold_){
							setback=true;
						}
					}
				}
				else if(strategy_==mhs_keepandweight){
					if(keepentries_[i][j]){
						if(fabs(orig_element-m_element)/1e-4 > threshold_){
							setback=true;
						}
					}
					else if(weightmatrix[i][j]* fabs(orig_element-m_element) > threshold_){
						setback=true;
					}
				}

				if(setback){
					allgood=false;
					absfullcorrection+=fabs(orig_element-m_element);
					m_element=orig_element;
					m_trans_element=orig_element;
				}


			}
		}
		if(allgood && call){
			break;
		}



		const TMatrixDEigen eigenvecs(m2);
		const TMatrixD eigenvals = eigenvecs.GetEigenValues();
#ifdef ARMADILLO
		arma::vec eigenvals;
		arma::mat eigenvecs;
		arma::eig_sym(eigenvals,eigenvecs,m2);
#endif
		for(int i=0;i<nrows;i++){
#ifdef ARMADILLO
			if(eigenvals.at(i)<=0){
				allgood=false;
			}
#else
			if(eigenvals[i][i]<=0){
				allgood=false;
			}
#endif
		}
		if(allgood){
			if(debug)
				std::cout << "matrix (including forced values) ok"  <<std::endl;
			break;
		}
		//decrease disturbances with number of calls
		double epsilon=1e-2;// * (call*10/ncalls_);//rather large to make it really positive definite
		//if(call<3)
		//	epsilon=1e-6;
		//if(call/ncalls_ <1e3) epsilon=1e-1;
		for(int i=0;i<nrows;i++){
#ifdef ARMADILLO
			if(eigenvals.at(i)<=0){
				arma::mat addm(nrows,nrows,arma::fill::zeros);
				for(int j=0;j<nrows;j++){
					addm.at(j,i)=eigenvecs.at(j,i);
				}

				arma::mat addfull=addm*addm.t();
				addfull*=-(eigenvals[i]-epsilon);

				m2+= addfull;

			}
#else
			if(eigenvals[i][i]<=0){
				TMatrixD addm(nrows,nrows);
				for(int j=0;j<nrows;j++){
					addm[j][i]=eigenvecs.GetEigenVectors () [j][i];
				}
				TMatrixD addmt=addm;
				addmt.T();

				TMatrixD addmatrix=addm;
				addmatrix*=addmt ;
				addmatrix *= -(eigenvals[i][i]-epsilon);
				m2+= addmatrix;
			}
#endif
		}
		//repeat until all positive

#ifdef ARMADILLO
		//make diagonal 1 again
		arma::mat diag(nrows,nrows,arma::fill::zeros);
		for(int l=0;l<nrows;l++){
			diag.at(l,l) = 1/sqrt(m2.at(l,l));
		}
#else
		TMatrixD diag(nrows,nrows);
		for(int l=0;l<nrows;l++){
			diag[l][l] = 1/sqrt(m2[l][l]);
		}

#endif

		m2*=diag;
		diag*=m2;
		m2=diag;

	}

	//end loop
	//feeed back to tmatrix
	if(allgood){
#ifdef ARMADILLO
		for(int i=0;i<nrows;i++){
			for(int j=0;j<nrows;j++){
				morig[i][j]=m2.at(i,j);
			}
		}
#else
		morig=m2;
#endif
	}

	return allgood;
}

