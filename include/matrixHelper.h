/*
 * matrixHelper.h
 *
 *  Created on: 14 Oct 2016
 *      Author: jkiesele
 */

#ifndef INCLUDE_MATRIXHELPER_H_
#define INCLUDE_MATRIXHELPER_H_

#include "triangularMatrix.h"
#include "TMatrixD.h"
#include <vector>
#include "TRandom3.h"


class matrixHelper{
public:
	matrixHelper(const TMatrixD &);

	enum matrixHelper_stategy {mhs_weighted, mhs_keepentries, mhs_keepandadapt, mhs_keepandweight};

	void setInvWeights(const std::vector<double>& w){ invweights_=w;}
	void setKeepEntries(const TMatrixD&);
	void setStrategy(matrixHelper_stategy s){strategy_=s;}\
	void setThreshold(const double& th){threshold_=th;}

	bool getPositiveDefinite( TMatrixD& m)const;
	bool getPositiveDefinite( TMatrixD& m, size_t & neededcalls)const;

	bool checkMatrixPosDefinite(const TMatrixD &m)const;

	void setRandom(TRandom3* rand){rand_=rand;}

	void setNCalls(size_t n){ncalls_=n;}

	void setAlwaysSetBack(double val){alwaysset_=val;}


	static bool debug;
private:
	std::vector<double> invweights_;
	matrixHelper_stategy strategy_;
	const TMatrixD& origm_;
	TMatrixD keepentries_;
	double threshold_;
	TRandom3* rand_;
	size_t ncalls_;
	double alwaysset_;

	bool getPositiveDefiniteDist( TMatrixD& m, size_t * neededcalls)const;
	bool getPositiveDefiniteMult( TMatrixD& m, size_t * neededcalls)const;


	bool getPositiveDefinitePriv( TMatrixD& m, size_t * neededcalls)const;

};

#endif /* INCLUDE_MATRIXHELPER_H_ */
