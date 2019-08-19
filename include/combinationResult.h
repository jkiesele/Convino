/*
 * combinationResult.h
 *
 *  Created on: 17 Aug 2016
 *      Author: jkiesele
 */

#ifndef COMBINATIONRESULT_H_
#define COMBINATIONRESULT_H_


#include "triangularMatrix.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>

class TH1;
class TGraphAsymmErrors;

class combiner;
class normaliser;

class combinationResult{
	friend class combiner;
	friend class normaliser;
public:
	combinationResult();

	combinationResult& operator=(const combinationResult& r){
	    copyFrom(r);
	    return *this;
	}
	/**
	 * Returns the number of combined quantities
	 */
	size_t getNCombined()const{return combined_.size();}

	/**
	 * Prints the combination results to a std::ostream
	 * object. This can be either cout or a file (fstream)
	 */
	void printResultOnly(std::ostream& out)const;

	/**
	 * Prints the combination results, and all
	 * correlations prior to the combination and
	 * after the combination to a std::ostream
	 * object. This can be either cout or a file (fstream)
	 */
	void printSimpleInfo(std::ostream& out)const;

	/**
	 * Prints the combination results, all
	 * correlations prior to the combination and
	 * after the combination, and all pulls and
	 * constraints on all parameters to a std::ostream
	 * object. This can be either cout or a file (fstream)
	 */
	void printFullInfo(std::ostream& out)const;



	void printImpactTable(std::ostream& out)const;

	/**
	 * Fills the combination result to a ROOT TH1
	 * object. This is useful for the combination of differential
	 * cross sections. This function cannot define the binning
	 * of the histogram. It must be set before, e.g. by copying
	 * it from the input histograms.
	 */
	void fillTH1(TH1*)const;

	/**
	 * Fills the combination result to a TGraphAsymmErrors.
	 * In case, the initial input had underflow or overflow bins,
	 * these can be removed from the visible area.
	 * The advantage of a TGraphAsymmErrors is that asymmetric uncertainties
	 * are properly displayed.
	 */
	void fillTGraphAsymmErrors(TGraphAsymmErrors*&)const;

	/**
	 * Returns the minimum of the combination likelihood
	 */
	double getChi2min() const {
		return chi2min_;
	}

	int getExcludeBin()const{
	    return excludebin_;
	}

	/**
	 * Returns a vector of the uncertainties on the
	 * combined values in "down" direction.
	 */
	const std::vector<double>& getCombErrdown() const {
		return comberrdown_;
	}

	/**
	 * Returns a vector of the uncertainties on the
	 * combined values in "up" direction.
	 */
	const std::vector<double>& getCombErrup() const {
		return comberrup_;
	}

	/**
	 * Returns a vector of the symmetrised uncertainties using the maximum of up and down variation
	 */
	std::vector<double> getCombSymmErr() const;

	/**
	 * Returns a vector of the combined values.
	 */
	const std::vector<double>& getCombined() const {
		return combined_;
	}

	/**
	 * Returns a vector of names of the combined values
	 */
	const std::vector<TString>& getCombNames() const {
		return combnames_;
	}

	/**
	 * Returns a vector of the constraints on all parameters
	 */
	const std::vector<double>& getConstraints() const {
		return constraints_;
	}

	/**
	 * Returns a vector of the pulls on all parameters
	 */
	const std::vector<double>& getPulls() const {
		return pulls_;
	}


	/**
	 * Returns a correlationMatrix object containing the
	 * input correlations between the estimates
	 */
	const correlationMatrix& getInputEstimateCorrelations() const {
		return orig_meas_correlations_;
	}

	/**
	 * Returns a correlationMatrix object containing the
	 * input correlation assumptions between the uncertainties.
	 */
	const correlationMatrix& getInputSysCorrelations() const {
		return orig_sys_correlations_;
	}

	/**
	 * Returns a correlationMatrix object containing the
	 * correlations between the combined quantities
	 */
	const correlationMatrix& getCombinedCorrelations() const {
		return post_meas_correlations_;
	}

    /**
     * Returns a correlationMatrix object containing the
     * correlations between the combined quantities
     */
    const triangularMatrix& getCombinedCovariance() const;

	/**
	 * Returns a correlationMatrix object containing the
	 * correlations between the uncertainties after the combination
	 */
	const correlationMatrix& getPostSysCorrelations() const {
		return post_sys_correlations_;
	}

    /**
     * Returns a correlationMatrix object containing the
     * correlations between the uncertainties and the combined values after the combination
     * (full matrix)
     */
    const correlationMatrix& getPostAllCorrelations() const {
        return post_all_correlations_;
    }

	/**
	 * Reduces the stored information to combined values, combined names
	 * and the uncertainties, only. Deletes all matrices, pulls and
	 * constraints on uncertainties.
	 */
	void reduce();

protected:
	std::vector<TString> combnames_;
	std::vector<double> combined_, comberrup_,comberrdown_;
	std::vector<std::pair< TString, std::vector<double> > > impacttable_;
	correlationMatrix orig_sys_correlations_;
	correlationMatrix orig_meas_correlations_;
	correlationMatrix post_sys_correlations_,post_meas_correlations_,post_all_correlations_;
	mutable triangularMatrix post_meas_covariance_;
	std::vector<double> pulls_;
	std::vector<double> constraints_;
	double chi2min_;
	bool isdifferential_;
	int excludebin_;

	void copyFrom(const combinationResult& r);

};


#endif /* CROSSSECTION_SUMMER16_COMBINATIONRESULT_H_ */
