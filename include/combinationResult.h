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

class TH1;
class TGraphAsymmErrors;

class combiner;

class combinationResult{
	friend class combiner;
public:
	combinationResult();
	size_t getNCombined()const{return combined_.size();}

	void printResultOnly(std::ostream& out)const;

	void printSimpleInfo(std::ostream& out)const;

	void printFullInfo(std::ostream& out)const;

	void fillTH1(TH1*)const;

	void fillTGraphAsymmErrors(TGraphAsymmErrors*, bool cutUFOF=false)const;

	double getChi2min() const {
		return chi2min_;
	}

	const std::vector<double>& getCombErrdown() const {
		return comberrdown_;
	}

	const std::vector<double>& getCombErrup() const {
		return comberrup_;
	}

	const std::vector<double>& getCombined() const {
		return combined_;
	}

	const std::vector<TString>& getCombNames() const {
		return combnames_;
	}

	const std::vector<double>& getConstraints() const {
		return constraints_;
	}

	const correlationMatrix& getInputEstimateCorrelations() const {
		return orig_meas_correlations_;
	}

	const correlationMatrix& getInputSysCorrelations() const {
		return orig_sys_correlations_;
	}

	const correlationMatrix& getCombinedCorrelations() const {
		return post_meas_correlations_;
	}

	const correlationMatrix& getPostSysCorrelations() const {
		return post_sys_correlations_;
	}

	const std::vector<double>& getPulls() const {
		return pulls_;
	}

protected:
	std::vector<TString> combnames_;
	std::vector<double> combined_, comberrup_,comberrdown_;
	correlationMatrix orig_sys_correlations_;
	correlationMatrix orig_meas_correlations_;
	correlationMatrix post_sys_correlations_,post_meas_correlations_;
	std::vector<double> pulls_;
	std::vector<double> constraints_;
	double chi2min_;
	bool isdifferential_,hasUF_,hasOF_;
};


#endif /* CROSSSECTION_SUMMER16_COMBINATIONRESULT_H_ */
