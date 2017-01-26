/*
 * measurement.h
 *
 *  Created on: 12 Dec 2016
 *      Author: jkiesele
 */

#ifndef INCLUDE_MEASUREMENT_H_
#define INCLUDE_MEASUREMENT_H_

#include <string>
#include <vector>
#include "triangularMatrix.h"
#include "parameters.h"
#include "TString.h"
#include "TH1D.h"
#include "namedMatrix.h"

class measurement{
public:
	measurement();
	measurement(const std::string& infile);
	~measurement();
	measurement& operator=(const measurement& r);
	measurement(const measurement& );

	///// test-based interface /////

	/**
	 * Reads the information from a measurement file.
	 * For its description, please refer to the manual
	 * or examples in the "examples" directory
	 */
	void readFromFile(const std::string& infile);

	///// ROOT Interface /////

	/**
	 * Set the measured values from a root TH1D.
	 * It must contain statistical uncertainties.
	 * If statistical uncertainties in the overflow
	 * or underflow bins are present, these bins are also
	 * considered for the combination.
	 *
	 * If this function is used, the estimates will
	 * be named automatically and later combined bin-
	 * by-bin (accounting for optional correlations)
	 *
	 * Combining only a subset of bins:
	 *   It is not possible to combine only a set of
	 *   bins directly. The number of bins must be the
	 *   same in all measurement objects.
	 *   If it is desired to only combine a subset of all
	 *   bins, the other bins must be defined, too. Their
	 *   impact can be made negligible, by assigning a
	 *   very large statistical uncertainty to them
	 */
	void setMeasured(const TH1D* h);

	/**
	 * Adds a systematic uncertainty with a name. The name
	 * must be unique. The TH1D histogram must contain the
	 * varied distribution (nominal values + uncertainties)
	 */
	void addSystematics(const TString& name, const TH1D* h);


	//////// C++ interface ///////

	/**
	 * Defines the measured values and their statistical
	 * uncertainties. Once this method is used, it is not
	 * possible to combine a different number of
	 * estimates in each measurement directly.
	 *
	 * If this is desired they should be provided with a hessian
	 * using the method setHessian(). However, a workaround
	 * is the following:
	 * If it is desired to only combine a subset of all
	 * estimates, the other estimates must be defined, too.
	 * Their impact can be made negligible, by assigning a
	 * very large statistical uncertainty to them.
	 */
	void setMeasured(const std::vector<double>& meas, const std::vector<double>& stat);

	/**
	 * Adds a systematic uncertainty which is uncorrelated
	 * to all others. Correlations between the uncertainties
	 * can be defined only by using the setHessian() method
	 * The vector of values must contain the varied entry
	 * (nominal + uncertainty)
	 */
	void addSystematics(const TString& name, std::vector<double> variedmeas);

	/**
	 * Defines the correlation between the estimates i and j.
	 * If the correlation was already defined by the Hessian,
	 * an exception is thrown.
	 */
	void setEstimateCorrelation(const size_t& i, const size_t& j, const double& val);


	/**
	 * Defines the Hessian of the measurement.
	 * Additional, non-correlated uncertainties
	 * can be added using the method addSystematics().
	 * If this is done, it is important that the
	 * order of the estimates matches.
	 */
	void setHessian(const triangularMatrix&h);

	/**
	 * Returns the Hessian matrix (if set)
	 */
	const triangularMatrix& getHessian(){return H_;}

	/**
	 * Defines the type of a parameter with index <idx>
	 * to be either:
	 *   parameter::para_estimate
	 *   parameter::para_unc_absolute
	 *   parameter::para_unc_relative
	 */
	void setParameterType(const size_t & idx, parameter::para_type type);

	/**
	 * Defines the type of a parameter with name <name>
	 * to be either:
	 *   parameter::para_estimate
	 *   parameter::para_unc_absolute
	 *   parameter::para_unc_relative
	 */
	void setParameterType(const TString & name, parameter::para_type type);

	/**
	 * Defines the nominal value of a parameter with index <idx>
	 * For uncertainties, the entry will be ignored
	 */
	void setParameterValue(const size_t & idx, const double& val);

	/**
	 * Defines the nominal value of a parameter with name <name>
	 * For uncertainties, the entry will be ignored
	 */
	void setParameterValue(const TString & name, const double& val);


	/////// interface to combined class. Should not be used by the user ///////

	void associateEstimate(const size_t & est_idx, const size_t & comb_idx);
	void associateEstimate(const TString & est_name, const size_t & comb_idx);
	void associateAllLambdas(const triangularMatrix& fullLambdaCorrs);
	const std::vector<parameter> & getParameters()const{return paras_;}
	void setup();
	double evaluate(const double* pars, double* df, const bool& pearson, const size_t& maxidx)const;
	const std::vector<parameter> & getLambdas()const{return lambdas_;}//for priors
	const std::vector<parameter> & getEstimates()const{return x_;}
	const parameter& getParameter(const TString& name)const;
	bool isDifferential()const{return isDifferential_;}
	const bool& hasUF()const{return hasUF_;}
	const bool& hasOF()const{return hasOF_;}

	static bool debug;

private:

	void copyFrom(const measurement& r);

	bool setup_;

	triangularMatrix H_;
	std::vector<parameter> paras_;

	/*
	 * Matrices from Hessian
	 */
	std::vector<std::vector<double> > M_, kappa_, tildeC_;

	/*
	 * Matrices to descrbe initial pseudo-likelihood
	 */
	std::vector<std::vector<double> > LM_, Lk_, LD_;


	std::vector<parameter> lambdas_;
	std::vector<parameter> x_;


	namedMatrix c_external_;

	/*
	 * Temps for root/c++ interface
	 */
	correlationMatrix est_corr_;

	bool hasUF_, hasOF_;

	static size_t nobjects_;

	bool isDifferential_;
};



#endif /* INCLUDE_MEASUREMENT_H_ */
