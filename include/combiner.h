/*
 * combiner.h
 *
 *  Created on: 30 Jun 2016
 *      Author: jkiesele
 */

#ifndef COMBINER_H_
#define COMBINER_H_



#ifndef NOCOMPILE

#include "triangularMatrix.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "parameters.h"
#include <algorithm>
#include <stdexcept>
#include "combinationResult.h"
#include "measurement.h"

class TH1;
class pseudoGenerator;
class fitfunctionBase;
/**
 * class to combine. Implement the function to be minimized here
 *
 *  add also possibility to input syst correlations and stat correlations between measurements separately
 *  (for diff xsec)
 *
 */
class combiner{
	friend class fitfunctionBase;
public:
	combiner():lh_mod_(lh_mod_neyman),npars_(0),nest_(0),lowestchi2_(1e19),isdifferential_(false),hasUF_(false),hasOF_(false){}
	~combiner(){clear();}
	enum lh_mod{lh_mod_neyman,lh_mod_pearson};

	/// ---- C++ interface -----
	/**
	 * Defines the mode, options:
	 * - combiner::lh_mod_neyman
	 * - combiner::lh_mod_pearson
	 */
	void setMode(lh_mod m){
		lh_mod_=m;
	}


	/**
	 * Adds a measurement object to the
	 * combiner such that it can be combined
	 * with another measurement
	 */
	void addMeasurement( measurement);


	/**
	 * Will be prepended to each output file
	 */
	void setOutputPrefix(const std::string& prefix){outprefix_=prefix+"_";}

	/**
	 * Set the correlation assumption between uncertainties
	 * of two measurements based on the uncertainty names.
	 * Use this function preferably to avoid
	 * wrongly-assigned correlations
	 */
	void setSystCorrelation(const TString & namea, const TString& nameb, const double& coeff);

	/**
	 * Set the correlation assumption between uncertainties
	 * of two measurements based on the uncertainty names.
	 * Use this function only if code performance is an issue.
	 */
	void setSystCorrelation(const size_t & idxa, const size_t& idxb, const double& coeff);

	/**
	 * Starts the combination procedure and returns a
	 * combinationResult object that stores the result of the
	 * combination and all input paramters.
	 */
	combinationResult combine()const;



	///// ----- helpers for text-based interface -----
	void readConfigFile(const std::string & filename);
	std::vector<std::vector<combinationResult> >
	scanCorrelationsIndep(std::ostream& out, const combinationResult& nominal, const std::string& outdir="")const;
	void printCorrelationMatrix()const;
	void associate(const TString & a, const TString& outname);

	static bool debug;
private:


	void createExternalCorrelations();
	void associatePriv(const TString & a, const TString& outname);

	combinationResult combinePriv();



	class correlationscan{
	public:
		correlationscan():
			idxa(std::string::npos),idxb(std::string::npos),nominal(0),low(0),high(0),steps(10){}
		correlationscan(const size_t& cidxa,const size_t& cidxb,const float& cnominal, const float& clow, const float& chigh,  size_t csteps=10):
			idxa(cidxa),idxb(cidxb),nominal(cnominal),low(clow),high(chigh),steps(csteps){}
		size_t idxa,idxb;
		float nominal;
		float low;
		float high;
		size_t steps;
	};

	/*
	 * set a measurement to be combined to <outname>
	 * e.g.
	 * associate(xsec_CMS_8TeV, xsec_8TeV);
	 * associate(xsec_ATLAS_8TeV, xsec_8TeV);
	 *
	 * associate(xsec_CMS_7TeV, xsec_7TeV);
	 * associate(xsec_ATLAS_7TeV, xsec_7TeV);
	 *
	 * etc..
	 */

	void addMeasurement(const std::string& infile);

	lh_mod lh_mod_;

	//for scans



	//mapping between fit indices and estimate_correlations_ indices
	//std::vector< std::pair<TString, std::vector<size_t> > > associations_;

	/// new input for simple delta cov delta combination


	void clear();

	size_t npars_,nest_;

	std::vector<correlationscan> syst_scanranges_;
	correlationMatrix external_correlations_;
	TMatrixD inv_priors_;
	std::vector<parameter> allsysparas_;

	std::vector<measurement> measurements_;
	std::vector<parameter>  allparas;
	std::vector<std::pair< TString, std::vector<TString> > > tobecombined_;

	//measurement allmeas_;

	//size_t npars_,nsys_,nest_;

	//bool measurements_readin_;
	std::string configfile_;
	double lowestchi2_;
	std::string outprefix_;

	//double resolvethresh_;

	bool isdifferential_,hasUF_,hasOF_;

	static const double maxcorr_;
};


#endif//nocompile


#endif /* CROSSSECTION_SUMMER16_COMBINER_H_ */
