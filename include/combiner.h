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
    combiner():lh_mod_(lh_mod_neyman),npars_(0),nest_(0),
    lowestchi2_(1e19),isdifferential_(false),
    normalised_input_(false),
    excludebin_(-1),
    normaliseinfit_(false),normalisewithtoys_(false),
    fastmode_(false),
    auglagrangemu_(1),auglagrangelambda_(1),
    allcontours_(true){}
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

    void setExcludeBin(int bin);

    int getExcludeBin()const{
        return excludebin_;
    }


    void setFastMode(bool fm);

    void setAllContours(bool set){allcontours_=set;}
    /**
     * Starts the combination procedure and returns a
     * combinationResult object that stores the result of the
     * combination and all input paramters.
     */
    combinationResult combine()const;



    ///// ----- helpers for text-based interface -----
    void readConfigFile(const std::string & filename);


    void addToImpactTable(const TString& impactdescription, const std::vector<TString>& contributing_uncertatinties);

    /**
     * can be used to read-in correlations from an external file for text based and C++ interface
     */
    void readCorrelationFile(const std::string & filename, bool requiremarkers=false);


    void readImpactFile(const std::string & filename, bool requiremarkers=false);

    std::vector<std::vector<combinationResult> >
    scanCorrelations(std::ostream& out, const combinationResult& nominal, const std::string& outdir="")const;

    std::vector<combinationResult>
    scanExcludeBins(std::ostream& out, const combinationResult& nominal)const;
    void printCorrelationMatrix()const;
    void associate(const TString & a, const TString& outname);

    static bool debug;

    void setNormaliseInFit(bool normit){
        normaliseinfit_=normit;
    }
    void setNormaliseWithToys(bool normit){
        normalisewithtoys_=normit;
    }

    bool isNormalisedDifferentialInput()const{
        return isdifferential_ && normalised_input_;
    }
private:

    void checkConsistency()const;

    void createExternalCorrelations();
    void associatePriv(const TString & a, const TString& outname);

    combinationResult combinePriv();



    class single_correlationscan{
    public:
        single_correlationscan():
            idxa(std::string::npos),idxb(std::string::npos),nominal(0),low(0),high(0){}
        single_correlationscan(const size_t& cidxa,const size_t& cidxb,const float& cnominal, const float& clow, const float& chigh,  size_t csteps=10):
            idxa(cidxa),idxb(cidxb),nominal(cnominal),low(clow),high(chigh){}

        float scanVal(size_t i)const{
            if(i>nPoints())
                throw std::out_of_range("single_correlationscan::scanVal");
            return low+ (double)i*(high-low)/(double)(nPoints()-1);
        }
        static size_t nPoints(){
            return 6;
        }
        size_t idxa,idxb;
        float nominal;
        float low;
        float high;
        TString namea, nameb;

    };
    //very simple
    class correlationscan{
    public:
        void setName(const std::string& n){
            name_=n;
        }
        const std::string& name()const{
            return name_;
        }
        void add(const single_correlationscan&s ){
            scans_.push_back(s);
        }
        const single_correlationscan& get(size_t i)const {
            return scans_.at(i);
        }
        void set(size_t i, const single_correlationscan& s) {
            scans_.at(i)=s;
        }
        size_t size()const{
            return scans_.size();
        }
        bool isSingle()const{
            return size()<2;
        }
        float getLowest()const{
            if(isSingle())
                return std::min(scans_.at(0).low,scans_.at(0).high);
            else
                return 0;
        }
        float getHighest()const{
            if(isSingle())
                return std::max(scans_.at(0).low,scans_.at(0).high);
            else
                return 1;
        }
    private:
        std::vector<single_correlationscan> scans_;
        std::string name_;
    };



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
    std::vector<std::pair< TString, std::vector<TString> > > impacttable_;


    //measurement allmeas_;

    //size_t npars_,nsys_,nest_;

    //bool measurements_readin_;
    std::string configfile_;
    double lowestchi2_;
    std::string outprefix_;

    //double resolvethresh_;

    bool isdifferential_, normalised_input_;
    int excludebin_;
    bool normaliseinfit_,normalisewithtoys_;

    static const double maxcorr_;

    static bool dummyrun_;//< for debugging purposes
    bool fastmode_;


    double auglagrangemu_,auglagrangelambda_;

    bool allcontours_;

};


#endif//nocompile


#endif /* CROSSSECTION_SUMMER16_COMBINER_H_ */
