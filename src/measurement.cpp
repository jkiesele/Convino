/*
 * measurement.cpp
 *
 *  Created on: 12 Dec 2016
 *      Author: jkiesele
 */


#include "measurement.h"
#include "TMatrixD.h"
#include "matrixHelper.h"
#include "helpers.h"
#include "fileReader.h"
#include <iostream>
#include <math.h>
#ifdef USE_MP
#include <omp.h>
#endif

bool measurement::debug=false;
size_t measurement::nobjects_=0;

measurement::measurement():setup_(false),isDifferential_(false),
        isnormalisedinput_(false),
        normalisation_(-1),
        normalise_(false),
        bypass_logic_check_(false){
    this_obj_counter_=nobjects_;
	nobjects_++;
	if(debug)
	    std::cout << "created measurement "<< this_obj_counter_<< std::endl;
}

measurement::measurement(const std::string& infile):setup_(false),isDifferential_(false),
        bypass_logic_check_(false){
	readFromFile(infile);
    this_obj_counter_=nobjects_;
	nobjects_++;
}
measurement::~measurement(){
}

measurement& measurement::operator=(const measurement& r){
	copyFrom(r);
	return *this;
}

measurement::measurement(const measurement&r ){
	copyFrom(r);
}

/////////// text based interface //////////

void measurement::readFromFile(const std::string& infile){
	//search for Hessian
	triangularMatrix hessian;

	fileReader fr;
	fr.setDelimiter(" ");
	fr.setComment("#");
	isDifferential_=false;




	if(fr.hasNonEmptyMarker(infile,"hessian")){
		hessian.readFromFile(infile, "hessian");
		setHessian(hessian);
	}
	else if(fr.hasNonEmptyMarker(infile,"correlation matrix")){

		std::vector<double> constr;
		triangularMatrix corr;

		corr.readFromFile(infile,"correlation matrix");

		setHessian(corr); //creates parameters

		fr.clear();
		fr.setStartMarker("[correlation matrix]");
		fr.setEndMarker("[end correlation matrix]");
		fr.readFile(infile);

		for(size_t i=0;i<fr.nLines();i++){
			size_t nentries = fr.nEntries(i);
			if(nentries>1 && nentries < 3)
				throw std::runtime_error("measurement::readFromFile: please provide the fitted constraints for all parameters when using a correlation matrix");

			TString sysname = fr.getData<TString>(i,0);

			TString constraintstr = fr.getData<TString>(i,1);

			/*
			 * convert to double
			 */
			constraintstr.ReplaceAll("(","");
			constraintstr.ReplaceAll(")","");
			double constraint = constraintstr.Atof();

			constr.push_back(constraint);

		}

		//create hessian from correlation matrix
		for(size_t i=0;i<corr.size();i++){
			for(size_t j=i;j<corr.size();j++){
				corr.setEntry(i,j, corr.getEntry(i,j) * constr.at(i) * constr.at(j));
			}
		}
		double det;
		corr.invert(det);
		if(det<0)
			throw std::runtime_error("measurement::readFromFile: covariance not positive definite.\ntry to give Hessian directly or increase precision.");
		H_=corr;

	}
	if(fr.hasNonEmptyMarker(infile,"not fitted")){
		//add externalised uncertainties
		fr.setStartMarker("[not fitted]");
		fr.setEndMarker("[end not fitted]");
		fr.readFile(infile);
		//read the input
		std::vector<TString> estnames, sysnames, allnames;
		std::vector<std::vector<uncertainty> > unc;
		std::vector<double> stat(fr.nLines()-1);
		bool hasstat=false;


		for(size_t i=0;i<fr.nEntries(0);i++){
			TString name=fr.getData<TString>(0,i);
			if(std::find(sysnames.begin(),sysnames.end(),name) != sysnames.end())
				throw std::runtime_error("measurement::readFromFile: wrong format in \"not fitted\": please provide unique names for systematics");
			if(name!="stat")
				sysnames.push_back(name);
			else
				hasstat=true;
		}

		size_t expectentries=sysnames.size()+1;
		if(hasstat) expectentries++;

		for(size_t l=1;l<fr.nLines();l++){
			if(fr.nEntries(l) != expectentries )
				throw std::runtime_error("measurement::readFromFile: wrong format in \"not fitted\"");

			std::vector<uncertainty> vunc;
			for(size_t e=0;e<fr.nEntries(l);e++){
			    TString thisentry = fr.getData<TString>(l,e);
				if(!e){
					estnames.push_back(fr.getData<TString>(l,e));
				}
				else{

					if(!hasstat || e<fr.nEntries(l)-1){
	                    uncertainty thisunc;
	                    thisunc.readFromString(fr.getData<std::string>(l,e));
						vunc.push_back(thisunc);
					}
					else{
						stat.at(l-1)=fr.getData<double>(l,e);}
				}
			}
			unc.push_back(vunc);
		}


		c_external_ = namedMatrix<uncertainty>(estnames,sysnames);
		for(size_t i=0;i<c_external_.nRows();i++){
			for(size_t j=0;j<c_external_.nCols();j++){
				c_external_.setEntry(i,j,unc.at(i).at(j));
			}
		}

		//check if already existing paras
		if(paras_.size()){
			for(const auto& e: estnames){
				bool foundest=false;
				for(const auto& p:paras_){
					if(e==p.name()){
						foundest=true;
						break;
					}
				}
				if(!foundest)
					throw std::runtime_error("measurement::readFromFile: wrong format in \"not fitted\": estimate not found");
			}

			//all estimates exists
			//check for ambig syst.
			for(const auto& p:paras_){
				if(std::find(sysnames.begin(),sysnames.end(),p.name())!=sysnames.end())
					throw std::runtime_error("measurement::readFromFile: wrong format in \"not fitted\": non-ambiguous uncertainty names");
			}
			for(const auto& s:sysnames){
				parameter np;
				np.setType(parameter::para_unc_absolute);
				np.setName(s);
				paras_.push_back(np);
			}
		}
		else{
			for(size_t e=0;e<estnames.size();e++){
				parameter pe;
				pe.setType(parameter::para_estimate);
				pe.setStat(stat.at(e));
				if(!stat.at(e))
				    throw std::logic_error("measurement::readFromFile: statistical uncertainty of at least one estimate 0. Check the input.");
				pe.setName(estnames.at(e));
				paras_.push_back(pe);
			}
			for(const auto& e:sysnames){
				parameter pe;
				pe.setName(e);
				pe.setType(parameter::para_unc_absolute);
				paras_.push_back(pe);
			}
			est_corr_ = correlationMatrix(estnames);
		}
	}


	if(paras_.size()<1)
		throw std::runtime_error("measurement::readFromFile: no description of the parameters found in "+infile);

	fr.setComment("#");
	fr.setStartMarker("[estimates]");
	fr.setEndMarker("[end estimates]");
	fr.readFile(infile);

	size_t n_est = fr.getValue<int>("n_estimates");
	for(size_t i=0;i<n_est;i++){
		std::string istr     = toString(i);
		TString estimatename = fr.getValue<TString>("name_"+istr);
		const double value   = fr.getValue<double>("value_"+istr);

		size_t idx=SIZE_MAX;
		for(size_t e=0;e<paras_.size();e++){
			if(paras_.at(e).name() == estimatename)
				idx= e;
		}
		if(idx==SIZE_MAX)
			throw std::runtime_error("measurement::readFromFile: \""+ (std::string)estimatename+"\"  not found in "+infile);

		paras_.at(idx).setNominalVal(value);
		paras_.at(idx).setType(parameter::para_estimate);

	}

	fr.clear();
	fr.setStartMarker("[systematics]");
	fr.setEndMarker("[end systematics]");
	fr.setDelimiter("=");
	fr.readFile(infile);

	for(size_t i=0;i<fr.nLines();i++){

		if(fr.nEntries(i)<2)continue;

		TString sysname = fr.getData<TString>(i,0);
		TString type    = fr.getData<TString>(i,1);
		size_t idx=SIZE_MAX;
		for(size_t e=0;e<paras_.size();e++){
			if(paras_.at(e).name() == sysname)
				idx= e;
		}
		if(idx==SIZE_MAX)
			throw std::runtime_error("measurement::readFromFile: \""+ (std::string)sysname+"\" uncertainty name not found in "+infile);


		if(type=="lognormal")
			paras_.at(idx).setType(parameter::para_unc_lognormal);
		else if(type=="absolute")
			paras_.at(idx).setType(parameter::para_unc_absolute);
		else if(type=="relative")
			paras_.at(idx).setType(parameter::para_unc_relative);
		else
			throw std::runtime_error(
					("combiner::measurement::readFromFile: setting systematics types: type "+
							type+" not recognized\noptions: lognormal, absolute, relative (default)\nthe latter two imply gaussian priors").Data());
	}



}



////////// C++ interface //////////

void measurement::setMeasured(const TH1D* h){

	int nbins=h->GetNbinsX()+1;
	int minbin=1,maxbin=nbins;
	std::vector<double> meas,stat;
	for(int i=minbin;i<maxbin;i++){
		meas.push_back(h->GetBinContent(i));
		stat.push_back(h->GetBinError(i));
	}
	setMeasured(meas,stat);
}

void measurement::setMeasured(const std::vector<double>& meas, const std::vector<double>& stat){
	if(paras_.size())
		throw std::logic_error("measurement::setMeasured: allowed to define only one measured distribution per instance");
	if(meas.size()!=stat.size())
		throw std::logic_error("measurement::setMeasured: measurement vector must be of same size as stat vector");
	isDifferential_=true;


	//remove exclude bin here!! //TBIExclude

	std::vector<TString> estnames=create_default_estnames(meas.size());
	std::vector<TString> used_estnames;
	for(size_t e=0;e<meas.size();e++){
		parameter pe;
		pe.setType(parameter::para_estimate);
		pe.setStat(stat.at(e));
		pe.setNominalVal(meas.at(e));
		TString name=estnames.at(e);
		pe.setName(name);
		paras_.push_back(pe);
		used_estnames.push_back(name);
	}
	est_corr_=correlationMatrix(used_estnames);

	c_external_=namedMatrix<uncertainty>(used_estnames,std::vector<TString>());
}

std::vector<TString> measurement::create_default_estnames(size_t nnames)const{

    std::vector<TString> estnames;
    for(size_t e=0;e<nnames;e++){
        TString name="__est_"+toTString(this_obj_counter_)+"_bin"+toTString(e);
        estnames.push_back(name);
    }
    return estnames;
}

void measurement::addSystematics(const TString& name, const TH1D* h){
	//check if there is overflow or underflow - for keeping track
	if(!est_corr_.size() && !H_.size())
		throw std::logic_error("measurement::addSystematics: first set measured values");

	int nbins=h->GetNbinsX()+1;
	int minbin=1,maxbin=nbins;

	std::vector<double> variedmeas;
	for(int i=minbin;i<maxbin;i++){
		variedmeas.push_back(h->GetBinContent(i));
	}
	addSystematics(name,variedmeas);
}
void measurement::addSystematics(const TString& name, std::vector<double> variedmeas){
    std::vector<uncertainty>  uncs(variedmeas.size());
    size_t i=0;
    for(const auto& p:paras_){
        if(p.getType()==parameter::para_estimate){
            variedmeas.at(i)-=p.getNominalVal();
            uncs.at(i)=uncertainty(variedmeas.at(i),-variedmeas.at(i));
            i++;
        }
    }
    addSystematics(name,uncs);
}

void measurement::addSystematics(const TString& name, double scaler){
    std::vector<uncertainty>  uncs;
    size_t i=0;
    for(const auto& p:paras_){
        if(p.getType()==parameter::para_estimate){
            double variedmeas = p.getNominalVal() * scaler - p.getNominalVal();
            uncs.push_back(uncertainty(variedmeas,-variedmeas));
            i++;
        }
    }
    addSystematics(name,uncs);
}


void measurement::addSystematics(const TString& name, std::vector<uncertainty> variedmeas){


	parameter newp;
	newp.setType(parameter::para_unc_absolute);
	newp.setName(name);
	paras_.push_back(newp);

	try{
		c_external_.addColumn(name, variedmeas);
	}catch(std::exception& e){
		std::cerr<< e.what()<<std::endl;
		throw std::out_of_range("in measurement::addSystematics");
	}
}

void measurement::setEstimateCorrelation(const size_t& i, const size_t& j, const double& val){
	if(H_.size())
		throw std::logic_error("measurement::setEstimateCorrelation: Correlation already defined by Hessian");
	est_corr_.setEntry(i,j,val);
}

void measurement::setHessian(const triangularMatrix&h){
	if(paras_.size() && !bypass_logic_check_)
		throw std::runtime_error("measurement::setHessian(const triangularMatrix&h): first set Hessian, then set additional parameters");
	setup_=false;
	TMatrixD H;
	h.toTMatrix(H);
	matrixHelper m(H);
	if(! m.checkMatrixPosDefinite(H))
		std::cout << "Warning: measurement::setHessian(const triangularMatrix&h): Hessian not positive definite. Fit might fail" << std::endl;;

	H_=h;

	if(paras_.size() && H_.size()<paras_.size()){
	    TString errormessage=(TString)"measurement::setHessian(const triangularMatrix&h): Hessian size "+
                H_.size()
                +(TString)" to small to describe all "+ paras_.size() +(TString)" parameters\n";
	    errormessage+="Possible reason: defined systematics before defining the estimate hessian/covariance leaves ambiguity w.r.t. the parameter sorting. Please define all estimates first and then add systematics (or start with the Hessian/Covariance).";
	    throw std::logic_error(errormessage.Data());
	}

	est_corr_.clear();
	if(!paras_.size()){
	    paras_.resize(H_.size());
	    for(size_t i=0;i<paras_.size();i++){
	        paras_.at(i).setName(H_.getEntryName(i));
	    }
	}
}

void measurement::setEstimateHessian(const TH2D&h){
    if(!paras_.size())
            throw std::runtime_error("measurement::setCovariance(const TH2&h): first set nominal input, then set esimate Hessian");

    triangularMatrix m(getEstimateNames());
    m.fillFromTH2(h);
    bypass_logic_check_=true;
    setHessian(m);
    bypass_logic_check_=false;
}

void measurement::setCovariance(const triangularMatrix&h){
    if(paras_.size() && !bypass_logic_check_)
        throw std::runtime_error("measurement::setCovariance(const triangularMatrix&h): first set Covariance, then set additional parameters");
    setup_=false;
    TMatrixD H;
    h.toTMatrix(H);
    matrixHelper m(H);
    if(! m.checkMatrixPosDefinite(H)){
        std::cout << "Warning: measurement::setCovariance(const triangularMatrix&h): Covariance not positive definite" << std::endl;
        std::cout << "Will try to construct Hessian nevertheless. But fit might fail. Check also debug output!" <<std::endl;
        triangularMatrix hes = h.createHessianFromCovariance();
        setHessian(hes);
    }
    else{
        triangularMatrix hcp=h;
        hcp.invert();
        setHessian(hcp);
    }

}

void measurement::setEstimateCovariance(const TH2D&h){
    if(!paras_.size())
        throw std::runtime_error("measurement::setCovariance(const TH2&h): first set nominal input, then set esimate Covariance");
    triangularMatrix m(getEstimateNames());
    m.fillFromTH2(h);
    bypass_logic_check_=true;
    setCovariance(m);
    bypass_logic_check_=false;
}

void measurement::setParameterType(const size_t & idx, parameter::para_type type){
	setup_=false;
	if(idx>=paras_.size()) throw std::out_of_range("measurement::setParameterType");
	if(type == parameter::para_unc_lognormal)
		throw std::runtime_error("measurement::setParameterType: log normal not fully supported yet");
	paras_.at(idx).setType(type);
}

void measurement::setParameterType(const TString & name, parameter::para_type type){
	setup_=false;
	//size_t idx=H_.getEntryIndex(name);
	//do not find the index based on the hessian in case it is not filled!
	for(auto& p: paras_){
		if(p.name() == name){
			p.setType(type);
			break;
		}
	}
}

void measurement::setParameterValue(const size_t & idx, const double& val){
	if(idx>=paras_.size()) throw std::out_of_range("measurement::setParameterValue");
	paras_.at(idx).setNominalVal(val);

}

void measurement::setParameterValue(const TString & name, const double& val){
	setup_=false;
	//size_t idx=H_.getEntryIndex(name);
	//do not find the index based on the hessian in case it is not filled!
	for(auto& p: paras_){
		if(p.name() == name){
			p.setNominalVal(val);
			break;
		}
	}
}





int measurement::getLeastSignificantBin()const{
    if(!isDifferential_)return -1;
    auto est = getEstimates();
    int lowestsign=-1;
    double temp=0;
    for(size_t i=0;i<est.size();i++){
        double nom = est.at(i).getNominalVal();
        double stat = 1/LM_.at(i).at(i);
        if(temp< stat/nom){
            temp=stat/nom;
            lowestsign=i;
        }
    }
    if(debug)
        std::cout << "measurement::getLeastSignificantBin: returned " << lowestsign <<std::endl;
    return lowestsign;

}
void measurement::setIsNormalisedInput(bool isn){
    std::cout << "measurement::setIsNormalisedInput: WARNING: experimental implementation" <<std::endl;
    isnormalisedinput_=isn;
}

bool measurement::isNormalisedInput()const{
    return isnormalisedinput_;
}

void measurement::setNormalise(bool norm){
    std::cout << "measurement::setIsNormalisedInput: WARNING: experimental implementation" <<std::endl;
    normalise_=norm;
}

////////// Interface to combiner ////////

void measurement::associateEstimate(const size_t & est_idx, const size_t & comb_idx){
	setup_=false;
	if(est_idx>=paras_.size()) throw std::out_of_range("measurement::associateEstimate");
	paras_.at(est_idx).setAsso(comb_idx);
	TString name=paras_.at(est_idx).name();
	for(auto& p:x_){
		if(p.name() == name){
			p.setAsso(comb_idx);
			break;
		}
	}
}

void measurement::associateEstimate(const TString & est_name, const size_t & comb_idx){

	for(auto& p:paras_){
		if(p.name() == est_name){
			p.setAsso(comb_idx);
			break;
		}
	}
	for(auto& p:x_){
		if(p.name() == est_name){
			p.setAsso(comb_idx);
			break;
		}
	}
}

void measurement::associateAllLambdas(const triangularMatrix& fullLambdaCorrs){
	if(debug)
		std::cout << "measurement::associateAllLambdas"  << std::endl;

	for(size_t i=0;i<fullLambdaCorrs.size();i++){
		const TString& name=fullLambdaCorrs.getEntryName(i);
		for(auto& p: paras_){
			if(name == p.name()){
				p.setAsso(i);
				if(debug)
					std::cout << p.name() << " " << i << std::endl;
			}
		}
		for(auto& p: lambdas_){
			if(name == p.name()){
				p.setAsso(i);
				if(debug)
					std::cout << p.name() << " " << i << std::endl;
			}
		}
	}
}


void measurement::setup(){
	if(debug)
		std::cout << "measurement::setup" <<std::endl;

	if(est_corr_.size() && H_.size())
		throw std::runtime_error("measurement::setup: hessian and direct estimate correlations provided. don't know what to chose");


	lambdas_.clear();
	x_.clear();
	M_.clear();
	kappa_.clear();
	tildeC_.clear();

	if(debug)
		std::cout << "measurement::setup: paras from Hessian" <<std::endl;

	//with Hessian paras first
	size_t nest=0,nlamb=0,nHest=0,nHlamb=0;
	for(const auto& p:paras_){
		if(p.getType() == parameter::para_estimate){
			nest++;
			if(H_.getEntryIndexUS(p.name())<SIZE_MAX){
				x_.push_back(p);
				nHest++;
			}
		}
		else{
			nlamb++;
			if(H_.getEntryIndexUS(p.name())<SIZE_MAX){
				lambdas_.push_back(p);
				nHlamb++;
			}
		}
	}


	M_.resize(nHest,std::vector<double>(nHest,0));
	kappa_.resize(nHest,std::vector<double>(nHlamb,0));
	tildeC_.resize(nHlamb,std::vector<double>(nHlamb,0));

	//order H in syst and estimates
	std::vector<TString> neworder;
	for(const auto& x:x_){
		neworder.push_back(x.name());
	}
	for(const auto& s:lambdas_){
		neworder.push_back(s.name());
	}

	if(debug)
		std::cout << "measurement::setup: reorder Hessian" <<std::endl; //DEBUG

	H_.reOrder(neworder);
	//fill Hessian parts
	for(size_t i=0;i<H_.size();i++){
		for(size_t j=0;j<H_.size();j++){
			if(i<nHest && j<nHest)
				M_.at(i).at(j)=H_.getEntry(i,j);
			if(i<nHest && j>=nHest)
				kappa_.at(i).at(j-nHest)=H_.getEntry(i,j);
			if(i>=nHest && j>=nHest)
				tildeC_.at(i-nHest).at(j-nHest)=H_.getEntry(i,j);
		}
	}


	if(debug)
		std::cout << "measurement::setup: add externalised" <<std::endl; //DEBUG
	//add rest
	for(const auto& p:paras_){
		if(p.getType() == parameter::para_estimate){
			if(H_.getEntryIndexUS(p.name())==SIZE_MAX)
				x_.push_back(p);
		}
		else{
			if(H_.getEntryIndexUS(p.name())==SIZE_MAX)
				lambdas_.push_back(p);
		}
	}





	//new approach. find paras in Hessian
	LM_.resize(nest,std::vector<double>(nest,0));
	Lk_.resize(nest,std::vector<uncertainty>(nlamb));
	LD_.resize(nlamb,std::vector<double>(nlamb,0));


	if(debug)
		std::cout << "measurement::setup: fill LM" <<std::endl; //DEBUG

	//fill LM_
	for(size_t i=0;i<x_.size();i++){
		for(size_t j=0;j<x_.size();j++){
			if(H_.size()){
				LM_.at(i).at(j)=H_.getEntry(i,j);
			}
			else{
				if(est_corr_.size()!=x_.size())
					throw std::runtime_error("internal error! est corr size");
				LM_.at(i).at(j)=est_corr_.getEntry(i,j)/(x_.at(i).stat()*x_.at(j).stat());
			}
		}
	}



	if(debug)
		std::cout << "measurement::setup: calculate Lk_" <<std::endl; //DEBUG
	//calculate for the inputs from Hessians
	/*
	 * *** calculate ks from initial Hessian part kappa
	 * *** if one is rel, leave abs 0 and vice versa
	 */
	TMatrixT< double> TM(M_.size(),M_.size());
	for(size_t i=0;i<M_.size();i++){
		for(size_t j=0;j<M_.size();j++){
			TM[i][j]=M_.at(i).at(j);
		}
	}
	double det=0;
	TM.Invert(&det);
	if(!det){
		throw std::runtime_error("Part of Hessian for estimates not invertible! Serious error, please check your input");
	}

	size_t kssize=0;
	if(kappa_.size())
		kssize=kappa_.at(0).size();
	// k_i = M^-1 kappa_i
	for(size_t i=0;i<kssize;i++){
		std::vector< double> thiskappas(M_.size());
		for(size_t mu=0;mu<M_.size();mu++){
			thiskappas.at(mu) = kappa_.at(mu).at(i);
		}
		TVectorT< double> vec(thiskappas.size(), &thiskappas.at(0));
		TVectorT< double> ki= TM * vec;
		for(size_t mu=0;mu<M_.size();mu++){
			Lk_.at(mu).at(i)=uncertainty(-ki[mu],ki[mu]);
		}
	}


	if(debug)
		std::cout << "measurement::setup: fill external Lk" <<std::endl; //DEBUG

	//fill the rest of Lks // LD remains 0 for them
	for(size_t i=kssize;i<nlamb;i++){ //USE THE MATRIX HERE> >>>> >>>>
		size_t exidxi= c_external_.getEntryIndexUSCol(lambdas_.at(i).name());
		for(size_t mu=0;mu<nest;mu++){
			size_t exidxmu = c_external_.getEntryIndexUSRow(x_.at(mu).name());
			Lk_.at(mu).at(i)=c_external_.getEntry(exidxmu,exidxi);
		}
	}


	if(debug)
		std::cout << "measurement::setup: create LD" <<std::endl; //DEBUG

	//fill only the first rows and cos of LD_. rest remains 0
	//all Lk_ here are symmetric by construction
	for(size_t i=0;i<tildeC_.size();i++){
		for(size_t j=0;j<tildeC_.size();j++){
			long double estsum_ij=0;
			for(size_t mu=0;mu<x_.size();mu++){
				for(size_t nu=0;nu<x_.size();nu++){
					estsum_ij +=  (long double)((long double)M_.at(mu).at(nu)/2.) * //+kr_delta(mu,nu))) *
							(long double)(Lk_.at(mu).at(i).symm()*Lk_.at(nu).at(j).symm() + Lk_.at(nu).at(i).symm()*Lk_.at(mu).at(j).symm());

				}
			}
			LD_.at(i).at(j)=  tildeC_.at(i).at(j) - estsum_ij - kr_delta(i,j) ;//fixed
		}
	}

	//make a check here
	for(size_t i=0;i<tildeC_.size();i++){
		if(LD_.at(i).at(i) +1 < 0){
			std::cout << LD_ << std::endl;
			throw std::runtime_error("LD_ diag +1 < 0");
		}
	}


	if(0){
		//for testing with neglecting correlations in a quick and dirty way.
		triangularMatrix Ltriag(LD_.size());
		for(size_t i=0;i<LD_.size();i++)
			for(size_t j=i;j<LD_.size();j++)
				Ltriag.setEntry(i,j,LD_.at(i).at(j) + kr_delta(i,j));

		double det;
		Ltriag.invert(det);
		if(!det)
			throw std::runtime_error("det 0");

		for(size_t i=0;i<LD_.size();i++)
			for(size_t j=i;j<LD_.size();j++){
				if(i!=j)
					Ltriag.setEntry(i,j,0);
			}
		Ltriag.invert(det);
		if(!det)
			throw std::runtime_error("det 0");

		for(size_t i=0;i<LD_.size();i++)
			for(size_t j=0;j<LD_.size();j++)
				LD_.at(i).at(j)=Ltriag.getEntry(i,j)-kr_delta(i,j);
	}

	triangularMatrix hinv=H_;
	hinv.invert();


	normalisation_=0;
	for(const auto& p:x_){
	    normalisation_+=p.getNominalVal();
	}
    //debug

	//DEBUG
	if(debug){
		std::cout << "H:\n" << H_<<std::endl;
		std::cout << "M: " << M_<<std::endl;
		std::cout << "kappa: " << kappa_<<std::endl;
		std::cout << "C_tilde: " << tildeC_<<'\n'<<std::endl;

		std::cout << "LM: " << LM_ << std::endl;
		std::cout << "k_tilde: " << Lk_<<std::endl;
		std::cout << "D: " << LD_<<std::endl;
	}

	setup_=true;
}



double measurement::evaluate(const double* pars, double* df, const bool& pearson, const size_t& maxidx)const{
	if(!setup_)
		throw std::runtime_error("measurement::evaluate: first call setup()");
	//cannot be done automatically since this function needs to be const

	const size_t nlamb=lambdas_.size();
	const size_t nx=x_.size();


	//std::cout << "eval" <<std::endl;


	// TBI implement gradient df! then, also the MP might become
	// reasonable. For now it does not increase perf.
	//df all set to zero before, only ADD here

	long double xsqpart=0;
	//x^2 part - symmetric
	double combsum=0;
	for(size_t mu=0;mu<nx;mu++){
	    combsum += pars[x_.at(mu).getAsso()];
	}
	//double rest_norm_comb = 1 - combsum;

	if(normalise_){
	    //calculate addition to scaler?... TBI
	    //and apply scaler to x_comb_mu/nu
	}

	for(size_t mu=0;mu<nx;mu++){


		double x_comb_mu = pars[x_.at(mu).getAsso()];
		/*
		if(excludebin_ == (int) mu){
		    if(isnormalisedinput_)
		        continue;
		    x_comb_mu=rest_norm_comb;
		}
        if(normalise_)
            x_comb_mu*=normalisation_;
		 */
		double x_meas_mu = x_.at(mu).getNominalVal();

		for(size_t i=0;i<nlamb;i++){
			if(! lambdas_.at(i).isRelative()) continue;
			double lambda_i = pars[lambdas_.at(i).getAsso()];
			x_comb_mu *= (1 - Lk_.at(mu).at(i).eval(lambda_i)/x_meas_mu);
		}
		for(size_t i=0;i<nlamb;i++){
			if(lambdas_.at(i).isRelative()) continue;
			double lambda_i = pars[lambdas_.at(i).getAsso()];
			x_comb_mu -= Lk_.at(mu).at(i).eval(lambda_i);
		}

		for(size_t nu=mu;nu<nx;nu++){

			double x_comb_nu = pars[x_.at(nu).getAsso()];
			double x_meas_nu = x_.at(nu).getNominalVal();

			/*
	        if(excludebin_ == (int) nu){
	            if(isnormalisedinput_)
	                continue;
	            x_comb_nu=rest_norm_comb;
	        }
	        if(normalise_)
	            x_comb_nu*=normalisation_;
			 */
			for(size_t i=0;i<nlamb;i++){
				if(! lambdas_.at(i).isRelative()) continue;
				double lambda_i = pars[lambdas_.at(i).getAsso()];
				x_comb_nu *= (1 - Lk_.at(nu).at(i).eval(lambda_i)/x_meas_nu);
			}
			for(size_t i=0;i<nlamb;i++){
				if(lambdas_.at(i).isRelative()) continue;
				double lambda_i = pars[lambdas_.at(i).getAsso()];
				x_comb_nu -= Lk_.at(nu).at(i).eval(lambda_i);
			}

			long double contribution = (x_meas_mu - x_comb_mu) * LM_.at(mu).at(nu) * (x_meas_nu - x_comb_nu);
			if(pearson){
				if(mu==nu){
					contribution*= (long double)x_meas_nu/x_comb_nu;
				}
				else{
					contribution*= (long double)std::sqrt(x_meas_nu/x_comb_nu * x_meas_mu/x_comb_mu);
				}
			}


			if(mu!=nu)
				contribution*=(long double)2.;
			xsqpart+=contribution;
		}
	}


	long double copart=0;
	//lambda^2 - symmetric in lambda
	for(size_t i=0;i<nlamb;i++){
		double lambda_i = pars[lambdas_.at(i).getAsso()];
		for(size_t j=i;j<nlamb;j++){


			double lambda_j = pars[lambdas_.at(j).getAsso()];

			long double contribution=lambda_i * lambda_j * LD_.at(i).at(j);

			if(i!=j)
				contribution*=2.;
			copart+=contribution;
		}
	}


	if(xsqpart != xsqpart  || copart!=copart){
		for(size_t i=0;i<lambdas_.size();i++){
			std::cout << lambdas_.at(i).name() <<" "<< pars[lambdas_.at(i).getAsso()] << std::endl;
		}
		for(size_t i=0;i<x_.size();i++){
			std::cout << x_.at(i).name() <<" "<< pars[x_.at(i).getAsso()]<< std::endl;
		}
		if(xsqpart != xsqpart )
			throw std::runtime_error("measurement::evaluate: nan in x^2 part");
		if(copart!=copart)
			throw std::runtime_error("measurement::evaluate: nan in lambda^2 part");
	}



	double out=xsqpart + copart;

	return out;
}

double measurement::getCombSum(const double * pars)const{
    const size_t nx=x_.size();
    double combsum=0;
    for(size_t mu=0;mu<nx;mu++){
        combsum += pars[x_.at(mu).getAsso()];
    }
    return combsum;
}



/////////// privates //////////

void measurement::copyFrom(const measurement& r){
	setup_=r.setup_;
	H_=r.H_;
	paras_=r.paras_;
	M_=r.M_;
	kappa_=r.kappa_;
	tildeC_=r.tildeC_;
	LM_=r.LM_;
	Lk_=r.Lk_;
	LD_=r.LD_;
	lambdas_=r.lambdas_;
	x_=r.x_;
	c_external_=r.c_external_;
	est_corr_=r.est_corr_;
	isDifferential_=r.isDifferential_;
}

std::vector<TString> measurement::getParameterNames()const{
    std::vector<TString> out;
    for(const auto& p: paras_)
        out.push_back(p.name());
    return out;
}
std::vector<TString> measurement::getEstimateNames()const{
    std::vector<TString> out;
    for(const auto& p: paras_){
        if(p.getType() != parameter::para_estimate) continue;
        out.push_back(p.name());
    }
    return out;

}



const parameter& measurement::getParameter(const TString& name)const{
	for(const auto& p:paras_)
		if(p.name() == name)return p;

	throw std::out_of_range(("measurement::getParameter "+name+" not found").Data());
}
