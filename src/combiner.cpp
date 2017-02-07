/*
 * combiner.cc
 *
 *  Created on: 30 Jun 2016
 *      Author: jkiesele
 */



#include "combiner.h"
#include "simpleFitter.h"
#include "fileReader.h"
#include "helpers.h"
#include "TFile.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TRandom3.h"
#include "matrixHelper.h"
#include "TCanvas.h"
#include "graphCosmetics.h"
#include "TStyle.h"
#include "fitfunction.h"

#ifdef USE_MP
#include <omp.h>
#endif


/*
 * Notes:
 *
 *   provide derivatives to minuit. should be simple for all abs unc. and neyman. otherwise more complicated
 *   use ROOT::Math::IMultiGradFunction inheritance
 *
 */


//small helper class
class scanresults{
public:
	void addResult(const combinationResult& r, float coeff){
		res_.push_back(r);
		coeffs_.push_back(coeff);
	}
	void sort(){
		if(coeffs_.size()<1)return;

		bool notdone=true;
		std::vector<combinationResult> sortedres;
		std::vector<float> sortedcoeffs;
		std::vector<size_t> used;
		size_t smallestidx=0;
		while(notdone){
			double smallest=1.1;
			for(size_t i=0;i<coeffs_.size();i++){
				if(std::find(used.begin(),used.end(),i)!=used.end())continue;
				if(coeffs_.at(i)<smallest){
					smallest=coeffs_.at(i);
					smallestidx=i;
				}
			}
			sortedres.push_back(res_.at(smallestidx));
			sortedcoeffs.push_back(coeffs_.at(smallestidx));
			used.push_back(smallestidx);
			if(used.size() == res_.size()) notdone=false;
		}
		res_=sortedres;
		coeffs_=sortedcoeffs;
	}
	void clear(){
		res_.clear();
		coeffs_.clear();
	}
	combinationResult& at(size_t i){return res_.at(i);}
	const combinationResult& at(size_t i)const{return res_.at(i);}

	float& C_at(size_t i){return coeffs_.at(i);}
	const float& C_at(size_t i)const{return coeffs_.at(i);}

	size_t size()const{return res_.size();}

	const std::vector<combinationResult>& getVector()const{return res_;}
private:
	std::vector<combinationResult> res_;
	std::vector<float> coeffs_;
};


bool combiner::debug=false;
bool combiner::dummyrun_=false;
const double combiner::maxcorr_=0.999;




void combiner::addMeasurement(const std::string& infile){
	if(debug)
		std::cout << "combiner::addMeasurement" <<std::endl;
	addMeasurement(measurement(infile));
}
void combiner::printCorrelationMatrix()const{
	if(debug)
		std::cout << "combiner::printCorrelationMatrix" <<std::endl;
}
void combiner::readConfigFile(const std::string & filename){
	isdifferential_=false;


	configfile_=filename;
	fileReader fr;
	fr.setComment("#");
	fr.setDelimiter(",");
	fr.setStartMarker("[inputs]");
	fr.setEndMarker("[end inputs]");
	fr.readFile(configfile_);
	int nmeas=fr.getValue<int>("nFiles");

	std::string configfiledir=textFormatter::getFileDir(configfile_);
	if(configfiledir.length()>0)
		configfiledir+="/";


	for(int i=0;i<nmeas;i++){
		std::string measfile=configfiledir+fr.getValue<std::string>("file"+toString(i));
		addMeasurement(measfile);
		//constraints

	}
	//std::cout << "all read in" <<std::endl; //DEBUG
	fr.clear();//will be used later

	fr.setStartMarker("[observables]");
	fr.setEndMarker("[end observables]");
	fr.setDelimiter("=");
	fr.readFile(configfile_);
	for(size_t i=0;i<fr.nLines();i++){
		if(fr.nEntries(i)<2)continue;
		TString merged=fr.getData<TString>(i,0);
		std::string tomerge=fr.getData<std::string>(i,1);
		textFormatter tf;
		tf.setDelimiter("+");
		std::vector<std::string> ins=tf.getFormatted(tomerge);
		for(size_t j=0;j<ins.size();j++){
			associatePriv(ins.at(j).data(),merged);
		}
	}
	///needs to be called to set up external correlation matrix


	fr.setStartMarker("[correlations]");
	fr.setEndMarker("[end correlations]");
	fr.setDelimiter("=");
	fr.readFile(configfile_);
	for(size_t i=0;i<fr.nLines();i++){
		if(fr.nEntries(i)<2)continue;
		correlationscan cscan;
		TString sysname=fr.getData<TString>(i,0);
		size_t idx=0;
		try{
			idx=external_correlations_.getEntryIndex(sysname);
		}catch(...){
			throw std::runtime_error(("combiner::setup: reading in correlation assumptions: failed to find systematic "+sysname).Data());
		}
		cscan.idxa=idx;
		cscan.nominal=0;
		cscan.high=0;
		cscan.low=0;
		cscan.steps=0;

		std::string corrinfo=fr.getData<std::string>(i,1);
		textFormatter tf;
		tf.setDelimiter("+");
		std::vector<std::string> contributions=tf.getFormatted(corrinfo);
		for(size_t j=0;j<contributions.size();j++){
			tf.setDelimiter(")");
			tf.setTrim(" (");
			std::vector<std::string> entry=tf.getFormatted(contributions.at(j));
			if(entry.size() < 1)continue;
			try{
				idx=external_correlations_.getEntryIndex(entry.at(1).data());
			}catch(...){
				throw std::runtime_error(("combiner::setup: reading in correlation assumptions: failed to find systematic "+(TString)entry.at(1)).Data());
			}
			cscan.idxb=idx;
			tf.setDelimiter("&");
			tf.setTrim(" ");
			std::vector<std::string> nominalandrange=tf.getFormatted(entry.at(0));
			if(nominalandrange.size()>0){
				cscan.nominal=atof(nominalandrange.at(0).data());
				cscan.high=cscan.nominal;
				cscan.low=cscan.nominal;
			}
			if(nominalandrange.size()>1){
				tf.setDelimiter(":");
				std::vector<std::string> range=tf.getFormatted(nominalandrange.at(1));
				if(range.size()<2)
					throw std::runtime_error(("combiner::setup: reading in correlation assumptions: failed to find scan range for "+sysname+": "+ entry.at(1)).Data());
				cscan.low=atof(range.at(0).data());
				cscan.high=atof(range.at(1).data());

				cscan.steps=4;// fixed. suffices
			}
			syst_scanranges_.push_back(cscan);
		}
	}
	//fill nominal assumptions
	for(size_t i=0;i<syst_scanranges_.size();i++){
		external_correlations_.setEntry(syst_scanranges_.at(i).idxa,syst_scanranges_.at(i).idxb, syst_scanranges_.at(i).nominal);
	}


}

void combiner::addMeasurement( measurement m){
	m.setup();
	if(measurements_.size()<1){
		isdifferential_=m.isDifferential();
		if(isdifferential_){
			hasUF_=m.hasUF();
			hasOF_=m.hasOF();
		}
	}
	else{
		if(isdifferential_ != m.isDifferential()
				|| (isdifferential_&& hasUF_!=m.hasUF())
				|| (isdifferential_&& hasOF_!=m.hasOF()))
			throw std::logic_error("combiner::addMeasurement: cannot combine measurements using root interface with measurements not using the root interface");
	}
	//check for ambiguities
	for(const auto& p: allparas){
		for(const auto& p1: m.getParameters()){
			if(p.name() == p1.name())
				throw std::logic_error("combiner::addMeasurement: uncertainties and estimates must have uniquie naming: "+(std::string)p.name().Data());
		}
	}

	measurements_.push_back(m);
	allparas.insert(allparas.end(),m.getParameters().begin(),m.getParameters().end());
	allsysparas_.insert(allsysparas_.end(), m.getLambdas().begin(), m.getLambdas().end());
	createExternalCorrelations();
}


void combiner::createExternalCorrelations(){
	//get associations
	std::vector<TString> names;
	for(const auto& p: allsysparas_)
		names.push_back(p.name());

	external_correlations_  = correlationMatrix(names);
	external_correlations_.setOffDiagonalZero();

	if(debug)
		std::cout << "setup done \n"<< external_correlations_<<std::endl;
}
void combiner::associatePriv(const TString & a, const TString& outname){
	if(debug)
		std::cout << "combiner::associate "<< a<<", "<<outname <<std::endl;

	//check
	for(size_t i=0;i<allparas.size()+1;i++){
		if(i==allparas.size())
			throw std::out_of_range(("combiner::associate "+a+" not found").Data());
		if(allparas.at(i).name() == a){
			break;
		}
	}


	for(auto& c: tobecombined_){
		if(c.first == outname){
			c.second.push_back(a);
			return;
		}
	}
	tobecombined_.push_back(std::pair< TString, std::vector<TString> > (outname, std::vector<TString> (1,a)));

}
void combiner::associate(const TString & a, const TString& outname){
	associatePriv(a,outname);
}



std::vector<std::vector<combinationResult> > combiner::scanCorrelationsIndep(std::ostream& out, const combinationResult& nominal, const std::string& outdir)const{

	/*
	 * First part:
	 *   Do the scan and save output to temp
	 */
	const size_t nscans=syst_scanranges_.size();
	if(nscans<1)
		return std::vector<std::vector<combinationResult> > ();



	//fully serialise for parallelisation
	std::vector<const correlationscan* > scanps;
	std::vector<double> scanvs;
	for(size_t ss=0;ss<nscans;ss++){
		const correlationscan* s= &syst_scanranges_.at(ss);
		double range=s->high-s->low;
		for(size_t i=0;i<=s->steps;i++){
			double value=s->low+ (double)i*range/(double)s->steps;
			scanps.push_back(s);
			scanvs.push_back(value);
		}
	}

	std::vector<combinationResult> ress(scanps.size());

	int ndone=0;
#ifdef USE_MP
#pragma omp parallel for
#endif
	for(size_t i=0;i<scanps.size();i++){
		bool success=true;
		combinationResult result;
		try{
			combiner combcp=*this;
			combcp.setSystCorrelation(scanps.at(i)->idxa,scanps.at(i)->idxb,scanvs.at(i));
			result=combcp.combine();
			if(dummyrun_)
				//sleep(5); //for MP debugging
			ndone++;
		}catch(std::exception& e){
#ifdef USE_MP
#pragma omp critical (cout)
#endif
			{
				std::cout << e.what() << std::endl;
				std::cout <<  external_correlations_.getEntryName(scanps.at(i)->idxa) << ", " << external_correlations_.getEntryName(scanps.at(i)->idxb) <<std::endl;
				std::cout << "scanning for one point (" << scanvs.at(i) << ") failed. Please check the scan plots carefully." ;
				success=false;
				ndone++;
			}
		}
		if(success){
#ifdef USE_MP
#pragma omp critical (scanpointresult)
#endif
			{
				ress.at(i)=result;
			}
#ifdef USE_MP
#pragma omp critical (cout)
#endif
			{
				std::cout << "\nscan progress: " << (int)((float)ndone/(float)scanvs.size()*100)<<"%" <<std::endl;
			}
		}
	}

	std::vector<scanresults >results;
	//order
	for(size_t ss=0;ss<nscans;ss++){
		const correlationscan* s= &syst_scanranges_.at(ss);
		scanresults resforthisscan;
		for(size_t i=0;i<scanps.size();i++){
			if(scanps.at(i) == s){
				if(ress.at(i).getNCombined()) //combination has not failed
					resforthisscan.addResult(ress.at(i),scanvs.at(i));
			}
		}
		resforthisscan.sort();
		out << "\nscanned " <<  external_correlations_.getEntryName(s->idxa) <<
				" || " << external_correlations_.getEntryName(s->idxb)  << "...\n"<<std::endl;
		for(size_t i=0;i<resforthisscan.size();i++){
			out << "rho = "<< resforthisscan.C_at(i) <<std::endl;
			resforthisscan.at(i).printResultOnly(out);
		}
		results.push_back(resforthisscan);
	}
	ress.clear();//free some mem before plotting


	//old

	const size_t nobs=tobecombined_.size();


	/*
	 * Second part:
	 *   Write output as graphs to TFile
	 */
	TFile * tfile=new TFile((TString)outprefix_+"scanPlots.root","RECREATE");

	for(size_t i=0;i<results.size();i++){
		const correlationscan * s= & syst_scanranges_.at(i);
		const TString & scanneda= external_correlations_.getEntryName(s->idxa);
		const TString & scannedb= external_correlations_.getEntryName(s->idxb);
		TString scannedacomp=textFormatter::makeCompatibleFileName(scanneda.Data() );
		TString scannedbcomp=textFormatter::makeCompatibleFileName(scannedb.Data() );
		TString scanned="#rho("+scanneda+", "+scannedb+")";
		for(size_t obs=0;obs<nobs;obs++){
			/*
			 * Prepare graph content input
			 */
			TString graphname=textFormatter::makeCompatibleFileName(( tobecombined_.at(obs).first+"_"+scanned).Data());
			std::vector<double> scan,zero,comb,errup,errdown;
			std::vector<double> minchi2;
			for(size_t j=0;j<results.at(i).size();j++){
				const combinationResult& result=results.at(i).at(j);

				if(result.getNCombined() != nobs){//sanity check for debugging
					TString errstr=nominal.getCombNames().at(obs)+"_"+scannedacomp+"_"+scannedbcomp ;
					throw std::logic_error("combiner::scanCorrelations: fatal error. not all (or too many) combined results: "+(std::string)errstr.Data());
				}
				comb.push_back(result.combined_.at(obs));
				errup.push_back(result.comberrup_.at(obs));
				errdown.push_back(-result.comberrdown_.at(obs));
				minchi2.push_back(result.getChi2min());
				scan.push_back(results.at(i).C_at(j));
				zero.push_back(0);
			}
			/*
			 * Create graph and add some cosmetics
			 */
			TGraphAsymmErrors * g=new TGraphAsymmErrors(comb.size(), &scan.at(0), &comb.at(0),
					&zero.at(0), &zero.at(0), &errup.at(0), &errdown.at(0));

			g->GetYaxis()->SetTitle(tobecombined_.at(obs).first);
			applyGraphCosmetics(g,gc_scancombined,s->low,s->high,graphname,scanned);
			g->SetName(graphname);
			g->SetTitle(graphname);
			g->Write();

			TGraphAsymmErrors * gminchi2=new TGraphAsymmErrors(comb.size(), &scan.at(0), &minchi2.at(0),
					&zero.at(0), &zero.at(0), &zero.at(0), &zero.at(0));
			applyGraphCosmetics(gminchi2,gc_minchi2,s->low,s->high, "min#chi^2",scanned);
			gminchi2->SetName("min#chi^{2}");
			gminchi2->SetTitle("min#chi^{2}");
			gminchi2->GetYaxis()->SetTitle("min#chi^{2}");
			if(!obs)
				gminchi2->Write();



			if(outdir.length()){
				system(("mkdir -p "+outdir).data());
				TCanvas cv;
				cv.SetLeftMargin(.15);
				cv.SetBottomMargin(.15);

				gStyle->SetOptTitle(0);

				double nomcorr= nominal.getInputSysCorrelations().getEntry(s->idxa,s->idxb);
				double nomerrd=-nominal.getCombErrdown().at(obs);
				TGraphAsymmErrors gnom(1, &nomcorr, &nominal.getCombined().at(obs),
						&zero.at(0), &zero.at(0), &nominal.getCombErrup().at(obs), &nomerrd);

				applyGraphCosmetics(&gnom,gc_nominal,s->low,s->high,graphname+"_nom",scanned);
				cv.Draw();
				g->Draw("Aa3pl");
				gnom.Draw("Pe");

				cv.Print((TString)outdir+"/"+nominal.getCombNames().at(obs)+"_"+scannedacomp+"_"+scannedbcomp +".pdf");

				cv.Clear();
				applyGraphCosmetics(g,gc_scancombinedUP,s->low,s->high,graphname,scanned,2.05);
				cv.Divide(1,2);
				TVirtualPad * pad=cv.cd(1);
				pad->SetBottomMargin(0.015);
				pad->SetLeftMargin(.15);
				pad->SetTopMargin(.189);
				g->Draw("Aa3pl");
				gnom.Draw("Pe");

				applyGraphCosmetics(gminchi2,gc_minchi2,s->low,s->high, "min#chi^2",scanned,2.05);
				pad=cv.cd(2);
				pad->SetBottomMargin(0.293);
				pad->SetTopMargin(0.025);
				pad->SetLeftMargin(.15);
				gminchi2->Draw("Al");



				cv.Print((TString)outdir+"/"+nominal.getCombNames().at(obs)+"_"+scannedacomp+"_"+scannedbcomp +"_detailed.pdf");

			}


			// g is owned by TFile, no delete necessary
		}
	}

	tfile->Close();
	delete tfile;
	std::vector<std::vector<combinationResult> > vout;
	for(const auto& r:results)
		vout.push_back(r.getVector());
	return vout;
}


combinationResult combiner::combine()const{
	combiner cp=*this;
	return cp.combinePriv();
}

combinationResult combiner::combinePriv(){
	if(debug)
		std::cout << "combiner::combine" <<std::endl;

	combinationResult out;
	/*
	 * for differential distributions associate automatically here
	 */
	if(tobecombined_.size()<1 && isdifferential_){ //set up assos automatically
		size_t psize=0;
		for(const auto& m: measurements_){
			const std::vector<parameter>& ps=m.getEstimates();
			if(!psize)psize=ps.size();
			if(psize!=ps.size())
				throw std::runtime_error("combiner::combinePriv: trying to combine differential distributions with a different number of bins. Check for overflow/underflow");
			for(size_t i=0;i<ps.size();i++)
				associate(ps.at(i).name(),"comb_"+toTString(i));
		}
	}



	const size_t nsys=external_correlations_.size();

	std::vector<double>fitparas(nsys, 0);


	for(size_t i=0;i<tobecombined_.size();i++){
		double mean=0;
		for(auto& m:measurements_){
			for(size_t j=0;j<tobecombined_.at(i).second.size();j++){
				double thisval=0;
				try{
					thisval=m.getParameter(tobecombined_.at(i).second.at(j)).getNominalVal();
					m.associateEstimate(tobecombined_.at(i).second.at(j),i+nsys);
				}catch(...){}
				if(thisval)mean+=thisval;
			}
		}
		mean/=(double)(tobecombined_.at(i).second.size());

		fitparas.push_back(mean);
	}

	for(auto& m:measurements_){
		m.associateAllLambdas(external_correlations_); //get lamda indecies right
		//m.setup();
	}


	npars_=fitparas.size();
	nest_=tobecombined_.size();
	std::vector<double> steps(fitparas.size(),1);//remove in the end


	std::vector<TString> names;
	for(size_t i=0;i<nsys;i++)
		names.push_back(external_correlations_.getEntryName(i));
	for(size_t i=0;i<tobecombined_.size();i++){
		names.push_back(tobecombined_.at(i).first);

	}

	external_correlations_.toTMatrix(inv_priors_);
	double det=0;
	inv_priors_.Invert(&det);
	if(!det || det<0){
		throw std::runtime_error("combiner::combinePriv: external correlations non invertible or not positive definite");
	}


	/*
	 *********   Configure the fitter and do the fit
	 */

	simpleFitter fitter;
	fitter.setOnlyRunDummy(dummyrun_);
	fitter.setParameters(fitparas,steps);
	for(size_t i=0;i<fitparas.size();i++){
		if(i>=nsys){
			fitter.setAsMinosParameter(i,true); //use minos for all quantities to be combined
			//	fitter.setParameterLowerLimit(i,0); //let in for now
		}
	}
	fitter.setParameterNames(names);
	fitfunctionGradient func(this);
	fitter.setMinFunction(func);

	if(debug)
		simpleFitter::printlevel=2;

	/*
	 * Make a first rough fit to get better starting values.
	 * Can help in case of convergence problems, usually not needed
	 */

	//throw std::runtime_error("abort"); //DEBUG
	/*
	 * Configure Minuit parameters and fit
	 */
	fitter.setStrategy(2);
	fitter.setTolerance(0.1);
	fitter.fit();
	if(!fitter.wasSuccess()){
		throw std::runtime_error("combiner::combine: fit not successful");
	}

	/*
	 *********  Save the results
	 */

	const size_t nsyst=nsys;
	const size_t ncomb=fitter.getParameters()->size()-nsyst ;
	std::vector<TString> combnames;
	for(size_t i=0;i<ncomb;i++){
		size_t idx=nsyst+i;
		combnames.push_back(fitter.getParameterNames()->at(idx));
	}

	out.orig_sys_correlations_=external_correlations_;
	//	out.meas_correlations_=estimate_correlations_;
	//	out.post_meas_correlations_=estimate_correlations_;
	//	out.post_sys_correlations_=syst_correlations_;
	out.chi2min_=fitter.getChi2Min();

	if(debug)
		std::cout << "combined (minL="<<fitter.getChi2Min()<<"):"<<std::endl;



	out.post_meas_correlations_=correlationMatrix(combnames);
	for(size_t i=0;i<ncomb;i++){
		size_t idx=nsyst+i;
		if(debug){
			std::cout << fitter.getParameterNames()->at(idx) << ": " << fitter.getParameter(idx)<<
					" +" <<  fitter.getParameterErrUp()->at(idx) << " "
					<<  fitter.getParameterErrDown()->at(idx) <<std::endl;
		}
		out.combined_.push_back(fitter.getParameter(idx));
		out.comberrup_.push_back(fitter.getParameterErrUp()->at(idx));
		out.comberrdown_.push_back(fitter.getParameterErrDown()->at(idx));
		out.combnames_.push_back(fitter.getParameterNames()->at(idx));
		for(size_t j=0;j<ncomb;j++){
			size_t idxb=nsyst+j;
			out.post_meas_correlations_.setEntry(i,j, fitter.getCorrelationCoefficient(idx,idxb));
		}
	}

	out.post_sys_correlations_=external_correlations_;
	for(size_t i=0;i<nsyst;i++){
		out.pulls_.push_back(fitter.getParameter(i));
		out.constraints_.push_back(fitter.getParameterErr(i));
		for(size_t j=0;j<=i;j++){
			out.post_sys_correlations_.setEntry(i,j,fitter.getCorrelationCoefficient(i,j));
		}
	}

	if(debug && nsyst<20)
		std::cout << "post combine systematics correlations: \n"<< out.post_sys_correlations_ << std::endl;



	if(debug){
		std::cout << "correlations between combined quantities" <<std::endl;
		std::cout << out.post_meas_correlations_ << std::endl;
	}

	/*
	 * Save the weights of each estimate for each combined value
	 */

	/* DEBUG FIXME
	out.estimateweights_.resize(ncomb);
	for(size_t i=0;i<associations_.size();i++){
		associations_.at(i).first ;
		double combined=out.combined_.at(i);
		for(size_t j=0;j<associations_.at(i).second.size();j++){
			double inputest=estimates_.at( associations_.at(i).second.at(j) ).getNominalVal();
		}
	}

	for one weight.. for more - not straightforward?
	c = w1 a + (1-w1) b
			c - b = w1 (a - b)
	 */

	/*
	 * **** check for deviations of parameters that were coupled with maximum correlation < 1
	 */
	for(size_t i=0;i<out.orig_sys_correlations_.size();i++){
		for(size_t j=0;j<out.orig_sys_correlations_.size();j++){
			const double& preCij=out.orig_sys_correlations_.getEntry(i,j);
			if(preCij == maxcorr_ || preCij == -maxcorr_){
				double sign = 1;
				if(preCij < 0) sign = -1;

				if( fabs(out.pulls_.at(i) - sign * out.pulls_.at(j)) > (1-maxcorr_)*30)
					std::cout << "WARNING: the parameters " <<  out.orig_sys_correlations_.getEntryName(i)
					<< " and " << out.orig_sys_correlations_.getEntryName(j)
					<< " have been set to be fully correlated. However, they show a difference"
					<< " of more than 3 per-mille (" <<out.pulls_.at(i) - sign * out.pulls_.at(j)  << ") after the combination. If these parameters have"
					<< " a non-negligible impact and the difference is large, please take the result with care." <<std::endl;
			}
		}
	}
	out.hasUF_=hasUF_;
	out.hasOF_=hasOF_;

	clear(); //to be sure this instance is not used anymore
	return out;
}




void combiner::setSystCorrelation(const TString & namea, const TString& nameb, const double& coeff){
	size_t idxa=external_correlations_.getEntryIndex(namea);
	size_t idxb=external_correlations_.getEntryIndex(nameb);
	setSystCorrelation(idxa,idxb,coeff);
}
void combiner::setSystCorrelation(const size_t & idxa, const size_t& idxb, const double& coeff){
	double usecoeff=coeff;
	if(usecoeff>maxcorr_)usecoeff=maxcorr_;
	if(usecoeff<-maxcorr_)usecoeff=-maxcorr_;
	external_correlations_.setEntry(idxa,idxb, usecoeff);
}





void combiner::clear(){
	external_correlations_.clear();
	tobecombined_.clear();
	measurements_.clear();
}

