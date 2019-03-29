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
#include <iomanip>

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
    auto m = measurement(infile);
    m.setIsDifferential(isdifferential_);
    if(norm_constraint_>0)
        m.setNormalisationTarget(1);
    else
        m.setNormalisationTarget(-1);
    addMeasurement(m);
}
void combiner::printCorrelationMatrix()const{
    if(debug)
        std::cout << "combiner::printCorrelationMatrix" <<std::endl;
}
void combiner::readConfigFile(const std::string & filename){


    configfile_=filename;
    fileReader fr;
    fr.setComment("#");
    fr.setDelimiter(",");
    fr.setStartMarker("[global]");
    fr.setEndMarker("[end global]");
    fr.readFile(configfile_);
    fr.setRequireValues(false);
    isdifferential_ = fr.getValue<bool>("isDifferential",false);
    float normconstr=0;
    if(fr.getValue<bool>("normalise",false))
        normconstr=100;
    norm_constraint_ = fr.getValue<float>("normConstraint",normconstr);
    normalised_input_ = fr.getValue<bool>("normalisedInput",false);


    fr.setRequireValues(true);
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

    if(debug)
        std::cout << "correlations" << std::endl;

    fr.setStartMarker("[correlations]");
    fr.setEndMarker("[end correlations]");
    fr.setDelimiter("=");
    fr.readFile(configfile_);
    for(size_t i=0;i<fr.nLines();i++){
        if(fr.nEntries(i)<2)continue;
        correlationscan cscans;
        TString sysname=fr.getData<TString>(i,0);
        size_t idxa=0;
        try{
            idxa=external_correlations_.getEntryIndex(sysname);
        }catch(...){
            throw std::runtime_error(("combiner::setup: reading in correlation assumptions: failed to find systematic "+sysname).Data());
        }

        std::string corrinfo=fr.getData<std::string>(i,1);
        textFormatter tf;
        tf.setDelimiter("+");
        std::vector<std::string> contributions=tf.getFormatted(corrinfo);
        std::string sysnames=",";
        for(size_t j=0;j<contributions.size();j++){
            single_correlationscan cscan;
            cscan.idxa=idxa;
            cscan.nominal=0;
            cscan.high=0;
            cscan.low=0;
            cscan.namea=sysname;

            size_t idxb=0;

            tf.setDelimiter(")");
            tf.setTrim(" (");
            std::vector<std::string> entry=tf.getFormatted(contributions.at(j));
            if(entry.size() < 1)continue;
            try{
                idxb=external_correlations_.getEntryIndex(entry.at(1).data());
            }catch(...){
                throw std::runtime_error(("combiner::setup: reading in correlation assumptions: failed to find systematic "+(TString)entry.at(1)).Data());
            }
            cscan.idxb=idxb;
            cscan.nameb=entry.at(1).data();
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


            }
            cscans.add(cscan);
            sysnames+=entry.at(1)+"+";
        }
        sysnames=std::string(sysnames.begin(),sysnames.end()-1);
        //get the name right right away
        cscans.setName(sysname.Data()+sysnames);
        syst_scanranges_.push_back(cscans);
    }
    //fill nominal assumptions
    for(size_t i=0;i<syst_scanranges_.size();i++){
        for(size_t j=0;j<syst_scanranges_.at(i).size();j++){
            setSystCorrelation(syst_scanranges_.at(i).get(j).idxa,syst_scanranges_.at(i).get(j).idxb,
                    syst_scanranges_.at(i).get(j).nominal);

        }
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

void combiner::checkConsistency()const{
    if(isdifferential_ && excludebin_>=0 && norm_constraint_ == 0)
        throw std::logic_error("Trying to exclude a bin without normalisation constraint on result will never converge.");
    if(!isdifferential_ && (excludebin_>=0 || norm_constraint_>0))
        throw std::logic_error("Input marked as not differential but trying to exclude a bin or applying norm constraint. Please check consistency of input configuration");
    if(!isdifferential_ && normalised_input_)
        throw std::logic_error("Input marked as not differential but as normalised. Please check consistency of input configuration");

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



std::vector<std::vector<combinationResult> > combiner::scanCorrelations(std::ostream& out, const combinationResult& nominal, const std::string& outdir)const{

    /*
     * First part:
     *   Do the scan and save output to temp
     */
    const size_t nscans=syst_scanranges_.size();
    if(nscans<1)
        return std::vector<std::vector<combinationResult> > ();

    std::vector<std::vector<combinationResult> > results;
    for(const auto& ss:syst_scanranges_){
        results.push_back(std::vector<combinationResult>());
        ss.size();//to avoid compiler warning
    }

    int ndone=0;

    //parallelise the main loop here - retains ordering, but not optimal parallelisation

#ifdef USE_MP
#pragma omp parallel for
#endif
    for(size_t i=0;i<syst_scanranges_.size();i++){

        std::vector<combinationResult> thisscan(single_correlationscan::nPoints());

        for(size_t step=0;step<single_correlationscan::nPoints();step++){
            combinationResult result;
            combiner combcp=*this;
            for(size_t j=0;j<syst_scanranges_.at(i).size();j++){
                combcp.setSystCorrelation(syst_scanranges_.at(i).get(j).idxa,
                        syst_scanranges_.at(i).get(j).idxb,
                        syst_scanranges_.at(i).get(j).scanVal(step));
            }
            combinationResult res;
            try{
                res=combcp.combine();

#ifdef USE_MP
#pragma omp critical (cout)
#endif
                {
                    std::cout << "performed "<< ndone << " of "<< syst_scanranges_.size()*single_correlationscan::nPoints() << " scans."<<std::endl;
                }
            }catch(...){
                std::cout << "one scan point for "<< syst_scanranges_.at(i).name() <<" failed: ";
                for(size_t j=0;j<syst_scanranges_.at(i).size();j++){
                    std::cout  << syst_scanranges_.at(i).get(j).scanVal(step) << " ";
                }
                std::cout  <<"ignoring" <<std::endl;
            }
            thisscan.at(step)=res;
            ndone++;
        }
#ifdef USE_MP
#pragma omp critical (scanpointresult)
#endif
        {
            results.at(i)=thisscan;
        }
    }



    //old

    const size_t nobs=tobecombined_.size();



    /*
     * Second part:
     *   Write output as graphs to TFile
     */
    TFile * tfile=new TFile((TString)outprefix_+"scanPlots.root","RECREATE");

    for(size_t i=0;i<results.size();i++){
        //const correlationscan * s= & syst_scanranges_.at(i);
        const TString & name= syst_scanranges_.at(i).name() ;
        TString filename=textFormatter::makeCompatibleFileName(name.Data());
        if(filename.Length()>100)
            filename=TString(filename,100);
        filename+="_";
        filename+=i;

        for(size_t obs=0;obs<nobs;obs++){
            /*
             * Prepare graph content input
             */
            TString graphname=textFormatter::makeCompatibleFileName(( tobecombined_.at(obs).first+"_"+filename).Data());
            std::vector<double> scan,zero,comb,errup,errdown;
            std::vector<double> minchi2;
            for(size_t j=0;j<results.at(i).size();j++){
                const combinationResult& result=results.at(i).at(j);
                if(result.getCombNames().size()<1) continue; //ignore failed scans

                if(result.getNCombined() != nobs){//sanity check for debugging
                    TString errstr=nominal.getCombNames().at(obs)+"_"+name ;
                    throw std::logic_error("combiner::scanCorrelations: fatal error. not all (or too many) combined results: "+(std::string)errstr.Data());
                }
                comb.push_back(result.combined_.at(obs));
                errup.push_back(result.comberrup_.at(obs));
                errdown.push_back(-result.comberrdown_.at(obs));
                minchi2.push_back(result.getChi2min());

                if(syst_scanranges_.at(i).isSingle())
                    scan.push_back(syst_scanranges_.at(i).get(0).scanVal(j));
                else
                    scan.push_back(1./(single_correlationscan::nPoints()-1)*(float)j);//no better way...
                zero.push_back(0);
            }
            /*
             * Create graph and add some cosmetics
             */
            TGraphAsymmErrors * g=new TGraphAsymmErrors(comb.size(), &scan.at(0), &comb.at(0),
                    &zero.at(0), &zero.at(0), &errup.at(0), &errdown.at(0));
            TGraphAsymmErrors * gline=new TGraphAsymmErrors(comb.size(), &scan.at(0), &comb.at(0),
                    &zero.at(0), &zero.at(0), &zero.at(0), &zero.at(0));

            g->GetYaxis()->SetTitle(tobecombined_.at(obs).first);
            if(syst_scanranges_.at(i).isSingle()){
                applyGraphCosmetics(g,gc_scancombined,syst_scanranges_.at(i).getLowest(),
                        syst_scanranges_.at(i).getHighest(),graphname,name);
                applyGraphCosmetics(gline,gc_scancombined,syst_scanranges_.at(i).getLowest(),
                        syst_scanranges_.at(i).getHighest(),graphname,name);
            }
            else{
                applyGraphCosmetics(g,gc_multiscan,syst_scanranges_.at(i).getLowest(),
                        syst_scanranges_.at(i).getHighest(),graphname,name);
                applyGraphCosmetics(gline,gc_multiscan,syst_scanranges_.at(i).getLowest(),
                        syst_scanranges_.at(i).getHighest(),graphname,name);

            }

            g->SetName(graphname);
            g->SetTitle(graphname);
            g->Write();

            TGraphAsymmErrors * gminchi2=new TGraphAsymmErrors(comb.size(), &scan.at(0), &minchi2.at(0),
                    &zero.at(0), &zero.at(0), &zero.at(0), &zero.at(0));
            applyGraphCosmetics(gminchi2,gc_minchi2,syst_scanranges_.at(i).getLowest(),
                    syst_scanranges_.at(i).getHighest(), "min#chi^2",name);
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

                double nomcorr=0;
                if(syst_scanranges_.at(i).isSingle())
                    nomcorr=nominal.getInputSysCorrelations().getEntry(syst_scanranges_.at(i).get(0).idxa,syst_scanranges_.at(i).get(0).idxb);
                double nomerrd=-nominal.getCombErrdown().at(obs);
                TGraphAsymmErrors gnom(1, &nomcorr, &nominal.getCombined().at(obs),
                        &zero.at(0), &zero.at(0), &nominal.getCombErrup().at(obs), &nomerrd);

                applyGraphCosmetics(&gnom,gc_nominal,syst_scanranges_.at(i).getLowest(),
                        syst_scanranges_.at(i).getHighest(),graphname+"_nom",name);
                cv.Draw();
                g->Draw("Aa3pl");
                gline->Draw("l");//again
                if(syst_scanranges_.at(i).isSingle())
                    gnom.Draw("Pe");


                cv.Print((TString)outdir+"/"+nominal.getCombNames().at(obs)+"_"+filename +".pdf");

                cv.Clear();
                applyGraphCosmetics(g,gc_scancombinedUP,syst_scanranges_.at(i).getLowest(),
                        syst_scanranges_.at(i).getHighest(),graphname,name,2.05);
                applyGraphCosmetics(gline,gc_scancombinedUP,syst_scanranges_.at(i).getLowest(),
                        syst_scanranges_.at(i).getHighest(),graphname,name,2.05);

                cv.Divide(1,2);
                TVirtualPad * pad=cv.cd(1);
                pad->SetBottomMargin(0.015);
                pad->SetLeftMargin(.15);
                pad->SetTopMargin(.189);
                g->Draw("Aa3pl");
                gline->Draw("l");//again
                if(syst_scanranges_.at(i).isSingle())
                    gnom.Draw("Pe");

                if(syst_scanranges_.at(i).isSingle())
                    applyGraphCosmetics(gminchi2,gc_minchi2,syst_scanranges_.at(i).getLowest(),
                            syst_scanranges_.at(i).getHighest(), "min#chi^2",name,2.05);
                else
                    applyGraphCosmetics(gminchi2,gc_minchi2multiscan,syst_scanranges_.at(i).getLowest(),
                            syst_scanranges_.at(i).getHighest(), "min#chi^2",name,2.05);

                pad=cv.cd(2);
                pad->SetBottomMargin(0.293);
                pad->SetTopMargin(0.025);
                pad->SetLeftMargin(.15);
                gminchi2->Draw("Al");



                cv.Print((TString)outdir+"/"+nominal.getCombNames().at(obs)+"_"+filename+"_detailed.pdf");

            }


            // g is owned by TFile, no delete necessary
        }
    }

    tfile->Close();
    delete tfile;

    for(const auto& sv:results){
        for(const auto& s:sv){
            s.printFullInfo(out);
        }
    }


    return results;
}


std::vector<combinationResult>
combiner::scanExcludeBins(std::ostream& out, const combinationResult& nominal)const{
    std::vector<combinationResult> output;
    if(!isdifferential_) return output;
    if(measurements_.size()<1){
        throw std::out_of_range("combiner::scanLeastExcludeBins: no measurements associated");
    }
    int nombin = measurements_.at(0).getLeastSignificantBin();
    int nbins = measurements_.at(0).getEstimates().size();

    for(int i=0;i<nbins;i++){
        combiner ccp=*this;
        if(i==nombin)continue;
        ccp.setExcludeBin(i);
        combinationResult res = ccp.combinePriv();
        output.push_back(res);
    }

    //make the output nicely readable: (only relative differnece per bin w.r.t. total unc and nominal value)
    const auto nomcomb = nominal.getCombined();
    const auto nomerr = nominal.getCombSymmErr();
    const auto names = nominal.getCombNames();

    out << "d/sigma: difference to nominal exclude bin result w.r.t. nominal total uncertainty\n";
    out << "dsigma/sigma: difference to nominal exclude bin error w.r.t. nominal error\n";
    out << "d/nom: difference to nominal exclude bin result w.r.t nominal value\n";

    for(size_t i=0;i<output.size();i++){
        out << "\n\nExclude bin " << output.at(i).getExcludeBin() << ":\n\n";
        out << "   bin     " << "  |  d/sigma [%]  |   dsigma/sigma [%] |   d/nom [%]\n";
        auto scanerr = output.at(i).getCombSymmErr();

        for(size_t j=0;j<nomcomb.size();j++){
            out << std::setw(11) <<  names.at(j) << "  | ";

            double diff = (output.at(i).getCombined().at(j) - nomcomb.at(j));
            double d_o_sigma = 100. * diff / nomerr.at(j);
            double dsigma_o_sigma = 100. * scanerr.at(j) / nomerr.at(j);
            double d_o_nom = 100. * diff / nomcomb.at(j);

            out << std::setprecision(3) << std::setw(11) << d_o_sigma << "   | ";
            out << std::setprecision(3) << std::setw(17) << dsigma_o_sigma << "  |  ";
            out << std::setprecision(3) << std::setw(5) << d_o_nom << "\n";
        }
    }

    out << "\n\n\n --- detailed information ---\n\n\n";
    out << "Nominal exclude bin "  << nominal.getExcludeBin() << "\n\n";
    nominal.printFullInfo(out);
    for(const auto & c: output){
        out << "\n\n\n --------- \n\n\n";
        out << "Exclude bin " << c.getExcludeBin() << ":\n\n";
        c.printFullInfo(out);
    }



    return output;
}


combinationResult combiner::combine()const{
    if(measurements_.size()<1)
        throw std::logic_error("combiner::combine: no measurement to combine associated");
    combiner cp=*this;
    return cp.combinePriv();
}

combinationResult combiner::combinePriv(){
    if(debug)
        std::cout << "combiner::combine" <<std::endl;

    checkConsistency();

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
                if(thisval){
                    thisval /= m.getNormalisationScaling();
                    mean+=thisval;

                }
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


    if(excludebin_ >= 0){
        for(auto & m: measurements_)
            m.setExcludeBin(excludebin_);
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

    out.post_all_correlations_ = fitter.getCorrelationMatrix();

    if(debug && nsyst<20)
        std::cout << "post combine systematics correlations: \n"<< out.post_sys_correlations_ << std::endl;



    if(debug){
        std::cout << "correlations between combined quantities" <<std::endl;
        std::cout << out.post_meas_correlations_ << std::endl;
    }


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
    out.excludebin_ = excludebin_;

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

void combiner::setExcludeBin(int bin){
    if((size_t)bin>=0 && measurements_.size() && (size_t)bin >= measurements_.at(0).getEstimates().size())
        throw std::out_of_range("combiner::setExcludeBin: bin out of range of input bins");

    excludebin_=bin;
}

void combiner::setExcludeBinAuto(){
    /*
     * check for normalised differential extras
     */
    if(isdifferential_ && normalised_input_ ){
        excludebin_ = measurements_.at(0).getLeastSignificantBin();
    }
}



void combiner::clear(){
    external_correlations_.clear();
    tobecombined_.clear();
    measurements_.clear();
}

