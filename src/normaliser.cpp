/*
 * normaliser.cpp
 *
 *  Created on: 5 Apr 2019
 *      Author: jkiesele
 */




#include "normaliser.h"
#include "multiGausGenerator.h"
#include "TH1D.h"
#include "TH2D.h"


void normaliser::setInput(const combinationResult& res){
    //do some consistency checks here maybe TBI

    in_res_=res;
}


void normaliser::setInput(const TH1D* h1d, const TH2D* h2d){
    in_res_ = combinationResult();
    in_res_.post_meas_covariance_ = triangularMatrix(h1d->GetNbinsX());
    in_res_.post_meas_covariance_.fillFromTH2(*h2d);


    in_res_.combined_.resize(h1d->GetNbinsX());
    in_res_.combnames_.clear();
    for(int i=0;i<h1d->GetNbinsX();i++){
        in_res_.combined_.at(i)=h1d->GetBinContent(i+1);
        TString n="comb_";
        n+=i;
        in_res_.combnames_.push_back(n);
    }
    setInput(in_res_);
}

combinationResult normaliser::getNormalised(double addUncertFraction, int bin)const{
    //set up inputs for normalised

    auto cov_in = in_res_.getCombinedCovariance();
    auto nominal = in_res_.getCombined();
    auto varied = nominal;
    double sum_nominal=0;
    for(const auto& v:nominal)
        sum_nominal+=v;
    const size_t nbins = cov_in.size();

    multiGausGenerator gen(seed_);
    gen.setCovariance(cov_in);
    double addbincontent=0;

    std::vector<std::vector<double> > covfilled(cov_in.size(),std::vector<double>(nbins,0));
    combinationResult out=in_res_;

    if(addbin_>=0){
        //generate expectation value for this bin

        for(size_t i=0;i<iterations_;i++){
            double sum_varied=0;
            auto vec = gen.generate();
            for(size_t i=0;i<nbins;i++){
                varied[i] = nominal[i] + vec[i];
                sum_varied+=varied[i];
            }
            addbincontent += 1 - sum_varied;
        }
        addbincontent /= (double)iterations_;
        covfilled.resize(cov_in.size()+1,std::vector<double>(nbins+1,0));
        sum_nominal=1;
    }



    for(size_t i=0;i<iterations_;i++){
        auto vec = gen.generate();

        //exploit loop here for sum
        double sum_varied=0;
        for(size_t i=0;i<nbins;i++){
            varied[i] = nominal[i] + vec[i];
            sum_varied+=varied[i];
        }
        double addbinthis = 1 - sum_varied;
        if(addbin_>=0)
            sum_varied=1;

        //fill cov and use norm, avoid additional loop / exploit symmetry
        size_t jcounter=0;
        size_t icounter=0;
        for(int i=0;i<(int)covfilled.size();i++){
            jcounter=0;
            double vari = varied[icounter]/sum_varied -  nominal[icounter]/sum_nominal;
            if(i==addbin_)
                vari = addbinthis - addbincontent;

            for(int j=0;j<=i;j++){
                double varj = varied[jcounter]/sum_varied -  nominal[jcounter]/sum_nominal;
                if(j==addbin_)
                    varj = addbinthis - addbincontent;
                covfilled[i][j] += vari*varj;

                if(j!=addbin_)
                    jcounter++;
            }
            if(i!=addbin_)
                icounter++;
        }

    }
    double sum_all = 0;
    for(const auto& v:nominal)
        sum_all+=v;
    sum_all+=addbincontent;



    triangularMatrix newcov = cov_in;
    if(addbin_>=0){
        auto names = cov_in.createNamesVector();
        std::vector<TString> newnames(names.size()+1);
        int counter=0;
        for(int i=0;i<(int)newnames.size();i++){
            if(addbin_!=i){
                newnames.at(i) = names.at(counter);
                counter++;
            }
            else{
                newnames.at(i) = addbin_name;
            }
        }
        newcov = triangularMatrix(newnames);
    }


    for(size_t i=0;i<covfilled.size();i++){
        for(size_t j=0;j<=i;j++){
            if((int)i==bin && (int)j==bin)
                newcov.setEntry(i,j, covfilled[i][j]*(1.+addUncertFraction) / (double)iterations_);
            else
                newcov.setEntry(i,j, covfilled[i][j] / (double)iterations_);
        }
    }


    std::vector<double> normed(newcov.size());
    int counter=0;
    for(int i=0;i<(int)normed.size();i++){
        if(addbin_!=i){
            normed.at(i) = nominal.at(counter)/sum_all;
            counter++;
        }
        else{
            normed.at(i) = addbincontent/sum_all;
        }
    }

    auto correlationm = newcov;
    correlationm.normalize();
    //exploit cov symmetry
    out.post_meas_covariance_=newcov;
    out.post_meas_correlations_ = correlationm;
    out.combined_ = normed;
    if(out.comberrup_.size()<newcov.size())
        out.comberrup_.resize(newcov.size());
    for(size_t i=0;i<newcov.size();i++)
        out.comberrup_.at(i) = std::sqrt(newcov.getEntry(i,i));
    out.comberrdown_=out.comberrup_;
    out.combnames_ = newcov.createNamesVector();
    out.post_all_correlations_.clear();//
    out.post_sys_correlations_.clear();
    out.constraints_.clear();
    out.pulls_.clear();
    out.chi2min_=-1;


    return out;//fornow
}


TH1D * normaliser::getNormalisedTH1D()const{
    if(out_res.combined_.size()<1)
        out_res=getNormalised();

    TH1D * h = new TH1D("normalised","normalised", out_res.combined_.size(), 0 , out_res.combined_.size());
    for(int i=0;i<(int)out_res.combined_.size();i++){
        h->SetBinContent(i+1, out_res.combined_.at(i));
        h->SetBinError(i+1 , out_res.comberrup_.at(i));
    }
    return h;

}
TH2D * normaliser::getNormalisedCovarianceTH2D(double addUncertFraction, int bin)const{
    if(out_res.combined_.size()<1)
        out_res=getNormalised(addUncertFraction,bin);

    TH2D * h = new TH2D("normalised_covariance","normalised_covariance", out_res.combined_.size(), 0 , out_res.combined_.size(),
             out_res.combined_.size(), 0 , out_res.combined_.size());
    for(int i=0;i<(int)out_res.combined_.size();i++){
        for(int j=0;j<(int)out_res.combined_.size();j++){
            h->SetBinContent(i+1,j+1,  out_res.post_meas_covariance_.getEntry(i,j));
            h->SetBinError(i+1 , j+1,  out_res.post_meas_covariance_.getEntry(i,j));
        }
    }
    return h;
}

double normaliser::normaliseVector(std::vector<double> & v)const{
    double sum=0;
    for(const auto& i:v)
        sum+=i;
    for(auto& i:v)
        i/=sum;
    return sum;
}
