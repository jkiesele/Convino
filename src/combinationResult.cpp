/*
 * combinationResult.cpp
 *
 *  Created on: 17 Aug 2016
 *      Author: jkiesele
 */
#include "combinationResult.h"
#include <iostream>
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "textFormatter.h"
#include "helpers.h"


combinationResult::combinationResult():chi2min_(0),isdifferential_(false),excludebin_(-1){}

void combinationResult::printResultOnly(std::ostream& out)const{
    out << "combined (minimum chi^2="<<chi2min_<<"):"<<std::endl;
    for(size_t i=0;i<combnames_.size();i++)
        out << combnames_.at(i) << ": " << combined_.at(i)
        << " +"<<comberrup_.at(i) << " " << comberrdown_.at(i) << std::endl;

    //calculate simple weight of each measurement
    post_meas_covariance_.clear();
}

void combinationResult::printSimpleInfo(std::ostream& out)const{
    out << "pre-combine systematics correlations:" <<std::endl;
    out << orig_sys_correlations_ << std::endl;

    out << "post-combine systematics correlations:" <<std::endl;
    out << post_sys_correlations_ << std::endl;

    out << "pre-combine estimate correlations:" <<std::endl;
    out << orig_meas_correlations_ << std::endl;

    out << "post-combine result correlations:" <<std::endl;
    out << post_meas_correlations_ << std::endl;

    out << "post-combine result covariance:" <<std::endl;
    out << getCombinedCovariance() << std::endl;

    printResultOnly(out);
}

void combinationResult::printFullInfo(std::ostream& out)const{
    printSimpleInfo(out);

    out << "full correlation matrix:" <<std::endl;
    out << post_all_correlations_ << std::endl;

    out << std::endl;
    out << textFormatter::fixLength("Name",11);
    out << textFormatter::fixLength("pull",7);
    out << "constraint\n";
    size_t maxlength=0;
    for(size_t i=0;i<post_sys_correlations_.size();i++){
        size_t tmp=post_sys_correlations_.getEntryName(i).Length();
        if(tmp>maxlength)maxlength=tmp;
    }
    for(size_t i=0;i<post_sys_correlations_.size();i++){
        out << textFormatter::fixLength(post_sys_correlations_.getEntryName(i).Data(),maxlength+1)<<" ";
        out << textFormatter::fixLength(toString(pulls_.at(i)),4) << "   ";
        out << textFormatter::fixLength(toString(constraints_.at(i)),4) << std::endl;
    }
    out << std::endl;
    printImpactTable(out);

}

void combinationResult::printImpactTable(std::ostream& out)const{
    out  << "Impact table: name, impact [%]" <<std::endl;
    out << textFormatter::fixLength(" ",15) <<" ";
    for(size_t c=0;c<combined_.size();c++){
        out << textFormatter::fixLength(combnames_.at(c).Data(), 15) <<" | ";
    }
    out << std::endl;
    for(const auto& i:impacttable_){
        out << textFormatter::fixLength(i.first.Data(),15) << " ";
        for(size_t c=0;c<combined_.size();c++){
            double rel = fabs(i.second.at(c) / combined_.at(c))*100.;
            out << textFormatter::fixLength(toString(rel),15) << " | ";
        }
        out << std::endl;
    }


}

void combinationResult::fillTH1(TH1*h)const{
    int bins=h->GetNbinsX();
    int start=1;
    int end=bins+1;
    /*
     * for future implementations: consder underflow/overflow
     */
    if((int)combined_.size() != end-start)
        throw std::out_of_range("combinationResult::fillTH1: bins don't match");


    auto err = getCombSymmErr();
    for(size_t i=0;i<combined_.size();i++){
        h->SetBinContent(i+start, combined_.at(i));
        h->SetBinError(i+start,err.at(i));
    }
}

std::vector<double> combinationResult::getCombSymmErr() const {
    std::vector<double> out;
    for(size_t i=0;i<combined_.size();i++){
        double maxerr=std::max(fabs(comberrup_.at(i)) ,fabs(comberrdown_.at(i)));
        out.push_back(maxerr);
    }
    return out;
}

void combinationResult::fillTGraphAsymmErrors(TGraphAsymmErrors*& h)const{
    int bins=0;
    bool isempty=false;
    TString name="combined";
    if(h){
        bins=h->GetN();
        if(bins<1)
            isempty=true;
        if(((TString)h->GetName()).Length())
            name = h->GetName();
    }
    else
        isempty=true;



    int start=0;
    int end=bins;

    if(bins==0){//empty graph
        start=0;
        end=combined_.size();
        bins=combined_.size();
        if(h)
            delete h;
        h= new TGraphAsymmErrors(bins);
        h->SetName(name);
    }

    if((int)combined_.size() != bins+start+(bins-end)){
        throw std::out_of_range("combinationResult::fillTGraphAsymmErrors: number of points doesn't match");
    }



    for(int i=0;i<end;i++){
        double pointx,pointy;
        h->GetPoint(i,pointx,pointy);
        if(isempty)
            pointx=i+1;
        h->SetPoint(i,pointx, combined_.at(i+start));
        h->SetPointEYlow(i, fabs(comberrdown_.at(i+start) ) );
        h->SetPointEYhigh(i, comberrup_.at(i+start));
    }

}

const triangularMatrix& combinationResult::getCombinedCovariance() const {
    if(post_meas_covariance_.size())
        return post_meas_covariance_;

    post_meas_covariance_.clear();
    post_meas_covariance_ = triangularMatrix(post_meas_correlations_.createNamesVector()); //could work

    auto symerr = getCombSymmErr();

    for(size_t i=0;i<post_meas_covariance_.size();i++){
        for(size_t j=0;j<=i;j++){
            post_meas_covariance_.setEntry(i,j,post_meas_correlations_.getEntry(i,j) * symerr.at(i)*symerr.at(j));
        }
    }

    return post_meas_covariance_;
}

void combinationResult::reduce(){

    orig_sys_correlations_.clear();
    orig_meas_correlations_.clear();
    post_sys_correlations_.clear();
    post_meas_correlations_.clear();
    post_meas_covariance_.clear();
    post_all_correlations_.clear();
    pulls_.clear();
    constraints_.clear();

}
void combinationResult::copyFrom(const combinationResult& r){
    combnames_=r.combnames_;
    combined_=r.combined_;
    comberrup_=r.comberrup_;
    comberrdown_=r.comberrdown_;
    orig_sys_correlations_=r.orig_sys_correlations_;
    orig_meas_correlations_=r.orig_meas_correlations_;
    post_sys_correlations_=r.post_sys_correlations_;
    post_meas_correlations_=r.post_meas_correlations_;
    post_all_correlations_=r.post_all_correlations_;
    post_meas_covariance_=r.post_meas_covariance_;
    pulls_=r.pulls_;
    constraints_=r.constraints_;
    chi2min_=r.chi2min_;
    isdifferential_=r.isdifferential_;
    excludebin_=r.excludebin_;
}
