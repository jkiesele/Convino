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

combinationResult::combinationResult():chi2min_(0),isdifferential_(false),hasUF_(false),hasOF_(false),excludebin_(-1){}

void combinationResult::printResultOnly(std::ostream& out)const{
	out << "combined (minimum chi^2="<<chi2min_<<"):"<<std::endl;
	for(size_t i=0;i<combnames_.size();i++)
		out << combnames_.at(i) << ": " << combined_.at(i)
		<< " +"<<comberrup_.at(i) << " " << comberrdown_.at(i) << std::endl;

	//calculate simple weight of each measurement

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


}

void combinationResult::fillTH1(TH1*h)const{
	int bins=h->GetNbinsX();
	int start=1;
	int end=bins+1;
	if(hasUF_)
		start--;
	if(hasOF_)
		end++;
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

void combinationResult::fillTGraphAsymmErrors(TGraphAsymmErrors* h, bool cutUFOF)const{
	int bins=h->GetN();

	int start=0;
	int end=bins;
	if(hasUF_ && cutUFOF)
		start++;
	if(hasOF_ && cutUFOF)
		end--;

	if((int)combined_.size() != bins+start+(bins-end)){
		if(!cutUFOF)
			throw std::out_of_range("combinationResult::fillTGraphAsymmErrors: number of points doesn't match. input could have underflow/overflow and corresponding bins are not foreseen in the output graph. Try the cutUFOF option.");
		throw std::out_of_range("combinationResult::fillTGraphAsymmErrors: number of points doesn't match");
	}



	for(int i=0;i<=end;i++){
		double pointx,pointy;
		h->GetPoint(i,pointx,pointy);
		h->SetPoint(i,pointx, combined_.at(i+start));
		h->SetPointEYlow(i, fabs(comberrdown_.at(i+start) ) );
		h->SetPointEYhigh(i, comberrup_.at(i+start));
	}

}


void combinationResult::reduce(){

	orig_sys_correlations_.clear();
	orig_meas_correlations_.clear();
	post_sys_correlations_.clear();
	post_meas_correlations_.clear();
	post_all_correlations_.clear();
	pulls_.clear();
	constraints_.clear();

}
