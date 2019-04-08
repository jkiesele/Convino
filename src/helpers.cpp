/*
 * helpers.cpp
 *
 *  Created on: 5 Apr 2019
 *      Author: jkiesele
 */

#include "helpers.h"
#include "TString.h"



TH2D removeOneBin(const TH2D& in, int bin){
    if(bin<0)return in;
    if(in.GetNbinsX() != in.GetNbinsY())
        throw std::out_of_range("removeOneBin: TH2D: only for symmetric inputs");
    TH2D h = TH2D(in.GetName()+(TString)"_rem",
            in.GetTitle()+(TString)"_rem",
             in.GetNbinsX()-1,0,in.GetNbinsX(),
             in.GetNbinsY()-1,0,in.GetNbinsY());

    int counteri=1;
    int counterj=1;
    for(int i=1;i<=in.GetNbinsX();i++){
        counterj=1;
        if(i == bin+1)continue;
        for(int j=1;j<=in.GetNbinsX();j++){
            if(j == bin+1)continue;
            h.SetBinContent(counteri,counterj,in.GetBinContent(i,j));
            counterj++;
        }
        counteri++;
    }
    return h;

}

TH1D removeOneBin(const TH1D& in, int bin){
    if(bin<0)return in;
    TH1D h = TH1D(in.GetName()+(TString)"_rem",
            in.GetTitle()+(TString)"_rem",
            in.GetNbinsX()-1,0,in.GetNbinsX());

    int counteri=1;
    for(int i=1;i<=in.GetNbinsX();i++){
        if(i == bin+1)continue;
        h.SetBinContent(counteri,in.GetBinContent(i));
        counteri++;
    }
    return h;
}
