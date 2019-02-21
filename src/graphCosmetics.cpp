/*
 * graphCosmetics.cpp
 *
 *  Created on: 28 Oct 2016
 *      Author: jkiesele
 */


#include "graphCosmetics.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "textFormatter.h"
#include "TStyle.h"

void applyGraphCosmetics(TGraphAsymmErrors* g, graphcosmetics c, double xlow, double xhigh, const TString & graphname, const TString & scanned,double scaler){



	g->GetXaxis()->SetLabelSize(scaler*0.05);
	g->GetXaxis()->SetTitleSize(scaler*0.06);
	g->GetXaxis()->SetTitleOffset(1.);
	g->GetYaxis()->SetLabelSize(scaler*0.05);
	g->GetYaxis()->SetTitleSize(scaler*0.06);
	g->GetYaxis()->SetTitleOffset(0.85/scaler);
	g->GetYaxis()->SetNdivisions(505);

	if(c==gc_nominal){

		g->SetMarkerSize(g->GetMarkerSize()+1);
		g->SetMarkerColor(kBlack);
		g->SetMarkerStyle(27);
		g->SetLineColor(kBlack);
		g->SetLineWidth(2);
	}
	else if(c==gc_scancombined || c==gc_scancombinedUP || c==gc_multiscan){
		g->GetYaxis()->SetTitleOffset(1/scaler);
		g->GetYaxis()->CenterTitle(true);
	//	if(gc_scancombinedUP)
	//		g->GetXaxis()->SetLabelSize(scaler*0.0);
		g->SetMarkerColor(kBlack);
		g->SetMarkerSize(0);
		g->SetMarkerStyle(20);
		g->SetFillColor(kGray);
		gStyle->SetHatchesLineWidth(1);
		gStyle->SetHatchesSpacing(0.5);
		g->SetFillStyle(3305);
		g->SetDrawOption("l");

	}
	else if(c==gc_minchi2|| c==gc_minchi2multiscan){
		g->GetYaxis()->SetTitleOffset(1/scaler);
		g->GetYaxis()->CenterTitle(true);


	}
	if(c==gc_multiscan || c==gc_minchi2multiscan){
	    //maybe additions
	}


	if(xlow<xhigh)
		g->GetXaxis()->SetRangeUser(xlow-0.1,xhigh+0.1);
	else
		g->GetXaxis()->SetRangeUser(xhigh-0.1,xlow+0.1);

	std::string xaxis;
	textFormatter tf;
	tf.setDelimiter(",");
	std::string first=tf.getFormatted(scanned.Data()).at(0);
	std::string second=tf.getFormatted(scanned.Data()).at(1);

	tf.setDelimiter("+");
	std::vector<std::string> fmt=tf.getFormatted(second);

	if(fmt.size()>1)
	    xaxis+="scan: ";
	for(const auto& s:fmt){
	    xaxis+="#rho(";
	    xaxis+=first+",";
	    xaxis+=s;
	    xaxis+=")&";
	}
	xaxis=std::string(xaxis.begin(),xaxis.end()-1);
	if(xaxis.length()>40)
	    g->GetXaxis()->SetTitleSize(g->GetXaxis()->GetTitleSize()*(40./xaxis.length()));

	g->GetXaxis()->SetTitle(xaxis.data());


}

