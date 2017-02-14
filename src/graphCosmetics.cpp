/*
 * graphCosmetics.cpp
 *
 *  Created on: 28 Oct 2016
 *      Author: jkiesele
 */


#include "graphCosmetics.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"

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
	else if(c==gc_scancombined || c==gc_scancombinedUP){
		g->GetYaxis()->SetTitleOffset(1/scaler);
		g->GetYaxis()->CenterTitle(true);
	//	if(gc_scancombinedUP)
	//		g->GetXaxis()->SetLabelSize(scaler*0.0);
		g->SetMarkerColor(kBlack);
		g->SetMarkerSize(0);
		g->SetMarkerStyle(20);
		g->SetFillColor(kBlack);
		g->SetFillStyle(3002);
		g->SetDrawOption("a3pl");

	}
	else if(c==gc_minchi2){
		g->GetYaxis()->SetTitleOffset(1/scaler);
		g->GetYaxis()->CenterTitle(true);


	}


	if(xlow<xhigh)
		g->GetXaxis()->SetRangeUser(xlow-0.1,xhigh+0.1);
	else
		g->GetXaxis()->SetRangeUser(xhigh-0.1,xlow+0.1);

	g->GetXaxis()->SetTitle(scanned);


}

