/*
 * differentialExample.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: jkiesele
 */


#include "combiner.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "fstream"
#include "measurement.h"

/*
 *  This is an example code for the combination of e.g. differential cross sections using
 *  the C++ interface directly to root histograms
 *
 *  To create own executables, this file can be used as template.
 *  All .cpp files located in the directory 'bin' will be compiled when running 'make'
 *  This feature can be used to easily create own executables
 */
TH1D createPseudoMeasurement(int nbins, TString name, bool withUFOF, int statistics){

	if(withUFOF)
		nbins+=2;
	TH1D h(name,name,nbins,-1,1);

	h.FillRandom("gaus",statistics);

	h.Scale(1000/statistics); //make them similar

	if(withUFOF){
		h.SetName("tmp");
		TH1D h2(name,name,nbins-2,-1,1);
		for(int i=0;i<=nbins-1;i++){
			h2.SetBinContent(i,h.GetBinContent(i+1));
			h2.SetBinError(i,h.GetBinError(i+1));
		}
		return h2;
	}
	else
		return h;
}


int main(){

	/*
	 * For reference create a TFile to save input and output
	 */

	TFile f("differentialExample.root","RECREATE");


	/*
	 * Create some pseudo measurements and a correlation matrix
	 */

	TH1::SetDefaultSumw2(true);

	const bool withUFOF=true;

	const size_t nbins=10;

	TH1D histo_meas1=createPseudoMeasurement(nbins,"histo_meas1",withUFOF,1000);
	histo_meas1.Write();

	/*
	 * create a pseudo systematic variation
	 * e.g. a variation of a global scale factor for histo_meas2
	 */
	TH1D syshisto_meas1=histo_meas1;
	syshisto_meas1.SetName("syshisto_meas1");
	syshisto_meas1.Scale(1.03);
	syshisto_meas1.Write();


	TH1D histo_meas2=createPseudoMeasurement(nbins,"histo_meas2",withUFOF,500);
	histo_meas2.Write();

	/*
	 * create another pseudo systematic variation
	 */
	TH1D syshisto_meas2=histo_meas2;
	syshisto_meas2.SetName("syshisto_meas2");
	syshisto_meas2.Scale(0.99);
	syshisto_meas2.Write();

	/*
	 * A correlation matrix between the bins
	 */
	int corrdimension=nbins;
	if(withUFOF)
		corrdimension+=2;
	double corrs[corrdimension][corrdimension];
	for(int i=0;i<corrdimension;i++){
		for(int j=0;j<corrdimension;j++){
			if(i==j)
				corrs[i][j]=1;
			else if(i+1==j || j+1==i)
				corrs[i][j]=0.3;
			else if(i+2==j || j+2==i)
				corrs[i][j]=0.1;
			else
				corrs[i][j]=0;
		}
	}



	/*
	 * Fill the information in the combiner.
	 * This is the 'real' part of the program
	 */

	combiner comb;

	measurement m1;
	m1.setMeasured(&histo_meas1);
	m1.addSystematics("syst_a",&syshisto_meas1);

	for(int i=0;i<corrdimension;i++){
		for(int j=0;j<corrdimension;j++){
			m1.setEstimateCorrelation(i,j,corrs[i][j]);
		}
	}

	comb.addMeasurement(m1);

	measurement m2;
	m2.setMeasured(&histo_meas2);
	/*
	 * The uncertainties MUST NOT have the same name, even if they are fully correlated. This case is defined later.
	 */
	m2.addSystematics("syst_1",&syshisto_meas2);

	//for simplicity the same correlationsa are used.
	for(int i=0;i<corrdimension;i++){
		for(int j=0;j<corrdimension;j++){
			m2.setEstimateCorrelation(i,j,corrs[i][j]);
		}
	}

	comb.addMeasurement(m2);

	/*
	 * Measurements are added. Define correlations between the systematics
	 */
	comb.setSystCorrelation("syst_a","syst_1",0.2);

	combinationResult comb_result=comb.combine();


	/*
	 * ****** safe and print the output
	 */
	comb_result.printFullInfo(std::cout);

	//can also be printed to an ofstream:
	std::ofstream output_textfile("differentialExample_output.txt");
	comb_result.printFullInfo(output_textfile);
	output_textfile.close();

	/*
	 *  ** safe the output to root a histogram and graph
	 */
	TH1D result("result","result",nbins,-1,1);

	comb_result.fillTH1(&result);
	result.Write();


	/*
	 * For asymmetric errors, the output can be written to a TGraphAsymmErrors
	 *
	 * The x positions of the points need to be the right ones
	 * therefore, the result TH1 is used as reference.
	 */
	TGraphAsymmErrors tg(&result);
	tg.SetName("result_tgraph_asymmerr");
	comb_result.fillTGraphAsymmErrors(&tg, true);
	tg.Write();
	f.Close();


	return 0;
}

