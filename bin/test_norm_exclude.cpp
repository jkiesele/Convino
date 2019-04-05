#include "combiner.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "fstream"
#include "measurement.h"
#include "math.h"
#include "normaliser.h"

/*
 *  This is an example code for the combination of e.g. differential cross sections using
 *  the C++ interface directly to root histograms
 *
 *  To create own executables, this file can be used as template.
 *  All .cpp files located in the directory 'bin' will be compiled when running 'make'
 *  This feature can be used to easily create own executables
 */


// The following functions are only helpers to create data for pseudo measurements.
// in reality, this would be the input from the measurements

const int nbins = 5;
bool useflat=false;

TH1D createPseudoMeasurement( TString name, int statistics);
TH1D createPseudoMeasurementFlat( TString name, int statistics);
TH1D createSystematicVariation(const TH1D& nominal, double scale, TString name);
TH2D createCovarianceMatrix(const TH1D& nominal, double correlation_strength);

int main(){

    /*
     * For reference create a TFile to save input and output
     */

    TFile f("differentialExample.root","RECREATE");


    /*
     * Create some pseudo measurements and a correlation matrix
     */

    TH1::SetDefaultSumw2(true);



    TH1D histo_meas1=createPseudoMeasurement("histo_meas1",50);
    histo_meas1.Scale(20);
   // TH1D syshisto_meas1=createSystematicVariation(histo_meas1,1.2,"syshisto_meas1");
    TH2D cov_1 = createCovarianceMatrix(histo_meas1,0.2);

    // create a second measurement
    TH1D histo_meas2=createPseudoMeasurement("histo_meas2",1100);
   // TH1D syshisto_meas2=createSystematicVariation(histo_meas2,1.5,"syshisto_meas2");
    TH2D cov_2 = createCovarianceMatrix(histo_meas2,0.2);


    combiner comb;

    measurement m1;
    m1.setMeasured(&histo_meas1);
   // m1.setEstimateCovariance(cov_1);
    comb.addMeasurement(m1);


    measurement m2;
    m2.setMeasured(&histo_meas2);
  //  m1.setEstimateCovariance(cov_2);
    comb.addMeasurement(m2);


    combinationResult comb_result=comb.combine();

    comb_result.printFullInfo(std::cout);
    //return 1;

    std::cout << ">>>>>>>>>>>>>>> normalise.." << std::endl;

    normaliser norm;
    norm.setInput(comb_result);
    auto normres = norm.getNormalised();
    normres.printFullInfo(std::cout);

    std::cout << ">>>>>>>>>>>>>>>>>>normalise first" << std::endl;


    norm.clear();
    norm.setInput(&histo_meas1,&cov_1);
    TH1D histo_meas1_norm = *norm.getNormalisedTH1D();
    TH2D cov_1_norm = *norm.getNormalisedCovarianceTH2D();

    norm.getNormalised().printFullInfo(std::cout);

    norm.setInput(&histo_meas2,&cov_2);
    TH1D histo_meas2_norm = *norm.getNormalisedTH1D();
    TH2D cov_2_norm = *norm.getNormalisedCovarianceTH2D();


    norm.getNormalised().printFullInfo(std::cout);


    std::cout << ">>>>>>>>>>>>>>>>>>adding to combiner" << std::endl;

    combiner comb_norm;
   // measurement::debug=true;
    measurement m1_norm;
    m1_norm.setExcludeBin(0);
    m1_norm.setMeasured(&histo_meas1_norm);
    m1_norm.setEstimateCovariance(cov_1_norm);
    comb_norm.addMeasurement(m1_norm);


    measurement m2_norm;
    m2_norm.setExcludeBin(0);
    m2_norm.setMeasured(&histo_meas2_norm);
    m2_norm.setEstimateCovariance(cov_2_norm);
    comb_norm.addMeasurement(m2_norm);

    combinationResult comb_result_norm=comb_norm.combine();

    comb_result_norm.printFullInfo(std::cout);

    norm.clear();
    norm.setInput(comb_result_norm);
    norm.getNormalised().printFullInfo(std::cout);


}




/*
 *
 *   THIS PART IS ONLY TO CREATE THE PSEUDO MEASUREMENTS ETC.
 *   THERE IS NO NEED TO INCLUDE ANY OF THIS IN AN ACTUAL COMBINATION OF REAL MEASUREMENTS
 *
 */


TH1D createPseudoMeasurement(TString name,  int statistics){

    if(useflat)
        return createPseudoMeasurementFlat(name,statistics);
    //if(withUFOF)
    //   nbins+=2;
    TH1D h(name,name,nbins,-2,2);

    h.FillRandom("gaus",statistics);

    return h;
}

TH1D createPseudoMeasurementFlat( TString name, int statistics){


    double nperbin = (double)statistics/(double)nbins;
    TH1D h(name,name,nbins,-2,2);
    for(int i=1;i<=nbins;i++){
        h.SetBinContent(i,nperbin);
        h.SetBinError(i,std::sqrt(nperbin));
    }
    return h;
}


TH1D createSystematicVariation(const TH1D& nominal, double scale, TString name){
    /*
     * create another pseudo systematic variation
     */
    TH1D syshisto_meas2=nominal;
    syshisto_meas2.SetName(name);
    syshisto_meas2.Scale(scale);
    syshisto_meas2.Write();
    return syshisto_meas2;
}

TH2D createCovarianceMatrix(const TH1D& nominal, double correlation_strength){
    TH2D  h(nominal.GetName() + (TString)"_covariance",nominal.GetName() + (TString)"_covariance",
            nominal.GetNbinsX(),0,nominal.GetNbinsX(),
            nominal.GetNbinsX(),0,nominal.GetNbinsX());

    for(int i=1; i<=nominal.GetNbinsX();i++){
        double stati=nominal.GetBinError(i);
        for(int j=1; j<=nominal.GetNbinsX();j++){
            double statj=nominal.GetBinError(j);
            double correlation = exp(-1/correlation_strength*(i-j)*(i-j));
            h.SetBinContent(i,j, stati*statj*correlation );
        }
    }
    return h;
}

