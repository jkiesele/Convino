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

const int nbins = 3;
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
    TGraphAsymmErrors * g = new TGraphAsymmErrors();

    /*
     * Create some pseudo measurements and a correlation matrix
     */

    TH1::SetDefaultSumw2(true);



    TH1D histo_meas1=createPseudoMeasurement("histo_meas1",10000);
    histo_meas1.Scale(1./100.);
    // TH1D syshisto_meas1=createSystematicVariation(histo_meas1,1.2,"syshisto_meas1");
    TH2D cov_1 = createCovarianceMatrix(histo_meas1,0.1);

    // create a second measurement
    TH1D histo_meas2=createPseudoMeasurement("histo_meas2",11000);
    histo_meas2.Scale(1./110.);
    // TH1D syshisto_meas2=createSystematicVariation(histo_meas2,1.5,"syshisto_meas2");
    TH2D cov_2 = createCovarianceMatrix(histo_meas2,0.1);





   // return 1;

    combiner comb;

   // comb.setExcludeBin(0);

    measurement m1;
   // m1.setExcludeBin(1);
    m1.setMeasured(&histo_meas1);

   // m1.setIsNormalisedInput(true);
    m1.setEstimateCovariance(cov_1);
    m1.addSystematics("scalesys1",1.1);
    comb.addMeasurement(m1);


    measurement m2;
  //  m2.setExcludeBin(1);
    m2.setMeasured(&histo_meas2);
  //  m2.setIsNormalisedInput(true);
    m2.setEstimateCovariance(cov_2);
    m2.addSystematics("scalesys2",1.1);
    comb.addMeasurement(m2);


    combinationResult comb_result=comb.combine();

    g->SetName("combined_non_norm");
    comb_result.fillTGraphAsymmErrors(g);
    g->Write();

    comb_result.printFullInfo(std::cout);
    //return 1;

    std::cout << ">>>>>>>>>>>>>>> normalise.." << std::endl;

    normaliser norm;
    //norm.setIterations(1e8);
    norm.setInput(comb_result);
    auto normres = norm.getNormalised();
    normres.printFullInfo(std::cout);
   // norm.setIterations(1e6);

    auto corr_combnorm = normres.getCombinedCorrelations();

    g = new TGraphAsymmErrors();
    g->SetName("first_com_then_norm");
    normres.fillTGraphAsymmErrors(g);
    g->Write();

    std::cout << ">>>>>>>>>>>>>>>>>>normalise first" << std::endl;


    norm.clear();
    norm.setInput(&histo_meas1,&cov_1);
    TH1D histo_meas1_norm = *norm.getNormalisedTH1D();
    TH2D cov_1_norm = *norm.getNormalisedCovarianceTH2D(0.2,2);

    std::cout << ">>>>>meas1 normed" << std::endl;
   // norm.getNormalised().printFullInfo(std::cout);

    norm.clear();
    norm.setInput(&histo_meas2,&cov_2);
    TH1D histo_meas2_norm = *norm.getNormalisedTH1D();
    TH2D cov_2_norm = *norm.getNormalisedCovarianceTH2D(0.2,2);


    std::cout << ">>>>>meas2 normed" << std::endl;
  //  norm.getNormalised().printFullInfo(std::cout);


    std::cout << ">>>>>>>>>>>>>>>>>>remove a bin" << std::endl;

    int removebin=1;

    histo_meas1_norm = removeOneBin(histo_meas1_norm,removebin);
    cov_1_norm = removeOneBin(cov_1_norm,removebin);
    histo_meas2_norm = removeOneBin(histo_meas2_norm,removebin);
    cov_2_norm = removeOneBin(cov_2_norm,removebin);




    std::cout << ">>>>>>>>>>>>>>>>>>adding to combiner" << std::endl;

    combiner comb_norm;

    measurement m1_norm;
    m1_norm.setMeasured(&histo_meas1_norm);
    m1_norm.setEstimateCovariance(cov_1_norm);
    m1_norm.addSystematics("scalesys",1.1);
    comb_norm.addMeasurement(m1_norm);


    measurement m2_norm;
    m2_norm.setMeasured(&histo_meas2_norm);
    m2_norm.setEstimateCovariance(cov_2_norm);
    comb_norm.addMeasurement(m2_norm);

    combinationResult comb_result_norm=comb_norm.combine();

    comb_result_norm.printFullInfo(std::cout);


    std::cout << ">>>>>>>>>>>>>>>>>>normalise again, addbin bin back" << std::endl;

    size_t iterations=1e6;

    norm.clear();
    norm.setIterations(iterations);
    norm.addFloatingBin(removebin,"added");
    norm.setInput(comb_result_norm);
    auto last_res = norm.getNormalised();
    last_res.printFullInfo(std::cout);

    auto corr_normcomb = last_res.getCombinedCorrelations();

   std::cout << "iterations " << iterations << " Freb dist " <<  corr_normcomb.FrobeniusDistance(corr_combnorm) << std::endl;

   corr_normcomb.removeEntry(removebin);
   corr_combnorm.removeEntry(removebin);
   std::cout << "only non removed part Freb dist " <<  corr_normcomb.FrobeniusDistance(corr_combnorm) << std::endl;



   g = new TGraphAsymmErrors();
    last_res.fillTGraphAsymmErrors(g);
    g->Write();
    f.Close();


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
    TH1D h(name,name,nbins,-5,5);

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
