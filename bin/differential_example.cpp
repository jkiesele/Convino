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

TH1D createPseudoMeasurement( TString name, bool withUFOF, int statistics, bool normalise=false);
TH1D createPseudoMeasurementFlat( TString name, bool withUFOF, int statistics, bool normalise=false);
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

    const bool withUFOF=false; //this is not fully validated, yet!

    const bool normalise_input=false;

    const bool flat=false;

    TH1D histo_meas1=createPseudoMeasurement("histo_meas1",withUFOF,1000,normalise_input);
    if(flat)
         histo_meas1=createPseudoMeasurementFlat("histo_meas1",withUFOF,200,normalise_input);
    histo_meas1.Write(); //just for bookkeeping

    /*
     * create a pseudo systematic variation
     */
    TH1D syshisto_meas1=createSystematicVariation(histo_meas1,1.2,"syshisto_meas1");
    syshisto_meas1.Write(); //just for bookkeeping


    /*
     * Create a covariance
     */

    TH2D cov_1 = createCovarianceMatrix(histo_meas1,0.02);
    cov_1.Write(); //just for bookkeeping



    // create a second measurement
    TH1D histo_meas2=createPseudoMeasurement("histo_meas2",withUFOF,1200,normalise_input);
    if(flat)
             histo_meas2=createPseudoMeasurementFlat("histo_meas2",withUFOF,400,normalise_input);
    histo_meas2.Write(); //just for bookkeeping

    /*
     * create a pseudo systematic variation for the second measurement
     */
    TH1D syshisto_meas2=createSystematicVariation(histo_meas2,1.5,"syshisto_meas2");
    syshisto_meas2.Write(); //just for bookkeeping

    /*
     * Create a covariance for measurement 2
     * with a bit stronger correlations
     */

    TH2D cov_2 = createCovarianceMatrix(histo_meas2,0.1);
    cov_2.Write(); //just for bookkeeping



    /*
     * Fill the information in the combiner.
     * This is the 'real' part of the program
     */

    //measurement::debug=true;

    std::cout << "here "<<std::endl;

    combiner comb;
    comb.setMode(combiner::lh_mod_neyman);

    measurement m1;
    m1.setMeasured(&histo_meas1);


    //first covariance, then systematics, otherwise association is ambiguous (and will lead to error messages)
    //m1.setEstimateCovariance(cov_1);
    m1.addSystematics("syst_a",&syshisto_meas1);

    comb.addMeasurement(m1);


    measurement m2;
    m2.setMeasured(&histo_meas2);
    //m2.setEstimateCovariance(cov_2);
    /*
     * The uncertainties MUST NOT have the same name, even if they are fully correlated. This case is defined later.
     */
    m2.addSystematics("syst_1",&syshisto_meas2);
    //m2.setParameterType("syst_1",parameter::para_unc_relative);
    comb.addMeasurement(m2);

    /*
     * Measurements are added. Define correlations between the systematics
     */
    //comb.setSystCorrelation("syst_a","syst_1",0.5);


    //normalise
   // comb.setNormConstraint(1e2);//norm pull always small...
   // comb.setExcludeBin(-1);


    //combiner::debug=true;


    combinationResult comb_result=comb.combine();





    std::cout << comb_result.getCombinedCovariance() << std::endl;

    std::cout << "done combining"<<std::endl;
    /*
     * ****** safe and print the output
     */
    comb_result.printFullInfo(std::cout);


    normaliser norm;
    norm.setInput(comb_result);
    auto normres = norm.getNormalised();
    normres.printFullInfo(std::cout);








    return 1;

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
    TGraphAsymmErrors * tg = new  TGraphAsymmErrors(&result);
    tg->SetName("result_tgraph_asymmerr");
    comb_result.fillTGraphAsymmErrors(tg, true);
    tg->Write();
    f.Close();


    return 0;
}




/*
 *
 *   THIS PART IS ONLY TO CREATE THE PSEUDO MEASUREMENTS ETC.
 *   THERE IS NO NEED TO INCLUDE ANY OF THIS IN AN ACTUAL COMBINATION OF REAL MEASUREMENTS
 *
 */


TH1D createPseudoMeasurement(TString name, bool withUFOF, int statistics, bool normalise){

    if(withUFOF) //in this function done
        throw std::runtime_error("inclusion of underflow and overflow is not implemented fully consistently everywhere yet. Please make sure your input does not include those.");

    //if(withUFOF)
     //   nbins+=2;
    TH1D h(name,name,nbins,-2,2);

    h.FillRandom("gaus",statistics);

    double integ=h.Integral();
    if(normalise)
        h.Scale(1/integ);


        return h;
}

TH1D createPseudoMeasurementFlat( TString name, bool withUFOF, int statistics, bool normalise){if(withUFOF) //in this function done
    if(withUFOF) //in this function done
        throw std::runtime_error("inclusion of underflow and overflow is not implemented fully consistently everywhere yet. Please make sure your input does not include those.");


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
            double correlation = exp(-correlation_strength*(i-j)*(i-j));
            h.SetBinContent(i,j, stati*statj*correlation );
        }
    }
    return h;
}

