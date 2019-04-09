#include "combiner.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "fstream"
#include "measurement.h"
#include "math.h"
#include "normaliser.h"
#include "TCanvas.h"
#include "TStyle.h"


/*
 *  This is an example code for the combination of e.g. differential cross sections using
 *  the C++ interface directly to root histograms
 *
 *  To create own executables, this file can be used as template.
 *  All .cpp files located in the directory 'bin' will be compiled when running 'make'
 *  This feature can be used to easily create own executables
 */

TH1D * getWithChecks(TFile&, const TString& name);
TH2D * getWithChecks2D(TFile&, const TString& name);
void divideByBinwidth(TH1D& h);
void divideCovByBinwidth(const TH1D& h, TH2D& cov);
measurement readInDataset(const TString& basedir, const TString& filename, const TString& prefix, const std::vector<TString>& sysnamesbare, TH1D& refhisto);

int main(){

    /*
     * For reference create a TFile to save input and output
     */
    std::vector<TString> distributions={
            "toppt","ttm","topy","topbarpt","tleppt","thadpt+ttm+cts","thardpt","tty","topbary","thadpt+ttm+tty","tsoftpt","cts","ttpt","thady","thadpt","tlepy"
    };

    gStyle->SetOptStat(0);

    TString basedir="/Users/jkiesele/temp/Convino/otto/";
    std::string correlationfile="/Users/jkiesele/temp/Convino/otto/convino_config/sys_correlations.txt";

    for(const auto& distribution : distributions){
        std::vector<TString> sysnames={
                "JER","JES","BDecay","FS","BTAG","PU","ISR","RS","BFrag","BKGSTNORM","BTAGL","LEP","FSR","BKG","HD","BKGSNORM","LUMI","TUNE","PDF","PDFAS","MTOP"
        };
        std::vector<TString> years={"18","17","16"};
        std::vector<measurement> measurements;

        TH1D refhisto; //for the binning

        for(const auto& year : years){

            TString readfile = "xsec_"+year+"/"+"xsec_"+distribution+".root";
            std::cout << readfile <<std::endl;
            measurement m = readInDataset(basedir,readfile,year+"_",sysnames,refhisto);
            measurements.push_back(m);
            if(year =="18")
                sysnames.push_back("PREFIRE");//use for 17 and 16
        }

        combiner comb;
        for(const auto& m: measurements)
            comb.addMeasurement(m);

        comb.readCorrelationFile(correlationfile);

        auto result = comb.combine();

        TFile fout("combined_"+distribution+".root","RECREATE");
        result.fillTH1(&refhisto);
        refhisto.SetName("combined");
        refhisto.Write();

        TCanvas cv;
        refhisto.Draw();
        cv.Print("combined_"+distribution+".pdf");

        std::ofstream output_textfile("combined_"+distribution+".txt");
        result.printFullInfo(output_textfile);
        output_textfile.close();

        fout.Close();

    }
}


void divideByBinwidth(TH1D& h){
    return ;
    for(int i=1;i<=h.GetNbinsX();i++){
        h.SetBinContent(i, h.GetBinContent(i)/h.GetBinWidth(i));
        h.SetBinError(i, h.GetBinError(i)/h.GetBinWidth(i));
    }

}
void divideCovByBinwidth(const TH1D& h, TH2D& cov){
    return ;
    for(int i=1;i<cov.GetNbinsX();i++){
        for(int j=0;j<cov.GetNbinsY();j++){
            cov.SetBinContent(i,j, cov.GetBinContent(i,j)/(h.GetBinWidth(i)*h.GetBinWidth(j)));
        }
    }
}

TH1D * getWithChecks(TFile& f, const TString&  name){
    TH1D * h = (TH1D*) f.Get(name);
    if(!h || h->IsZombie())
        throw std::runtime_error((name+" not found in "+f.GetPath()).Data());
    return h;
}

TH2D * getWithChecks2D(TFile& f, const TString&  name){
    TH2D * h = (TH2D*) f.Get(name);
    if(!h || h->IsZombie())
        throw std::runtime_error((name+" not found in "+f.GetPath()).Data());
    return h;
}

measurement readInDataset(const TString& basedir, const TString& filename, const TString& prefix, const std::vector<TString>& sysnamesbare, TH1D& refhisto){

    TFile fin(basedir+"/"+filename);
    TString add="_merged";
    measurement m;

    TH1D * nominal = getWithChecks(fin, "xsec"+add);
    divideByBinwidth(*nominal);
    m.setMeasured(nominal);
    refhisto = *nominal;

    TH2D * cov = getWithChecks2D(fin,"cov_stat"+add);
    divideCovByBinwidth(*nominal,*cov);
    m.setEstimateCovariance(*cov);

    for(const auto& sys: sysnamesbare){
       // std::cout << "adding systematics " << sys << " to " << prefix <<" --> "<< prefix+sys << std::endl;
        TH1D * sysh = getWithChecks(fin, (TString)(sys+"_up"+add));
        sysh->Add(nominal);
        m.addSystematics(prefix+sys,sysh);
    }

    return m;

}

