/*
 * main.cc
 *
 *  Created on: 30 Jun 2016
 *      Author: jkiesele
 */


#include "combiner.h"
#include <iostream>
#include "TMatrixD.h"
#include <fstream>
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH2D.h"


void coutHelp(){
	std::cout << "\nUSAGE: convino <options> <inputfile>\n";
	std::cout << "OPTIONS:\n";
	std::cout << "\n-h --help       display this help message\n";
	/*std::cout << "\n-e              make an inaccurate fit first before\n";
	std::cout << "                the final fit to get better starting\n";
	std::cout << "                parameters\n";
	std::cout << "                (slower, can help in case of convergence problems)\n"; */
    std::cout << "\n-c              creates contour plots for all";
    std::cout << "                 combined values w.r.t. all uncertainties\n";
	std::cout << "\n-d              switch on debug mode (more printout)\n";
    std::cout << "\n-e              enables scan of different exclude bins for\n";
    std::cout << "                normalised differential inputs\n";
	std::cout << "\n-s              enable scan of correlation assumptions\n";
	std::cout << "\n-p              create plots for the correlation scans in .pdf format\n";
    std::cout << "                (a directory will be created)\n";
    std::cout << "\n-r              create root output file <results.root> with result (TGraph and TH2D for matrices)\n";
	std::cout << "\n--prefix <pref> define a prefix for the output files/directories\n";
    std::cout << "\n--excludebin <bin> force specific exclude bin\n";
    std::cout << "\n--neyman        use neyman chi2 (faster)\n";
    std::cout << "\n--fast          switches on fast mode (no covariance output)\n";
	std::cout << std::endl;
	std::cout << "EXAMPLE: convino -sp examples/exampleconfig.txt" <<std::endl;
	std::cout << std::endl;
}


int main(int argc, char* argv[]){

	combiner comb;
	comb.setMode(combiner::lh_mod_pearson);
	std::string outputprefix;
	bool writeScanPlotPdfs=false;
	if(argc<2)
		return -1;
	std::string infile;
	TString roothistoinfile="";
	TString roothistoin="";
	bool doscan=false;
	bool doexcludebinscan=false;
	bool rootoutput=false;
	bool contours=false;
	int excludebin=-1;
	//very simple parsing
	for(int i=1;i<argc;i++){
		TString targv=argv[i];
		if(targv.BeginsWith("--")){
		    if(targv == ("--help")){
		                    coutHelp();
		                    exit(0);
		                }
		    else if(targv == ("--neyman"))
				comb.setMode(combiner::lh_mod_neyman);
            else if(targv == ("--fast"))
                comb.setFastMode(true);
		    else if(targv == ("--prefix")){
				if(i+1>=argc || ((TString)argv[i+1]).BeginsWith("-")){
					std::cerr << "please specify a valid prefix, e.g. --prefix <prefix>" <<std::endl;
					exit(-1);
				}
				outputprefix=argv[++i];
				comb.setOutputPrefix(outputprefix);
			}
			else if(targv == ("--excludebin")){
                if(i+1>=argc || ((TString)argv[i+1]).BeginsWith("-")){
                    std::cerr << "please specify a valid input, e.g. --excludebin 3" <<std::endl;
                    exit(-1);
                }
                excludebin = atoi(argv[++i]);
            }
            else if(targv == ("--ri")){ //TBI upon request
                if(i+2>=argc || ((TString)argv[i+1]).BeginsWith("-") || ((TString)argv[i+2]).BeginsWith("-")){
                    std::cerr << "please specify a valid input, e.g. --ri myfile.root myhisto" <<std::endl;
                    exit(-1);
                }
                roothistoinfile = argv[++i];
                roothistoin = argv[++i];
            }

		}
		else if(targv.BeginsWith("-")){ //simple option
			/*if(targv.Contains("e"))
				comb.setEstimateFirst(true);*/
			if(targv.Contains("d"))
				combiner::debug=true;
			if(targv.Contains("s"))
				doscan=true;
			if(targv.Contains("p"))
				writeScanPlotPdfs=true;
			if(targv.Contains("e"))
			    doexcludebinscan=true;
            if(targv.Contains("r"))
                rootoutput=true;
            if(targv.Contains("c"))
                contours=true;
			if(targv.Contains("h")){
				coutHelp();
				exit(0);
			}
			if(targv.Contains("t") && i+1<argc){
			//	comb.setResolveThreshold(atof(argv[i+1]));
			}
		}
		else
			infile=argv[i];
	}
	if(infile.length()<1){
		coutHelp();
		return -1;
	}
	if(outputprefix.length())
		outputprefix+="_";
	if(!rootoutput && roothistoinfile.Length()>1){
	    std::cout << "Warning: The input root file and histograms are only used for binning the output graph in x. If no root output is requested, they are ignored."<<std::endl;
	}

	if(excludebin>=0){
	    comb.setExcludeBin(excludebin);
	}

	comb.readConfigFile(infile);

	comb.setAllContours(contours);

	//comb.setup();

	std::cout << "setup done, processing..."<<std::endl;



	combinationResult result=comb.combine();

	result.printResultOnly(std::cout);

	std::ofstream logfile((outputprefix+"result.txt").data());
	result.printFullInfo(logfile);
	logfile.close();

	std::cout << "done. saved output to "<< outputprefix+"result.txt"<<std::endl;


	if(doscan){
		std::ofstream logfilescan((outputprefix+"scan_result.txt").data());
		if(writeScanPlotPdfs)
			comb.scanCorrelations(logfilescan,result,(outputprefix+"scan_results").data());
		else
			comb.scanCorrelations(logfilescan,result);
		logfilescan.close();
	}
	if(doexcludebinscan){
	    if(!comb.isNormalisedDifferentialInput()){
	        std::cout << "exclude bin scan does not make sense if input is not normalised differential. skipping" <<std::endl;
	    }
	    else{
	        std::ofstream logfilescan((outputprefix+"excludebin_scan_result.txt").data());
	        comb.scanExcludeBins(logfilescan,result);
	        logfilescan.close();
	    }
	}
	if(rootoutput){
	    TGraphAsymmErrors * g=0;

	    //TBI upon request
	    if(roothistoinfile){
	        TFile fin(roothistoinfile,"READ");
	        if(fin.IsZombie()){
	            std::cout << "Could not read input root file " + roothistoinfile << std::endl;
	            exit(-2);
	        }

	    }


	    TFile f("result.root","RECREATE");

	    result.fillTGraphAsymmErrors(g);
	    g->SetMarkerStyle(21);
	    g->SetMarkerColor(kBlack);
        g->SetName("combined_result");
	    g->Write();

	    //result.getCor

	    f.Close();

	}
	if(contours){
	    result.writeAllContourPlots(outputprefix+"contours.root");
	    result.saveAllContourPlots(outputprefix+"contours");
	}

	return 0;
}
