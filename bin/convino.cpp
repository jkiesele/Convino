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


void coutHelp(){
	std::cout << "\nUSAGE: convino <options> <inputfile>\n";
	std::cout << "OPTIONS:\n";
	std::cout << "\n-h --help       display this help message\n";
	/*std::cout << "\n-e              make an inaccurate fit first before\n";
	std::cout << "                the final fit to get better starting\n";
	std::cout << "                parameters\n";
	std::cout << "                (slower, can help in case of convergence problems)\n"; */
	std::cout << "\n-d              switch on debug mode (more printout)\n";
	std::cout << "\n-s              enable scan of correlation assumptions\n";
	std::cout << "\n-p              create plots for the correlation scans in .pdf format\n";
	std::cout << "                (a directory will be created)\n";
	std::cout << "\n--neyman        use neyman chi2 (faster)\n";
	std::cout << "\n--prefix <pref> define a prefix for the output files/directories\n";
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
	bool doscan=false;
	//very simple parsing
	for(int i=1;i<argc;i++){
		TString targv=argv[i];
		if(targv.BeginsWith("--")){
			if(targv == ("--neyman"))
				comb.setMode(combiner::lh_mod_neyman);
			if(targv == ("--prefix")){
				if(i+1>=argc || ((TString)argv[i+1]).BeginsWith("-")){
					std::cerr << "please specify a valid prefix, e.g. --prefix <prefix>" <<std::endl;
					exit(-1);
				}
				outputprefix=argv[++i];
				comb.setOutputPrefix(outputprefix);
			}
			if(targv == ("--help")){
				coutHelp();
				exit(0);
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

	comb.readConfigFile(infile);

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

	return 0;
}
