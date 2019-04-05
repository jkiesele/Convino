/*
 * root_to_measured.cpp
 *
 *  Created on: 27 Mar 2019
 *      Author: jkiesele
 */

#include <iostream>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include <algorithm>

void coutHelp(){
    std::cout << "program to convert a root th1x to measurement description in text file format for convino. " ;
    std::cout << "works only for orthogonal input systematics. Otherwise the full hessian needs to be parsed. The latter is not supported by this small tool and needs to be done by the user.";
    std::cout << "\nUSAGE:\root_to_measured <options> <inputfile nominal> <histogram name nominal> ";
    std::cout <<  "<syst a name> <inputfile syst a> <histogram name syst a> ";
    std::cout <<  "<syst b name> <inputfile syst b> <histogram name syst b> ...\n";
    std::cout << "best pipe the output to a configuration file\n";
    std::cout << "OPTIONS:\n";
    std::cout << "\n-h              display this help message\n";

    std::cout << "\n-p              measurement prefix name (will be added to \"binX\")";
    std::cout << "\n-a              add to nominal histogram to get syst. varied histogram";
    std::cout << "\n--noufof        do not consider under and overflow bins";

    std::cout << std::endl;
    std::cout << "EXAMPLE: root_to_matrix  -o cov my_file.root my_histo_2d" <<std::endl;
    std::cout << std::endl;
}

std::vector<double> extractHisto(const TString& file, const TString& histo, bool ignoreufof, std::vector<double>& stat){

    TFile f(file);
    auto* o=f.Get(histo);
    if(!o){
        std::cout << "object "<< histo << " not found" <<std::endl;
        exit(-1);
    }
    TString classname = o->ClassName();
    if(!classname.Contains("TH1")){
        std::cout << histo << " is not a TH1x histogram" <<std::endl;
        exit(-1);
    }
    TH1 * h = (TH1*)o;

    int xbins = h->GetNbinsX();

    TString output="";
    int maxbin=xbins+2;
    int minbin=0;

    if(!xbins){
        std::cout << "histogram does not contain any bin"<<std::endl;
        exit(-1);
    }

    if(!ignoreufof &&  (h->GetBinContent(0,0)==0 && h->GetBinError(0,0)==0)){
        std::cout << "histogram "<<file <<" :/ "<< histo<<"does not contain any underflow bin, please run with --noufof\n";
        std::cout << "also, please notice that migrations out of the phase space cannot be taken into account this way."<<std::endl;
        exit(-1);
    }
    if(ignoreufof){
        maxbin--;
        minbin++;
    }
    std::vector<double> out;
    stat.clear();

    for(int i=minbin;i<maxbin;i++){
        out.push_back(h->GetBinContent(i));
        stat.push_back(h->GetBinError(i));
    }


    return out;

}

int main(int argc, char* argv[]){



    TString prefix="";
    std::vector<TString> opts;

    bool ignoreufof=false;
    bool addtonominal=false;

    for(int i=1;i<argc;i++){
        TString targv=argv[i];
        if(targv == "--noufof"){
            ignoreufof=true;
        }
        else if(targv.BeginsWith("-")){

            if(targv.Contains("h")){
                coutHelp();
                exit(0);
            }
            else if(targv == "-p"){
                if(i+1>=argc || ((TString)argv[i+1]).BeginsWith("-")){
                    std::cerr << "please specify a valid option argument" <<std::endl;
                    coutHelp();
                    exit(-1);
                }
                TString next=argv[++i];
                prefix=next;
            }

            else if(targv.Contains("a")){
                addtonominal=true;
            }
            else{

                std::cerr << "please specify a valid option argument" <<std::endl;
                coutHelp();
                exit(-1);
            }

        }
        else{
            opts.push_back(targv);
        }
    }
    if(opts.size()<2 || (opts.size() - 2 )% 3){
        coutHelp();
        exit(-1);
    }

    //get nominal
    std::vector<double> stat ,dummy;
    std::vector<double>  nominal = extractHisto(opts.at(0),opts.at(1),ignoreufof,stat);
    std::vector<std::pair<TString,std::vector<double> > > systs;
    for(size_t i=2;i<opts.size();i++){
        TString sysname = opts.at(i);
        std::cout << sysname << std::endl;
        std::vector<double> sys = extractHisto(opts.at(i+1),opts.at(i+2),ignoreufof,dummy);
        i+=2;
        systs.push_back(std::pair<TString,std::vector<double> > (sysname,sys));
    }

    TString out="[estimates]\n  n_estimates = ";

    out+=nominal.size();
    out+="\n";
    for(size_t i=0;i<nominal.size();i++){
        out+="\n  name_";
        out+=i;
        out+=" = "+prefix+"bin";
        out+=i;
        out+="\n  value_";
        out+=i;
        out+=" = ";
        out+=nominal.at(i);
        out+="\n";
    }
    out+="[end estimates]";
    std::cout << out << std::endl;


    out="[not fitted]\n";
    for(size_t i=0;i<(size_t)prefix.Length()+8;i++)
        out+=" ";

    for(size_t i=0;i<systs.size();i++){
        out+=systs.at(i).first+" ";
    }
    out+=" stat";
    out+="\n";
    for(size_t i=0;i<nominal.size();i++){
        out+=prefix+"bin";
        out+=i;
        out+="  ";
        double nominal_val=nominal.at(i);
        if(addtonominal)
            nominal_val=0;
        for(size_t j=0;j<systs.size();j++){
            out+=systs.at(j).second.at(i)-nominal_val;
            out+=" ";
        }
        out+=stat.at(i);
        out+="\n";
    }


    out+="[end not fitted]\n";
    std::cout << out << std::endl;

}





