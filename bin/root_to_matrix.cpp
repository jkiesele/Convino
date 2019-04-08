/*
 * root_to_matrix.cpp
 *
 *  Created on: 27 Mar 2019
 *      Author: jkiesele
 */

#include <iostream>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include "triangularMatrix.h"
#include <sstream>

void coutHelp(){
    std::cout << "program to convert a root th2x to a correlation matrix, hessian or covariance matrix to convino configuration text file format." ;
    std::cout << "\nUSAGE:\root_to_matrix <options> <inputfile> <histogram name>\n";
    std::cout << "best pipe the output to a configuration file\n";
    std::cout << "OPTIONS:\n";
    std::cout << "\n-h              display this help message\n";

    std::cout << "\n-i              input format.\n";
    std::cout << "                  Options: cov(ariance) (default, will be converted to hessian),\n";
    std::cout << "                  hes(sian), cor(relation) matrix\n";
    std::cout << "\n-p              measurement prefix name (will be added to \"binX\")";

    std::cout << std::endl;
    std::cout << "EXAMPLE: root_to_matrix  -o cov my_file.root my_histo_2d" <<std::endl;
    std::cout << std::endl;
}

TString extractMatrix(const TString& file, const TString& histo, TString prefix , bool convertcovtohes){

    TFile f(file);
    auto* o=f.Get(histo);
    if(!o){
        std::cout << "object "<< histo << " not found" <<std::endl;
        exit(-1);
    }
    TString classname = o->ClassName();
    if(!classname.Contains("TH2")){
        std::cout << histo << " is not a TH2x histogram" <<std::endl;
        exit(-1);
    }
    TH2 * h = (TH2*)o;

    int xbins = h->GetNbinsX();
    int ybins = h->GetNbinsY();
    if(xbins!=ybins || !xbins){
        std::cout << "input matrix is not symmetric - please check format"<<std::endl;
        exit(-1);
    }
    TString output="";
    int maxbin=xbins+1;
    int minbin=1;


    if(!convertcovtohes){

        int bincounter=0;
        for(int i=minbin;i<maxbin;i++){
            output+=prefix+"bin";
            output+=bincounter;
            bincounter++;
            output+="  ";
            for(int j=minbin;j<=i;j++){
                output += h->GetBinContent(i,j);
                output+=" ";
            }
            output+="\n";
        }
    }
    else{
        std::vector<TString> names;
        int bincounter=0;
        for(int i=minbin;i<maxbin;i++){
            TString n=prefix+"bin";
            n+=bincounter;
            bincounter++;
            names.push_back(n);
        }
        int bincounteri=0;
        int bincounterj=0;
        triangularMatrix m(names);
        for(int i=minbin;i<maxbin;i++){
            bincounterj=0;
            for(int j=minbin;j<=i;j++){
                m.setEntry(bincounteri,bincounterj,h->GetBinContent(i,j));
                bincounterj++;
            }
            bincounteri++;
        }
        double det=1;
        m.invert(det);
        if(!det){
            std::cout << "covariance could not be inverted. Please check for input errors or if the matrix is ill-defined." << std::endl;
            exit(-1);
        }
        for(size_t i=0;i<names.size();i++){
            output+=names.at(i)+"  ";
            for(size_t j=0;j<=i;j++){
                output+=m.getEntry(i,j);
                output+=" ";
            }
            output+="\n";
        }

    }


    return output;

}

int main(int argc, char* argv[]){


    TString histoname;
    TString outputname="covariance";
    TString prefix="";
    std::vector<TString> opts;


    for(int i=1;i<argc;i++){
        TString targv=argv[i];
        if(targv.BeginsWith("-")){

            if(targv.Contains("h")){
                coutHelp();
                exit(0);
            }
            if(targv == "-i"){
                if(i+1>=argc || ((TString)argv[i+1]).BeginsWith("-")){
                    std::cerr << "please specify a valid option argument" <<std::endl;
                    exit(-1);
                }
                TString next=argv[++i];
                if(next == "cov"){
                    outputname="covariance";
                }
                else if(next == "hes"){
                    outputname="hessian";
                }
                else if(next == "cor"){
                    outputname="correlation";
                }
                else{
                    std::cerr << "please specify a valid option argument" <<std::endl;
                    coutHelp();
                    exit(-1);
                }
            }
            else if(targv == "-p"){
                if(i+1>=argc || ((TString)argv[i+1]).BeginsWith("-")){
                    std::cerr << "please specify a valid option argument" <<std::endl;
                    exit(-1);
                }
                TString next=argv[++i];
                prefix=next;
            }

        }
        else{
            opts.push_back(targv);
        }
    }
    if(opts.size()!=2){
        coutHelp();
        exit(-1);
    }

    bool convertcovtohes=false;
    if(outputname=="covariance")
        convertcovtohes=true;
    outputname="hessian";

    TString all ="[";
    all+=outputname;
    all+="]\n\n";
    all += extractMatrix(opts.at(0),opts.at(1),prefix,convertcovtohes);
    all+="\n";
    all+="[end ";
    all+=outputname;
    all+="]\n";

    std::cout << all << std::endl;


}
