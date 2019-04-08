/*
 * gentest.cpp
 *
 *  Created on: 5 Apr 2019
 *      Author: jkiesele
 */


#include "multiGausGenerator.h"
#include "triangularMatrix.h"
#include "math.h"
#include "TH1D.h"
#include "textFormatter.h"
#include "helpers.h"
#include "normaliser.h"

int nbins=5;

triangularMatrix createCovarianceMatrix(const TH1D& nominal, double correlation_strength);
TH1D createPseudoMeasurement(TString name, bool withUFOF, int statistics, bool normalise=false, bool flat=false);


std::ostream& operator<<(std::ostream& os, const std::vector<double>& m);
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double> >& m);


int main(){

    TH1D meas = createPseudoMeasurement("histo", false, 2000, false,true);

    triangularMatrix m = createCovarianceMatrix(meas,1);

    std::cout << m << std::endl;







    return 1;


    multiGausGenerator gen;

    gen.setCovariance(m);

    std::vector<std::vector<double> > covfilled(m.size(),std::vector<double>(m.size(),0));

    int igen=5000000;

    for(int i=0;i<igen;i++){
        auto vec = gen.generate();

       // std::cout << vec << std::endl;

        for(size_t i=0;i<m.size();i++){
            for(size_t j=0;j<m.size();j++){
                covfilled[i][j] += (vec[i]-0.)*(vec[j]-0.); //expectation value is zero
            }
        }

    }
    for(auto& i:covfilled)
        for(auto& j:i)
            j/=(double)igen;

std::cout << "cov\n" << covfilled << std::endl;

}








triangularMatrix createCovarianceMatrix(const TH1D& nominal, double correlation_strength){


    triangularMatrix m(nbins);

    for(size_t i=0; i<m.size();i++){
        double stati=nominal.GetBinError(i+1);
        for(int j=0; j<=i;j++){
            double statj=nominal.GetBinError(j+1);
            double correlation = exp(-1/correlation_strength*(i-j)*(i-j));
            m.setEntry(i,j, stati*statj*correlation );
        }
    }
    return m;
}

TH1D createPseudoMeasurement(TString name, bool withUFOF, int statistics, bool normalise, bool flat){

    if(withUFOF) //in this function done
        throw std::runtime_error("inclusion of underflow and overflow is not implemented fully consistently everywhere yet. Please make sure your input does not include those.");

    //if(withUFOF)
     //   nbins+=2;
    TH1D h(name,name,nbins,-2,2);

    h.FillRandom("gaus",statistics);
    if(flat)
        for(int i=0;i<=h.GetNbinsX();i++)
            h.SetBinContent(i, (double)statistics/(double)nbins);

    double integ=h.Integral();
    if(normalise)
        h.Scale(1/integ);

    return h;
}



std::ostream& operator<<(std::ostream& os, const std::vector<double>& m){
    std::streamsize save=os.width();
    os.width(4);

    for(int j=0;j<m.size();j++){
        double entrd=round(m[j] ,0.0001);
        std::string entr=toString(entrd);
        os << textFormatter::fixLength(entr,6);
        os <<" ";
    }

    os.width(save);
    return os;

}

std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double> >& m){
    std::streamsize save=os.width();
    os.width(4);

    for(size_t i=0;i<m.size();i++){
        for(int j=0;j<m.at(i).size();j++){
            double entrd=round(m[i][j] ,0.0001);
            std::string entr=toString(entrd);
            os << textFormatter::fixLength(entr,6);
            os <<" ";
        }
        os << '\n';
    }

    os.width(save);
    return os;

}
