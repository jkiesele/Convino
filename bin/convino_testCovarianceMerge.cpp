
#include "triangularMatrix.h"



int main(){

    triangularMatrix m;
    m.readFromFile("examples/testCovariance.txt","covariance");

    std::vector<std::pair<TString, std::vector<TString> > > which=
    { { "Msys1" , {"systA", "systC", "systD"} },
            { "Msys2" , {"systB", "systF"} }
    };


    auto merged = m.mergeCovarianceEntries(which);


    std::cout << "original" <<std::endl;
    std::cout << m << std::endl;

    std::cout << "merged" <<std::endl;
    std::cout << merged << std::endl;

    return 0;
}
