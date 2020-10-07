

#include <iostream>
#include "triangularMatrix.h"
#include "TString.h"
#include <fstream>

void coutHelp(){
    std::cout << "\nUSAGE: convino_text_to_pdf_matrix <OPTIONS> <input text file> <output tex file>\n";
    std::cout << "OPTIONS:\n";
    std::cout << "\n-h     display this help message\n";
    std::cout << "\n-i     identifier in the text file such as theMatrix, omitting the brackets\n";
    std::cout << "\n-t     threshold to plot a number in bold (default none)\n";
    std::cout << "\n-p     precision (default 0.01)\n";
    std::cout << "\n--cov  converts a correlation matrix with constraints into a covariance matrix\n";
    std::cout << std::endl;
    std::cout << "EXAMPLE: convino_text_to_pdf_matrix -i myMatrix result.txt myMatrix.tex" <<std::endl;
}




int main(int argc, char* argv[]){

    if(argc<2+1 || argc >9+1){
        coutHelp();
        return -1;
    }
    std::string id,infile,outfile;
    float mathbfthresh=-1, precision=0.01;
    bool converttocov=false;
    for(int i=1;i<argc;i++){
        TString targv=argv[i];
        if(targv.BeginsWith("-")){
            if(targv == "--cov"){
                converttocov=true;
                continue;
            }
            if(targv.Contains("h")){
                coutHelp();
                return 0;
            }
            if(targv.Contains("i") &&  i+1<argc){
                id = argv[i+1];
                i++;
                continue;
            }
            if(targv.Contains("t") &&  i+1<argc){
                mathbfthresh = atof(argv[i+1]);
                i++;
                continue;
            }
            if(targv.Contains("p") &&  i+1<argc){
                precision = atof(argv[i+1]);
                i++;
                continue;
            }
            continue;
        }
        if(!infile.length()){
            infile=targv;
            continue;
        }
        if(!outfile.length()){
            outfile=targv;
            continue;
        }
    }

    triangularMatrix m;
    m.readFromFile(infile,id,converttocov);

    std::cout << "writing to " << outfile << " precision " << precision <<  std::endl;
    std::ofstream out(outfile);

    m.printToStream(out, true,mathbfthresh,precision);

    out.close();

    return 0;
}
