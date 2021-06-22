/*
 * corrMatrix.cc
 *
 *  Created on: Apr 13, 2015
 *      Author: kiesej
 */


#include "triangularMatrix.h"
#include "multiGausGenerator.h"
#include "fileReader.h"
#include "helpers.h"
#include "textFormatter.h"
#include <unistd.h>

const double triangularMatrix::one_=1.;
triangularMatrix::triangularMatrix(const size_t &tsszie, double def){
	for(size_t i=0;i<tsszie;i++){
		TString tmp=toTString(i);
		names_.push_back(tmp);
	}
	const size_t & N=names_.size();
	size_t size=(N+1)*N/2;
	entries_.resize(size,def);

}
triangularMatrix::triangularMatrix(const std::vector<TString>& names, double def){
	for(size_t i=0;i<names.size();i++){
		names_.push_back(names.at(i));
	}
	const size_t & N=names_.size();
	size_t size=(N+1)*N/2;
	entries_.resize(size,def);
}

triangularMatrix::triangularMatrix(const std::vector<std::vector<double> > & v){

	if(v.size()<1)
		return;
	size_t tsszie=v.size();
	if(v.size()!=v.at(0).size())
		throw std::runtime_error("triangularMatrix::triangularMatrix: input vector<vector> must be symmetric");
	for(size_t i=0;i<tsszie;i++){
		TString tmp=toTString(i);
		names_.push_back(tmp);
	}
	const size_t & N=names_.size();
	size_t size=(N+1)*N/2;
	entries_.resize(size);

	for(size_t i=0;i<tsszie;i++){
		for(size_t j=i;j<tsszie;j++)
			setEntry(i,j,v.at(i).at(j));
	}

}

bool triangularMatrix::operator ==(const triangularMatrix& rhs)const{
	return !(*this!=rhs);
}
bool triangularMatrix::operator !=(const triangularMatrix& rhs)const{
	if(entries_!=rhs.entries_) return true;
	if(names_!=rhs.names_) return true;
	return false;
}
triangularMatrix triangularMatrix::operator *(double scalar)const{
    auto out = *this;
    for(size_t i=0;i<size();i++){
        for(size_t j=0;j<size();j++){
            out.setEntry(i,j, scalar * getEntry(i,j));
        }
    }
    return out;
}


triangularMatrix triangularMatrix::operator +(const triangularMatrix& rhs)const{
    if(size() != rhs.size())
        throw std::runtime_error("triangularMatrix::operator + sizes don't match");
    auto out = *this;
    for(size_t i=0;i<size();i++){
        for(size_t j=0;j<size();j++){
            out.setEntry(i,j, getEntry(i,j) + rhs.getEntry(i,j));
        }
    }
    return out;
}


const TString& triangularMatrix::getEntryName(const size_t &idx)const{
	return names_.getData(idx);
}
const size_t& triangularMatrix::getEntryIndex(const TString& name)const{
	size_t idx=names_.getIndex(name);
	if(idx>=names_.size())
		throw std::out_of_range(("triangularMatrix::getEntryIndex: name \""+name+"\" not found").Data());
	return names_.getIndex(name);
}

size_t triangularMatrix::getEntryIndexUS(const TString& name)const{
	size_t idx=names_.getIndex(name);
	if(idx>=names_.size())
	    return SIZE_MAX;
	return idx;
}
std::vector<TString> triangularMatrix::createNamesVector()const{
    std::vector<TString> out(size());
    for(size_t i=0;i<size();i++){
        out[i]=names_.getData(i);
    }
    return out;
}

void triangularMatrix::removeEntries(const std::vector<size_t>& idxs){
	std::vector<TString> namesnew;
	for(size_t i=0;i<names_.size();i++){
		if(std::find(idxs.begin(),idxs.end(),i)!=idxs.end()) continue;
		namesnew.push_back(names_.at(i));
	}
	triangularMatrix out(namesnew);
	size_t ni=0,nj=0;
	for(size_t i=0;i<size();i++){
		nj=0;
		if(std::find(idxs.begin(),idxs.end(),i)!=idxs.end())  continue;
		for(size_t j=0;j<=i;j++){
			if(std::find(idxs.begin(),idxs.end(),j)!=idxs.end())  continue;
			out.setEntry(ni,nj,getEntry(i,j));
			nj++;
		}
		ni++;
	}
	*this=out;
}

void triangularMatrix::invert(double& det){
	TMatrixD m;
	toTMatrix(m);
	m.Invert(&det);
	if(!det)
		throw std::runtime_error("triangularMatrix::invert: determinant 0");
	importTMatrix(m);
}
void triangularMatrix::invert(){
	double det;
	invert(det);
}

triangularMatrix triangularMatrix::createHessianFromCovariance()const{
    if(determinant()>0){
        triangularMatrix cp=*this;
        cp.invert();
        return cp;
    }
    triangularMatrix cp=*this;
    //get sigmas
    std::vector<double> sigmasq(size());
    for(size_t i=0;i<size();i++)
        sigmasq.at(i)=getEntry(i,i);

    for(size_t i=0;i<size();i++){
        for(size_t j=0;j<=i;j++){
            cp.setEntry(i,j,getEntry(i,j)/(sigmasq.at(i)*sigmasq.at(j)));
        }
    }
    return cp;
}


triangularMatrix triangularMatrix::mergeCovarianceEntries(
        std::vector<std::pair<TString, std::vector<TString> > > which,
        const size_t stat_prec) const{

    //create a new covariance
    auto names = createNamesVector();
    std::vector<TString> newnames;
    //map the 1:1 associations
    std::map<size_t,std::vector<size_t> > asso_reduced_to_full;
    for(size_t i=0;i<names.size();i++){
        const auto& n = names.at(i);
        auto it = std::find_if(
                which.begin(), which.end(),
                [&n](const std::pair<TString, std::vector<TString> > & x)
                { for(const auto& ww:x.second)
                    if(n == ww)
                        return true;
                return false;});

        if(it == which.end()){
            asso_reduced_to_full[newnames.size()]={i};
            newnames.push_back(n);
        }
    }

    for (const auto &w : which){
        std::vector<size_t>  ids;
        for(const auto& sw: w.second){
            ids.push_back(getEntryIndex(sw));
        }
        asso_reduced_to_full[newnames.size()]=ids;
        newnames.push_back(w.first);
    }

    triangularMatrix out(newnames);
    multiGausGenerator gen;
    gen.setCovariance(*this);

    for(size_t toy=0;toy<stat_prec;toy++){

        auto toyvals =gen.generate();

        for(size_t i=0;i<out.size();i++){
            double xi=0;
            for(const auto& ii: asso_reduced_to_full[i])
                xi += toyvals.at(ii);

            for(size_t j=0;j<=i;j++){ //optimise
                double xj=0;
                for(const auto& jj: asso_reduced_to_full[j])
                    xj += toyvals.at(jj);

                double contrib = xi*xj/((double)(stat_prec));
                //if(i != j)
                //    contrib*=2;
                out.setEntry(i,j, out.getEntry(i,j) +  contrib);

            }

        }


    }
    return  out;
}

double triangularMatrix::determinant()const{
	TMatrixD m;
	toTMatrix(m);
	return m.Determinant();
}

void triangularMatrix::removeEntry(const size_t& idx){
	std::vector<size_t > idxs;
	idxs.push_back(idx);
	removeEntries(idxs);
	return;
	/*
	std::vector<TString> namesnew;
	for(size_t i=0;i<names_.size();i++){
		if(i==idx) continue;
		namesnew.push_back(names_.at(i));
	}
	triangularMatrix out(namesnew);
	size_t ni=0,nj=0;
	for(size_t i=0;i<size();i++){
		nj=0;
		if(i==idx) continue;
		for(size_t j=0;j<=i;j++){
			if(j == idx) continue;
			out.setEntry(ni,nj,getEntry(i,j));
			nj++;
		}
		ni++;
	}
	 *this=out;*/
}

void triangularMatrix::clear(){
	names_.clear();
	entries_.clear();
}
void triangularMatrix::removeSmallEntries(double threshold){

	for(size_t i=0;i<size();i++){
		for(size_t j=0;j<=i;j++){
			if(std::abs(getEntry(i,j))<=threshold)
				setEntry(i,j,0);
		}
	}


}

void triangularMatrix::normalize(){
	TMatrixD m2;
	toTMatrix(m2);
	TMatrixD diag(size(),size());
	for(int l=0;l<(int)size();l++){
		double entr=getEntry(l,l);
		if(entr>0)
			diag[l][l] = 1/sqrt(getEntry(l,l));
		else
			throw std::runtime_error("triangularMatrix::normalize: only for >0 diagonal elements");
	}
	m2*=diag;
	diag*=m2;
	m2=diag;
	importTMatrix(m2);
}

double triangularMatrix::maxDifference(const triangularMatrix& rhs)const{
	double out=0;
	for(size_t i=0;i<size();i++){
		for(size_t j=i;j<size();j++){
			double diff=fabs(getEntry(i,j)-rhs.getEntry(j,i));
			if(diff>out)out=diff;
		}
	}
	return out;
}
double triangularMatrix::FrobeniusDistance(const triangularMatrix& rhs)const{
	double out=0;
	for(size_t i=0;i<size();i++){
		for(size_t j=i;j<size();j++){
			const double& et=getEntry(i,j);
			const double& er=rhs.getEntry(j,i);
			if(i!=j)
				out+=2* (et-er)*(et-er);
			else
				out+=(et-er)*(et-er);
		}
	}
	return sqrt(out);
}
void triangularMatrix::getMinMaxEntry(double& min, double& max, bool nozeros)const{
	if(size()<1){
		min=0;max=0;return;
	}
	min=entries_.at(0);
	max=entries_.at(0);
	for(const auto& e:entries_){
		if(nozeros&&!e)continue;
		if(e<min)min=e;
		if(e>max)max=e;
	}
}
void triangularMatrix::getAbsMinMaxEntry(double& min, double& max, bool nozeros)const{
	if(size()<1){
		min=0;max=0;return;
	}
	min=std::abs(entries_.at(0));
	max=std::abs(entries_.at(0));
	for(const auto& e:entries_){
		double abse=std::abs(e);
		if(nozeros&&!abse)continue;
		if(abse<min)min=abse;
		if(abse>max)max=abse;
	}
}
void triangularMatrix::getAbsMinMaxDiagEntry(double& min, double& max, bool nozeros)const{
	if(size()<1){
		min=0;max=0;return;
	}
	min=std::abs(entries_.at(0));
	max=std::abs(entries_.at(0));
	for(size_t i=0;i<size();i++){
		double en=getEntry(i,i);
		if(nozeros && ! en)continue;
		if(en<min)min=en;
		if(en>max)max=en;
	}
}
void triangularMatrix::readFromFile(const std::string& in, const std::string& marker, bool convertToCov){
	fileReader fr;
	fr.setComment("#");
	fr.setDelimiter(" ");
	fr.setStartMarker((std::string)"["+marker+"]");
	fr.setEndMarker((std::string)"[end "+marker+"]");
	fr.readFile(in);
	size_t lastlinelength=1;
	size_t matrixoffset=0;
	bool check_triangular=true;
	bool iscorrtocov=false;
	if(fr.nLines() && (fr.nLines() == fr.nEntries(0)-1 ||  fr.nLines() == fr.nEntries(0)-2))
	    check_triangular=false;



	std::vector<TString> names;
	for(size_t i=0;i<fr.nLines();i++){
		size_t thislinelength=fr.nEntries(i);
		if(!thislinelength)
		    continue;
		if(!i && thislinelength > 2){ //check once: this could be additional correlation info
			TString entry=fr.getData<TString>(i,1);
			if(entry.BeginsWith("(") && entry.EndsWith(")")){
				//this is a constraint, handle later but allow here
				lastlinelength++;
				matrixoffset++;
				iscorrtocov=true;
			}
		}
		//either triangular format or square
		if(thislinelength!=lastlinelength+1 && check_triangular){
			std::string entrystr="";
			if(fr.nEntries(i)>0)
				entrystr=fr.getData<std::string>(i,0);
			entrystr="triangularMatrix::readFromFile: format wrong (see docu) "+entrystr;
			entrystr+=". expected ";
			entrystr+=toString(lastlinelength+1)+" found " + toString(thislinelength);
			throw std::runtime_error(entrystr);
		}
		lastlinelength=thislinelength;
		names.push_back(fr.getData<TString>(i,0));
	}
	*this=triangularMatrix(names);
	size_t im=0;
	for(size_t i=0;i<fr.nLines();i++){
        if(!fr.nEntries(i))
            continue;
		for(size_t j=1;j<fr.nEntries(i)-matrixoffset;j++){
			setEntry(i,j-1,fr.getData<double>(im,j+matrixoffset));
		}
		im++;
	}

	if(convertToCov && iscorrtocov){

	    std::vector<double> constr;
	    for(size_t i=0;i<fr.nLines();i++){
	        size_t nentries = fr.nEntries(i);
	        if(nentries>1 && nentries < 3)
	            continue; //allow for empty lines
	        TString constraintstr = fr.getData<TString>(i,1);

	        constraintstr.ReplaceAll("(","");
	        constraintstr.ReplaceAll(")","");
	        double constraint = constraintstr.Atof();

	        constr.push_back(constraint);
	    }

	    if(constr.size() != size())
	        throw std::runtime_error("triangularMatrix::readFromFile: text file format has issues, please check.");

	    auto corr=*this;
	    for(size_t i=0;i<corr.size();i++){
	        for(size_t j=i;j<corr.size();j++){
	            corr.setEntry(i,j, getEntry(i,j) * constr.at(i) * constr.at(j));
	        }
	    }
	    *this=corr;
	}

}

size_t triangularMatrix::linearize(const int & col,const int & row)const{
	size_t idx=0;
	if(col<=row)
		idx= col + (row+1)*row/2;
	else
		idx= row + (col+1)*col/2;
	if(idx>=entries_.size())
		throw std::out_of_range("triangularMatrix::linearize: row or col out of range");
	return idx;
}

std::pair <size_t,size_t> triangularMatrix::deLinearize(const size_t & in)const{
	if(in>=entries_.size())
		throw std::out_of_range("triangularMatrix::deLinearize");
	for(size_t idxa=0;idxa<size();idxa++){
		for(size_t idxb=0;idxb<size();idxb++){
			if(in == linearize(idxa,idxb))
				return std::pair <size_t,size_t> (idxa,idxb);
		}
	}
	throw std::out_of_range("triangularMatrix::deLinearize");//never happens
}

const double& triangularMatrix::getEntry(const size_t & row, const size_t & col)const{
	if(col > entries_.size() || row> entries_.size()){
		throw std::out_of_range("corrMatrix::getEntry: index out of range");
	}
	return entries_.at(linearize(row,col));
}


void  triangularMatrix::setEntry(const size_t & row,const size_t & col,const double & entr){
	if(col > entries_.size() || row> entries_.size()){
		throw std::out_of_range("corrMatrix::setEntry: index out of range");
	}
	entries_.at(linearize(row,col))=entr;
}
void triangularMatrix::toTMatrix(TMatrixD& out)const{

	out.ResizeTo(size(),size());
	for(size_t i=0;i<size();i++){
		for(size_t j=0;j<size();j++){
			out[i][j]=getEntry(i,j);
		}
	}
}

void triangularMatrix::importTMatrix(const TMatrixD& in){
	if(in.GetNcols() != in.GetNrows()
			|| in.GetNcols() != (int)size())
		throw std::out_of_range("triangularMatrix::importTMatrix: sizes don't match");

	for(size_t i=0;i<size();i++){
		for(size_t j=0;j<size();j++){
			setEntry(i,j,in[i][j]);
		}
	}

}

void triangularMatrix::fillFromTH2(const TH2D& h){
    if(h.GetNbinsX() != h.GetNbinsY())
        throw std::out_of_range("triangularMatrix::fillFromTH2: input not symmetric");
    size_t start=1,end=h.GetNbinsX()+1;

    if(end-start != size())
        throw std::out_of_range(((TString)"triangularMatrix::fillFromTH2: input size "+ (end-start) + (TString)" does not match matrix size "+
                size()).Data());

    for(size_t i=0;i<size();i++){
        for(size_t j=i;j<size();j++){
            setEntry(i,j,h.GetBinContent(i+start,j+start));
        }
    }
}


void triangularMatrix::reOrder(const std::vector<TString>& neworder){
	bool failed=false;
	if(neworder.size()!=names_.size()){
		failed=true;
	}
	if(!failed){
		size_t count=0;
		for(size_t i=0;i<names_.size();i++){
			if(std::find(neworder.begin(),neworder.end(), names_.at(i))!=neworder.end())
				count++;
		}
		if(count!=names_.size())
			failed=true;
	}
	if(failed)
		throw std::runtime_error("triangularMatrix::reOrder: input wrong");

	triangularMatrix newm(neworder);
	for(size_t i=0;i<newm.size();i++){
		size_t oi=getEntryIndex(newm.getEntryName(i));
		for(size_t j=i;j<newm.size();j++){
			size_t oj=getEntryIndex(newm.getEntryName(j));
			newm.setEntry(i,j, getEntry(oi,oj));
		}
	}
	*this=newm;
}


void triangularMatrix::append(const triangularMatrix& rhs){
	size_t oldsize=size();
	for(size_t i=0;i<rhs.names_.size();i++){
		size_t oldnamess=names_.size();
		if(names_.push_back(rhs.names_.getData(i))<oldnamess)
			throw std::logic_error("triangularMatrix::append: try to append at least one entry with same name");
	}
	std::vector<double> newentries=entries_;
	const size_t & N=names_.size();
	size_t newsize=(N+1)*N/2;
	entries_.resize(newsize,0);
	//entries_.insert(entries_.end(),rhs.entries_.begin(),rhs.entries_.end());
	for(size_t i=0;i<size();i++){
		for(size_t j=oldsize;j<size();j++){
			if(i>j)continue; //symmetric
			if(i<oldsize){ setEntry(i,j,0.);continue;}
			setEntry(i,j, rhs.getEntry(i-oldsize,j-oldsize));

		}
	}
}

void triangularMatrix::copyEntries(const triangularMatrix&rhs){
	if(rhs.size()!=size())
		throw std::out_of_range("triangularMatrix::copyEntries: sizes don't match");

	for(size_t i=0;i<size();i++){
		for(size_t j=i;j<size();j++){
			setEntry(i,j,rhs.getEntry(i,j));
		}
	}
}


std::vector<std::pair <size_t,size_t> > triangularMatrix::getSortedList()const{
	//find smalles value;
	if(size()<1) return std::vector<std::pair <size_t,size_t> >();

	std::vector<std::pair <size_t,size_t> > allidxs(entries_.size());

	for(size_t i=0;i<size();i++){
		std::cout << i << std::endl;
		for(size_t j=i;j<size();j++){
			allidxs.at(linearize(i,j))=std::pair<size_t,size_t>(i,j);
		}
	}

	std::vector<double> serntries=entries_;
	std::vector<size_t> sidxs=retsort(serntries.begin(),serntries.end());
	std::vector<std::pair <size_t,size_t> > out(entries_.size());
	for(size_t i=0;i<sidxs.size();i++){
		std::cout << i << std::endl;
		out.at(i)=allidxs.at(sidxs.at(i));
	}
	return out;
}

void triangularMatrix::splitIntoBlocks(size_t splitIdxToRight,
		triangularMatrix& firstblock,
		triangularMatrix& secondblock,
		TMatrixD& offdiagonal)const{

	if(&firstblock == this || &secondblock == this)
		throw std::logic_error("triangularMatrix::splitIntoBlocks cannot use itself as reference");

	std::vector<TString> namesfirst(splitIdxToRight);
	std::vector<TString> namessec  (size()-splitIdxToRight);

	offdiagonal.ResizeTo(size()-splitIdxToRight,splitIdxToRight);

	for(size_t i=0;i<size();i++){
		if(i<splitIdxToRight)
			namesfirst.at(i) = getEntryName(i);
		else
			namessec.at(i-splitIdxToRight) = getEntryName(i);
	}
	firstblock=triangularMatrix(namesfirst);
	secondblock=triangularMatrix(namessec);
	for(size_t i=0;i<size();i++){
		for(size_t j=i;j<size();j++){
			if(i<splitIdxToRight && j<splitIdxToRight){
				firstblock.setEntry(i,j, getEntry(i,j));
			}
			else if(i>=splitIdxToRight && j>=splitIdxToRight){
				secondblock.setEntry(i-splitIdxToRight,j-splitIdxToRight, getEntry(i,j) );
			}
			else{//off diagonal
				offdiagonal[j-splitIdxToRight][i]=getEntry(i,j);
			}
		}
	}


}






void triangularMatrix::printToStream(std::ostream& os, bool texFormatting,
        float mathbfthresh, float precision)const{
    double min,max;
    getMinMaxEntry(min,max,false);
    double delta = max-min;
    double scaler=1;
    while(!texFormatting && delta && delta*scaler < 1)scaler*=10;
    while(!texFormatting && delta && delta*scaler > 10)scaler*=0.1;

    std::string separator=" ";
    std::string newline="\n";
    if(texFormatting){
        separator=" & ";
        newline="\\\\ \\hline \n";

        //print tex header
        os << "\\begin{tabular}{| l | ";
        for(size_t i=0;i<size();i++)
            os << " c |";
        os << "}\n\\hline\n";

    }

    if(scaler <= 0.01 || scaler >= 100){
        os << "multiplied by: " << scaler <<'\n';}
    else{
        scaler=1;}

    std::streamsize save=os.width();
    if(!texFormatting)
        os.width(4);
    size_t maxnamewidth=0;
    for(size_t i=0;i<size();i++){
        size_t s=getEntryName(i).Length();
        if(s>maxnamewidth)maxnamewidth=s;
    }
    int number_length=4;
    int no_exp=0;
    if(precision>0){
        float preccp=precision;
        number_length=0;
        while(preccp <= (float)1){
            number_length++;
            no_exp++;
            preccp*=10;
        }
        if(!number_length)
            while(preccp > (float)1){
                number_length++;
                no_exp++;
                preccp/=10;
            }
        number_length++;
    }


    for(size_t i=0;i<size();i++){
        os.width(maxnamewidth+1);
        auto ename = getEntryName(i);
        if(texFormatting)
            ename = textFormatter::makeTexCompatible(ename.Data());
        os << std::left << ename;
        if(texFormatting)
            os << " & ";
        for(size_t j=0;j<size();j++){
            if(texFormatting){
                float entry=getEntry(i,j);
                if(precision>0)
                    entry = round(entry,precision);
                bool plotbf=false;
                if(mathbfthresh>0 && fabs(entry)>mathbfthresh)
                    plotbf=true;
                os << textFormatter::toScientificTex(entry,number_length,no_exp,true,plotbf);
            }
            else{
                double entrd=round(getEntry(i,j)*scaler,0.000001);
                os << textFormatter::fixLength(entrd,9);
            }
            if(j<size()-1)
                os << separator;
        }
        if(i<size()-1 || texFormatting)
            os<<newline;
    }
    if(texFormatting)
        os << "\\end{tabular}\n";
    os.width(save);
}




void correlationMatrix::setOffDiagonalZero(){
	entries_ = std::vector<double> (entries_.size(),0);
}

void correlationMatrix::readFromFile(const std::string& in, const std::string& marker){
	triangularMatrix::readFromFile(in, marker);

	//check entries if they make sense
	for(size_t i=0;i<entries_.size();i++)
		if(fabs(entries_.at(i)) > 1){
			clear();
			throw std::runtime_error("correlationMatrix::readFromFile: at least one entry > +-1 please check input.");
		}
}

const double& correlationMatrix::getEntry(const size_t & r, const size_t & c)const{
	if(r==c)return one_;
	return triangularMatrix::getEntry(r,c);
}



std::ostream& operator<<(std::ostream& os, const triangularMatrix& m)
{
    m.printToStream(os,false);
    return os;
}
