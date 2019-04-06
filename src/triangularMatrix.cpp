/*
 * corrMatrix.cc
 *
 *  Created on: Apr 13, 2015
 *      Author: kiesej
 */


#include "triangularMatrix.h"
#include "fileReader.h"
#include "helpers.h"
#include "textFormatter.h"

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
void triangularMatrix::readFromFile(const std::string& in, const std::string& marker){
	fileReader fr;
	fr.setComment("#");
	fr.setDelimiter(" ");
	fr.setStartMarker((std::string)"["+marker+"]");
	fr.setEndMarker((std::string)"[end "+marker+"]");
	fr.readFile(in);
	size_t lastlinelength=1;
	size_t matrixoffset=0;
	std::vector<TString> names;
	for(size_t i=0;i<fr.nLines();i++){
		size_t thislinelength=fr.nEntries(i);
		if(!i && thislinelength > 2){ //check once: this could be additional correlation info
			TString entry=fr.getData<TString>(i,1);
			if(entry.BeginsWith("(") && entry.EndsWith(")")){
				//this is a constraint, handle later but allow here

				lastlinelength++;
				matrixoffset++;
			}
		}
		if(thislinelength!=lastlinelength+1){
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
	for(size_t i=0;i<fr.nLines();i++){
		for(size_t j=1;j<fr.nEntries(i)-matrixoffset;j++){
			setEntry(i,j-1,fr.getData<double>(i,j+matrixoffset));
		}
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

void triangularMatrix::fillFromTH2(const TH2D& h, bool includeufof){
    if(h.GetNbinsX() != h.GetNbinsY())
        throw std::out_of_range("triangularMatrix::fillFromTH2: input not symmetric");
    size_t start=1,end=h.GetNbinsX()+1;
    if(includeufof){
        start=0;
        end++;
    }
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
    double min,max;
    m.getMinMaxEntry(min,max,false);
    double delta = max-min;
    double scaler=1;
    while(delta && delta*scaler < 1)scaler*=10;
    while(delta && delta*scaler > 10)scaler*=0.1;

    if(scaler <= 0.01 || scaler >= 100)
        os << "multiplied by: " << scaler <<'\n';
    else
        scaler=1;
	std::streamsize save=os.width();
	os.width(4);
	size_t maxnamewidth=0;
	for(size_t i=0;i<m.size();i++){
		size_t s=m.getEntryName(i).Length();
		if(s>maxnamewidth)maxnamewidth=s;
	}
	for(size_t i=0;i<m.size();i++){
		os.width(maxnamewidth+1);
		os << std::left << m.getEntryName(i);
		for(size_t j=0;j<m.size();j++){
			double entrd=round(m.getEntry(i,j)*scaler,0.001);
			std::string entr=toString(entrd);
			os << textFormatter::fixLength(entr,6);
			os <<" ";
		}
		os<<'\n';
	}

	os.width(save);
	return os;
}
