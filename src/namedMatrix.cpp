/*
 * matrix.cpp
 *
 *  Created on: 25 Jan 2017
 *      Author: jkiesele
 */


#include "../include/namedMatrix.h"
#include "helpers.h"

namedMatrix::namedMatrix(const std::vector<TString>& ni,const std::vector<TString>& nj, double def){
	entries_.resize(ni.size(),std::vector<double>(nj.size(),def));
	for(const auto& n:ni)
		namesi_.push_back(n);
	for(const auto& n:nj)
		namesj_.push_back(n);
}


const TString& namedMatrix::getEntryNameRow(const size_t & idx)const{
	return namesi_.getData(idx);
}
//return SIZE_MAX if not found
size_t namedMatrix::getEntryIndexUSRow(const TString&name)const{
	size_t idx=namesi_.getIndex(name);
	if(idx>=namesi_.size())
		return SIZE_MAX;
	return idx;
}

const TString& namedMatrix::getEntryNameCol(const size_t &idx)const{
	return namesj_.getData(idx);
}

//return SIZE_MAX if not found
size_t namedMatrix::getEntryIndexUSCol(const TString& name)const{
	size_t idx=namesj_.getIndex(name);
	if(idx>=namesj_.size())
		return SIZE_MAX;
	return idx;
}


void namedMatrix::addColumn(const TString& name, const std::vector<double> vals){
	if(namesj_.getIndex(name) < namesj_.size())
		throw std::runtime_error("namedMatrix::addColumn: names must be unique");

	if(vals.size()!=namesi_.size())
		throw std::runtime_error("namedMatrix::addColumn: size does not match: "+toString(vals.size())+"/"+toString(namesi_.size()));

	for(size_t j=0;j<vals.size();j++){
		std::vector<double>& e=entries_.at(j);
		size_t oldsize=e.size();
		e.resize(e.size()+1);
		e.at(oldsize) = vals.at(j);
	}
	namesj_.push_back(name);

}

std::ostream& operator<<(std::ostream& os,  namedMatrix& m)
{
	os<<"\n             ";
	for(size_t i=0;i<m.nCols();i++)
		os << m.getEntryNameCol(i) << "  ";


	for(size_t i=0;i<m.nRows();i++){
		os <<'\n'<< m.getEntryNameRow(i) << "   ";
		for(size_t j=0;j<m.nCols();j++){
			os << m.getEntry(i,j) << " ";
		}
	}

	return os;
}
