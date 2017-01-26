/*
 * matrix.h
 *
 *  Created on: 25 Jan 2017
 *      Author: jkiesele
 */

#ifndef CONVINO_INCLUDE_NAMEDMATRIX_H_
#define CONVINO_INCLUDE_NAMEDMATRIX_H_

#include "indexMap.h"
#include <vector>
#include "TString.h"

//small helper
//not as sophisticated as triang. Make base class at some point maybe..
class namedMatrix{
public:

	namedMatrix(const std::vector<TString>&,const std::vector<TString>&, double def=0);
	namedMatrix(){}

	size_t nRows()const{return namesi_.size();}
	size_t nCols()const{return namesj_.size();}

	const TString& getEntryNameRow(const size_t &)const;
	 size_t getEntryIndexUSRow(const TString&)const;

	const TString& getEntryNameCol(const size_t &)const;
	 size_t getEntryIndexUSCol(const TString&)const;


	void setEntry(const size_t & i,const size_t & j, const double & val){
		entries_.at(i).at(j)=val;
	}
	const double& getEntry(const size_t & i, const size_t & j)const{
		return entries_.at(i).at(j);
	}

	void addColumn(const TString& name, const std::vector<double> vals);

private:
	std::vector<std::vector<double> >entries_;
	indexMap<TString> namesi_,namesj_;
};


std::ostream& operator<<(std::ostream& os,  namedMatrix& m);


#endif /* CONVINO_INCLUDE_NAMEDMATRIX_H_ */
