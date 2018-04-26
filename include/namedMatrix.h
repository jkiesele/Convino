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
#include "helpers.h"

//small helper
//not as sophisticated as triang. Make base class at some point maybe..
template <class T>
class namedMatrix{
public:

	namedMatrix(const std::vector<TString>&,const std::vector<TString>&, T def=T());
	namedMatrix(){}

	size_t nRows()const{return namesi_.size();}
	size_t nCols()const{return namesj_.size();}

	const TString& getEntryNameRow(const size_t &)const;
	 size_t getEntryIndexUSRow(const TString&)const;

	const TString& getEntryNameCol(const size_t &)const;
	 size_t getEntryIndexUSCol(const TString&)const;


	void setEntry(const size_t & i,const size_t & j, const T & val){
		entries_.at(i).at(j)=val;
	}
	const T& getEntry(const size_t & i, const size_t & j)const{
		return entries_.at(i).at(j);
	}

	void addColumn(const TString& name, const std::vector<T> vals);

private:
	std::vector<std::vector<T> >entries_;
	indexMap<TString> namesi_,namesj_;
};

template<class T>
std::ostream& operator<<(std::ostream& os,  namedMatrix<T>& m);

template<class T>
namedMatrix<T>::namedMatrix(const std::vector<TString>& ni,const std::vector<TString>& nj, T def){
    entries_.resize(ni.size(),std::vector<T>(nj.size(),def));
    for(const auto& n:ni)
        namesi_.push_back(n);
    for(const auto& n:nj)
        namesj_.push_back(n);
}

template<class T>
const TString& namedMatrix<T>::getEntryNameRow(const size_t & idx)const{
    return namesi_.getData(idx);
}
//return SIZE_MAX if not found
template<class T>
size_t namedMatrix<T>::getEntryIndexUSRow(const TString&name)const{
    size_t idx=namesi_.getIndex(name);
    if(idx>=namesi_.size())
        return SIZE_MAX;
    return idx;
}
template<class T>
const TString& namedMatrix<T>::getEntryNameCol(const size_t &idx)const{
    return namesj_.getData(idx);
}

//return SIZE_MAX if not found
template<class T>
size_t namedMatrix<T>::getEntryIndexUSCol(const TString& name)const{
    size_t idx=namesj_.getIndex(name);
    if(idx>=namesj_.size())
        return SIZE_MAX;
    return idx;
}

template<class T>
void namedMatrix<T>::addColumn(const TString& name, const std::vector<T> vals){
    if(namesj_.getIndex(name) < namesj_.size())
        throw std::runtime_error("namedMatrix::addColumn: names must be unique");

    if(vals.size()!=namesi_.size())
        throw std::runtime_error("namedMatrix::addColumn: size does not match: "+toString(vals.size())+"/"+toString(namesi_.size()));

    for(size_t j=0;j<vals.size();j++){
        std::vector<T>& e=entries_.at(j);
        size_t oldsize=e.size();
        e.resize(e.size()+1);
        e.at(oldsize) = vals.at(j);
    }
    namesj_.push_back(name);

}

template<class T>
std::ostream& operator<<(std::ostream& os,  namedMatrix<T>& m)
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


#endif /* CONVINO_INCLUDE_NAMEDMATRIX_H_ */
