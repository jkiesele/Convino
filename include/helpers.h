/*
 * helpers.h
 *
 *  Created on: 1 Jul 2016
 *      Author: jkiesele
 */

#ifndef HELPERS_H_
#define HELPERS_H_

#include "TMatrixT.h"
#include "TVectorT.h"
#include <algorithm>
#include <iterator>
#include <sstream>
#include "textFormatter.h"
#include "TH2D.h"

/**
 * no overflow underflow
 * <bin> number is c++ like (starting at 0) not root like!
 */
TH2D removeOneBin(const TH2D& in, int bin);

inline double kr_delta(const size_t& i, const size_t& j){
	if(i==j)return 1.;
	return 0;
}


template <class T>
T round(T f,float pres)
{
	return (T) (floor(f*(1.0f/pres) + 0.5)/(1.0f/pres));
}


template<class t>
std::string toString(t in) {
	std::ostringstream s;
	s << in;
	std::string out = s.str();
	return out;
}


template<class t>
TString toTString(t in) {
	std::ostringstream s;
	s << in;
	TString out = s.str();
	return out;
}

template <class T>
std::ostream& operator<<(std::ostream& os,  const TMatrixT<T>& m)
{

	for(size_t i=0;i<(size_t)m.GetNrows();i++){
		for(size_t j=0;j<(size_t)m.GetNcols();j++){
			double entr=round(m[i][j],0.0001);
			os << std::left << textFormatter::fixLength(toString(entr), 6)<<" ";
		}
		os<<'\n';
	}

	return os;
}


template <class T>
std::ostream& operator<<(std::ostream& os,  const TVectorT<T>& m)
{
	std::streamsize save=os.width();
	os.width(4);
	std::streamsize savep=os.precision();
	os.precision(3);
	for(size_t i=0;i<m.GetNrows();i++){
		os.width(4);
		os.precision(2);
		os << std::left << toString(m[i])<<" ";
	}
	os<<'\n';

	os.width(save);
	os.precision(savep);
	return os;
}


template <class T>
std::ostream& operator<<(std::ostream& os,  const std::vector<T>& m)
{
	std::streamsize save=os.width();
	os.width(4);
	std::streamsize savep=os.precision();
	os.precision(3);
	for(size_t i=0;i<m.size();i++){
		os.width(4);
		os.precision(2);
		os << std::left << toString(m[i])<<" ";
	}
	os.width(save);
	os.precision(savep);
	return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os,  const std::vector< std::vector<T> >& m)
{
	os<<'\n';
	for(size_t i=0;i<m.size();i++){
		os << m.at(i) <<std::endl;
	}

	return os;
}



template <class T>
T FrobeniusDistance(TMatrixT<T> a, TMatrixT<T> b){
	if(a.GetNrows()!=b.GetNrows() || a.GetNcols() != b.GetNcols())
		return T(0);
	T out=0;
	for(int i=0;i<(int)a.GetNrows();i++){
		for(int j=0;j<(int)a.GetNcols();j++){
			out+=(a[i][j]-b[i][j])*(a[i][j]-b[i][j]);
		}
	}
	return sqrt(out);
}
template <class T>
T WeightedDistance(TMatrixT<T> a, TMatrixT<T> b, TMatrixT<T> weight){
	if(a.GetNrows()!=b.GetNrows() || a.GetNcols() != b.GetNcols())
		return T(0);
	T out=0;
	for(int i=0;i<(int)a.GetNrows();i++){
		for(int j=0;j<(int)a.GetNcols();j++){
			out+=(a[i][j]-b[i][j])*(a[i][j]-b[i][j])* weight[i][j]* weight[i][j];
		}
	}
	return sqrt(out);
}
template <class T>
T MaxDiff(TMatrixT<T> a, TMatrixT<T> b){
	if(a.GetNrows()!=b.GetNrows() || a.GetNcols() != b.GetNcols())
		return T(0);
	T out=0;
	for(int i=0;i<(int)a.GetNrows();i++){
		for(int j=0;j<(int)a.GetNcols();j++){
			T diff=fabs((a[i][j]-b[i][j]));
			if(diff>out)out=diff;
		}
	}
	return (out);
}


/**
 * works like std::sort but returns a vector of indecies that has the following form:
 * vector.at(unsrotedindex) == sortedindex
 * switched off in ROOT!!!
 */
template<typename _RandomAccessIterator, typename _Compare>
inline std::vector<size_t>
retsort(_RandomAccessIterator __first, _RandomAccessIterator __last,
		_Compare __comp);
/**
 * works like std::sort but returns a vector of indecies that has the following form:
 * vector.at(unsrotedindex) == sortedindex
 */
template<typename _RandomAccessIterator>
inline std::vector<size_t>
retsort(_RandomAccessIterator __first, _RandomAccessIterator __last);





template<typename _RandomAccessIterator, typename _Compare>
inline std::vector<size_t>
retsort(_RandomAccessIterator __first, _RandomAccessIterator __last,
		_Compare __comp){
	typedef typename std::iterator_traits<_RandomAccessIterator>::value_type
			_ValueType;
	std::vector<_ValueType> copy(__first,__last); //copy
	std::vector<size_t> sortedilo;
	std::sort(copy.begin(),copy.end(),__comp);

	for(_RandomAccessIterator it=copy.begin();it!=copy.end();++it){
		//get the position in input
		size_t pos=std::find(__first,__last,*it)-__first;
		sortedilo.push_back(pos);
	}
	std::copy(copy.begin(),copy.end(),__first);
	return sortedilo;
}


template<typename _RandomAccessIterator>
inline std::vector<size_t>
retsort(_RandomAccessIterator __first, _RandomAccessIterator __last){
	typedef typename std::iterator_traits<_RandomAccessIterator>::value_type
			_ValueType;
	std::vector<_ValueType> copy(__first,__last); //copy
	std::vector<size_t> sortedilo;
	std::sort(copy.begin(),copy.end());

	for(_RandomAccessIterator it=copy.begin();it!=copy.end();++it){
		//get the position in input
		size_t pos=std::find(__first,__last,*it)-__first;
		while(std::find(sortedilo.begin(),sortedilo.end(),pos)!=sortedilo.end())
			pos=std::find(__first+pos+1,__last,*it)-__first;
		sortedilo.push_back(pos);
	}
	std::copy(copy.begin(),copy.end(),__first);
	return sortedilo;
}


#endif /* CROSSSECTION_SUMMER16_HELPERS_H_ */
