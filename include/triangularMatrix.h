/*
 * corrMatrix.h
 *
 *  Created on: Apr 13, 2015
 *      Author: kiesej
 */

#ifndef TRIANGULARMATRIX_H_
#define TRIANGULARMATRIX_H_

#include <vector>
#include "TString.h"
#include "indexMap.h"
#include "TMatrixD.h"
#include <iostream>
#include <string>
#include "TH2.h"

/*
 *
 * Use as text-level interface. For the actual calculations export to TMatrix
 * It contains an unambiguous mapping of a name (e.g. uncertainty) to a row / column
 *
 */
class correlationMatrix;

class triangularMatrix{
public:
	triangularMatrix(const size_t &, double def=0);
	triangularMatrix(const std::vector<TString>&, double def=0);
	triangularMatrix(const std::vector<std::vector<double> > &);
	triangularMatrix(){}

	bool operator ==(const triangularMatrix&)const;
	bool operator !=(const triangularMatrix&)const;

	triangularMatrix operator *(double)const;
    triangularMatrix operator +(const triangularMatrix&)const;

	virtual ~triangularMatrix(){}

	virtual const double& getEntry(const size_t & , const size_t &)const;
	void setEntry(const size_t &,const size_t &, const double &);

	void toTMatrix(TMatrixD& out)const;
	void importTMatrix(const TMatrixD& in);

	void fillFromTH2(const TH2D&);
    void createFromTH2(const TH2D&);

	const size_t& size()const{return names_.size();}

	void reOrder(const std::vector<TString>& neworder);

	const TString& getEntryName(const size_t &)const;
	const size_t& getEntryIndex(const TString&)const;
	//return SIZE_MAX if not found
	 size_t getEntryIndexUS(const TString&)const;
	 std::vector<TString> createNamesVector()const;

	void removeEntries(const std::vector<size_t>& idxs);
	void removeEntry(const size_t& idx);

	void invert(double& det);
	void invert();

	/**
	 * this creates a Hessian from a covariance,
	 * NOT assuming it to be positive definite,
	 * or trying to invert it. BUT assuming everything to
	 * be Gaussian etc.
	 */
	triangularMatrix createHessianFromCovariance()const;

	double determinant()const;


	void clear();

	void removeSmallEntries(double threshold=0.1);

	/**
	 * converts e.g. a covariance matrix to a correlation matrix
	 */
	void normalize();

	double maxDifference(const triangularMatrix& rhs)const;
	double FrobeniusDistance(const triangularMatrix& rhs)const;

	void getMinMaxEntry(double& min, double& max, bool nozeros)const;
	void getAbsMinMaxEntry(double& min, double& max, bool nozeros)const;
	void getAbsMinMaxDiagEntry(double& min, double& max, bool nozeros)const;

	/**
	 * reads matrices from files
	 * The format is
	 * <any text>
	 * [matrix]
	 *
	 * #comment can be anywhere
	 * <name0> 1
	 * <name1> rho 1
	 * <name2> rho rho 1
	 *
	 * [end matrix]
	 */
	void readFromFile(const std::string& in, const std::string& marker="matrix");

	/**
	 * appends a block matrix after the last row/column
	 */
	void append(const triangularMatrix&);

	/**
	 * copies entries only, names are untouched
	 */
	void copyEntries(const triangularMatrix&);

	/**
	 * gives back an index vector starting from the smallest entry.
	 */
	std::vector<std::pair <size_t,size_t> > getSortedList()const;

	void splitIntoBlocks(size_t splitIdxToRight,
			triangularMatrix& firstblock,
			triangularMatrix& secondblock,
			TMatrixD& offdiagonal)const;

	virtual size_t indepDim()const{return entries_.size();}

	void printToStream(std::ostream& os, bool texFormatting=false,
	        float mathbfthresh=-1, float precision=-1)const;

protected:

	static const double one_;
	std::vector<double> entries_;
	indexMap<TString> names_;
private:
	size_t linearize(const int &,const int &)const;
	std::pair <size_t,size_t> deLinearize(const size_t &)const;
};

class correlationMatrix: public triangularMatrix{
public:
	correlationMatrix(): triangularMatrix(){}
	correlationMatrix(const size_t &s, double def=0): triangularMatrix(s,def){
	}
	correlationMatrix(const std::vector<TString>&n, double def=0):triangularMatrix(n,def){}

	correlationMatrix& operator=(const triangularMatrix& r){
		triangularMatrix::operator = (r);
		//entries_=r.entries_;
		//names_=r.names_;
		return *this;
	}

	void setOffDiagonalZero();

	size_t indepDim()const{return entries_.size()-size();}

	void readFromFile(const std::string& in, const std::string& marker="matrix");

	const double& getEntry(const size_t & , const size_t &)const;
private:

};



std::ostream& operator<<(std::ostream& os, const triangularMatrix& m);



#endif /* CORRMATRIX_H_ */
