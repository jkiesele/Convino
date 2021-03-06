/*
 * fileReader.h
 *
 *  Created on: Oct 10, 2013
 *      Author: kiesej
 */

#ifndef FILEREADER_H_
#define FILEREADER_H_


#include <string>
#include <vector>
#include <sstream>
#include "TString.h"
#include <stdexcept>
#include <iostream>
#include "textFormatter.h"



class fileReader : public textFormatter{
public:
	fileReader(): textFormatter(), start_(""),end_(""),blindmode_(false),requirevalues_(true){}
	~fileReader(){clear();}
	fileReader(const fileReader&rhs){
	    copyFrom(rhs);
	}
	fileReader & operator =(const fileReader& rhs){
	    copyFrom(rhs);
	    return *this;
	}

	/**
	 * ignores any white spaces!!
	 */
	void setStartMarker( std::string d);
	/**
		 * ignores any white spaces!!
		 */
	void setEndMarker( std::string d);

	void setBlindMode(bool blind){blindmode_=blind;}
	void setRequireValues(bool req) {requirevalues_=req;}

	/*
    std::string & ltrim(std::string & str) const;
    std::string & rtrim(std::string & str) const;
    std::string & trim(std::string & str) const;
    std::string & trimcomments(std::string & str) const;
	 */

	void readFile(const std::string& filename, bool allowempty=false);

	void addFromFile(const std::string& filename);


	// "[" are added automatically
	bool hasNonEmptyMarker(const std::string& filename, const std::string& startmarker, const std::string& endmarker)const;
	bool hasNonEmptyMarker(const std::string& filename, const std::string& startmarker)const;

	bool isEmpty()const{return lines_.size()<1;}
	size_t nLines()const{return lines_.size();}
	size_t nEntries(const size_t &line)const{return getData(line).size();}
	const std::vector<std::string>& getData(const size_t &line) const;
	template<class T>
	T getData(const size_t &line,const size_t &entry) const{
		if(line>=lines_.size()){
			throw std::out_of_range("fileReader::getData: line out of range");
		}
		if(entry>=lines_.at(line).size()){
			throw std::out_of_range("fileReader::getData: entry out of range");
		}
		T out;
		std::stringstream ss(getData(line).at(entry));
		ss >> out;
		return out;
	}

	std::string getReJoinedLine(size_t line) const;
	/**
	 * returns additional information connected to a marker formatted:
	 *
	 * [<markername> - <markervalue>] // additional spaces are allowed everywhere
	 *
	 */
	std::vector<std::string> getMarkerValues(const std::string& markername)const;

	std::string getValueString(const std::string &str, bool checkdoubles=true);

	template<class T>
	T getValue(const std::string & str){
		T out;
		std::string s(getValueString(str,true));
		if(s.length()<1){
			std::string exstring="fileReader::getValue: value for " + str + " not found in " + tempinfilename_;
			throw std::runtime_error(exstring);
		}
		std::stringstream ss(s);
		ss >> out;
		return out;
	}
	template<class T>
	T getValue(const std::string & str, T def_val){
		if(requirevalues_)
			return getValue<T>(str);
		T out;
		std::string s(getValueString(str,true));
		if(s.length()>0){
			std::stringstream ss(s);
			ss >> out;
			return out;
		}
		else{
			return def_val;
		}
	}

	void clear(){lines_.clear();}
	static bool debug;

	/**
	 * returns path to file
	 */
	std::string dumpFormattedToTmp()const;

	std::vector<std::string> readList(const std::string& infile,
			const std::string& startmarker, const std::string& endmarker,
			 const std::string& comment="", const std::string& delimiter="")const;

private:

	void copyFrom(const fileReader& rhs);

	enum read_return {read_success,read_nomarker,read_failed};

	read_return readFilePriv(const std::string& filename);

	std::string start_,end_;
	std::vector<std::vector<std::string> > lines_;
	std::vector<std::string> externals_;
	bool blindmode_;
	bool requirevalues_;
	std::string tempinfilename_;


	//make inherited private
	std::vector<std::string> getFormatted(const std::string& in)const{return std::vector<std::string>();}

};


template<>
inline bool fileReader::getData<bool>(const size_t &line,const size_t &entry) const{
	bool out;
	std::stringstream ss(getData(line).at(entry));
	ss >> std::boolalpha >> out;
	return out;
}
template<>
inline bool fileReader::getValue<bool>(const std::string & str){
	bool out;
	std::string s(getValueString(str,true));
	if(s.length()<1){
		std::string exstring="fileReader::getValue: value for " + str + " not found in " + tempinfilename_;
		throw std::runtime_error(exstring);
	}

	std::stringstream ss(s);
	ss >> std::boolalpha >> out;
	return out;
}
template<>
inline bool fileReader::getValue<bool>(const std::string & str, bool def_val){
	if(requirevalues_)
		return getValue<bool>(str);
	bool out;
	std::string s(getValueString(str,true));
	if(s.length()>0){
		std::stringstream ss(s);
		ss >> std::boolalpha >>out;
		return out;
	}
	else{
		return def_val;
	}
}
template<>
inline std::string fileReader::getData<std::string>(const size_t &line,const size_t &entry) const{
	return (getData(line).at(entry));
}
template<>
inline std::string fileReader::getValue<std::string>(const std::string & str){
	std::string s(getValueString(str));
	if(s.length()<1){
		std::string exstring="fileReader::getValue: value for " + str + " not found in " + tempinfilename_;
		throw std::runtime_error(exstring);
	}
	return s;
}
template<>
inline std::string fileReader::getValue<std::string>(const std::string & str,  std::string   def_val){
	if(requirevalues_)
		return getValue<std::string>(str);
	std::string s(getValueString(str));
	if(s.length()<1)
		return def_val;
	return s;
}
template<>
inline TString fileReader::getData<TString>(const size_t &line,const size_t &entry) const{
	return (TString)(getData(line).at(entry));
}
template<>
inline TString fileReader::getValue<TString>(const std::string & str){
	std::string s(getValueString(str));
	if(s.length()<1){
		std::string exstring="fileReader::getValue: value for " + str + " not found in " + tempinfilename_;
		throw std::runtime_error(exstring);
	}
	return (TString)s;
}
template<>
inline TString fileReader::getValue<TString>(const std::string & str,  TString  def_val){
	if(requirevalues_)
		return getValue<TString>(str);
	std::string s(getValueString(str));
	if(s.length()<1)
		return def_val;
	return (TString)s;
}





#endif /* FILEREADER_H_ */
