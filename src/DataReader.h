/*
 * DataReader.h
 *
 *  Created on: 5 mar. 2012
 *      Author: kopp
 */

#ifndef DATAREADER_H_
#define DATAREADER_H_

#include <StringTools.h>
#include <Log.h>

//#include <gsl/gsl_math.h>

#include <fstream>
#include <map>
#include <algorithm>

class DataReader
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m) :
				msg("DataReader::" + m){}
		~Exception() throw (){}
		const char* what() const throw (){return msg.c_str();}
	private:
		std::string msg;
	};
	typedef std::string ColumnName;
	typedef std::vector<double> ColumnType;
	typedef std::map<ColumnName, double> RowType;
	typedef std::pointer_to_unary_function<double, double> TransformType;
	class TCondition: public std::unary_function<const RowType&, bool>
	{
	public:
		TCondition() {}
		virtual ~TCondition() {}
		virtual bool operator()(const RowType& arg) const = 0;
	};
public:
	DataReader();
	virtual ~DataReader();
	void readFile(std::string filename);
	std::string getFilename() const {return filename;}

	void getColumn(DataReader::ColumnType& col, const ColumnName& name) const;
	void getColumn(DataReader::ColumnType& col, const ColumnName& name,
			const TransformType& apply) const;
	void filter(const TCondition * apply = NULL);
private:
	std::string filename;

	std::map<ColumnName, ColumnType> data;
	std::map<ColumnName, ColumnType> dataf;
	std::vector<ColumnName> columnNames;
	RowType row;
	size_t nbRows;

	bool isComment(const std::string& str);
	void readHeader(const std::string& str);
};

#endif /* DATAREADER_H_ */
