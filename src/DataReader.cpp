/*
 * DataReader.cpp
 *
 *  Created on: 5 mar. 2012
 *      Author: kopp
 */

#include "DataReader.h"

DataReader::DataReader()
{
	nbRows = 0;
}

DataReader::~DataReader()
{
}

void DataReader::readFile(std::string filename)
{
	this->filename=filename;

	std::ifstream fin(filename.c_str());
	if (!fin)
		throw Exception("File has not been opened:\t" + filename);

	std::string line;
	std::istringstream is;
	size_t ln;
	std::vector<std::string> values;
	double value;

	//here we read a header line
	getline(fin, line);
	readHeader(line);

	ln = 1;
	while (!fin.eof())
	{
		getline(fin, line);
		++ln;
		//skip the comment
		if (isComment(line))
			continue;

		//split line into columns
		values = split(line, " \t\n \r");
		//check if nb of columns corresponds to number of headers
		if(values.size() != columnNames.size())
		{
			/*LOG(logWARNING) << "Line: " << ln << " of "<< filename << ":\n"
							<< "\t nb columns expected:\t" << columnNames.size() << "\n"
							<< "\t eb columns found:\t" << values.size();*/

			continue;
		}

		for(size_t icol = 0; icol < columnNames.size(); icol++)
		{
			is.str(values[icol]);
			is >> value;
			is.clear();

			data[columnNames[icol]].push_back(value);
			dataf[columnNames[icol]].push_back(value);

		}
		++nbRows;

	}
	fin.close();
}

bool DataReader::isComment(const std::string& str)
{
	//empty string
	if(str.empty()) return true;
	//comment
	if(str[0]=='#') return true;
	return false;
}

void DataReader::readHeader(const std::string& str)
{
	//no header line of type: % colName1 colName2...
	if(str[0] != '#')
		throw Exception("No header line:\t" + filename);

	//split header line into column names with delimiters "% \t\n"
	columnNames = split(str, "# \t\n\r");

	if(columnNames.empty())
		throw Exception("Empty header:\t" + filename);
}

void DataReader::getColumn(ColumnType& col, const ColumnName& name) const
{
	col.clear();

	std::map<ColumnName, ColumnType>::const_iterator it_map = dataf.find(name);
	if(it_map != dataf.end())
	{
		for(size_t i = 0; i < it_map->second.size(); ++i)
			col.push_back(it_map->second.at(i));
	}
	else
		throw Exception("Column \""+ name +"\" has not been read from" + filename);
}

void DataReader::getColumn(ColumnType& col, const ColumnName& name,
		const TransformType& apply) const
{
	getColumn(col, name);

	std::transform(col.begin(), col.end(), col.begin(), apply);
}

void DataReader::filter(const TCondition * apply)
{
	std::map<ColumnName, ColumnType>::const_iterator it_data;
	RowType::const_iterator it_row;

	if (apply == NULL)
	{
		dataf = data;
	}
	else if (!data.empty())
	{
		dataf.clear();
		for (size_t ip = 0; ip < nbRows; ip++)
		{
			for (it_data = data.begin(); it_data != data.end(); ++it_data)
			{
				row[it_data->first] = it_data->second[ip];
			}
			if (apply->operator ()(row))
			{
				for (it_row = row.begin(); it_row != row.end(); ++it_row)
				{
					dataf[it_row->first].push_back(it_row->second);
				}
			}
		}
	}
}
