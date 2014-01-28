/*
 * UnitCell.cpp
 *
 *  Created on: 2 nov. 2012
 *  Changed on: 22 jan. 2013
 *      Author: kopp
 */

#include "UnitCell.h"

UnitCell UnitCellReader::readFromUCFile(std::string filename)
{
	const size_t nbColumnsExpected = 6;
	std::ifstream fin;
	std::string line;
	std::istringstream is;
	std::vector<std::string> values;
	UnitCellAtom atom;
	UnitCell cell;
	size_t ln;

	ln = 0;
	fin.open(filename.c_str());
	while (!fin.eof())
	{
		//read line
		getline(fin, line);
		//increment line number
		++ln;
		//skip the comments and empty lines
		if (isEmpty(line))
			continue;
		if (isUCComment(line))
			continue;

		//split line into columns
		values = split(line, " \t\n");

		//check if nb of columns corresponds to 6 (name, x, y, z, occ, sigma^2)
		if(values.size() != nbColumnsExpected)
		{
			LOG(logWARNING) << "Line: " << ln << " of "<< filename << " has been skipped:\n"
							<< "\t nb columns expected:\t" << nbColumnsExpected << "\n"
							<< "\t nb columns found:\t" << values.size();

			continue;
		}

		is.str(line);
		is >> atom._element >> atom._x >> atom._y >> atom._z >> atom._occupation >> atom._sigmasq;
		is.clear();

		cell.push_back(atom);

	}
	fin.close();

	return cell;
}

UnitCell UnitCellReader::read(std::string filename, FileType ft)
{
	UnitCell cell;

	switch(ft)
	{
	case ftXYZ:
		break;
	case ftCIF:
		break;
	case ftUC:
	default:
		cell = readFromUCFile(filename);
		break;
	}

	return cell;
}

bool UnitCellReader::isUCComment(const std::string& str)
{
	//comment-line starts with #
	if(str[0]=='#') return true;
	return false;
}

bool UnitCellReader::isEmpty(const std::string& str)
{
	//comment starts with #
	if(str.empty()) return true;
	return false;
}
