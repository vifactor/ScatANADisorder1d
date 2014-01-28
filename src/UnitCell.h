/*
 * StructureXYZ.h
 *
 *  Created on: 2 nov. 2012
 *      Author: kopp
 */

#ifndef STRUCTUREXYZ_H_
#define STRUCTUREXYZ_H_

#include <StringTools.h>
#include <Log.h>
#include <fstream>

class UnitCellAtom
{
public:
	UnitCellAtom(std::string el = "", double x = 0.0, double y = 0.0, double z = 0.0, double p = 1.0, double s = 0)
	{
		_element = el;
		_x = x;
		_y = y;
		_z = z;
		_occupation = p;
		_sigmasq = s;
	}
	//coordinates
	double _x, _y, _z;
	//occupation probability
	double _occupation;
	//name
	std::string _element;
	//Debye-Waller scattering factor
	double _sigmasq;
};

class UnitCell : public std::vector<UnitCellAtom>
{
public:
	UnitCell(){}
	~UnitCell() {}
protected:
};

class UnitCellReader
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "::UnitCellReader" + m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};
	enum FileType {ftUC, ftXYZ, ftCIF};
	UnitCellReader() {}
	~UnitCellReader() {}
	UnitCell read(std::string filename, FileType ft = ftUC);
protected:
	bool isEmpty(const std::string& str);
	bool isUCComment(const std::string& str);
	UnitCell readFromUCFile(std::string filename);
	FileType defineFileType(std::string filename);
};

#endif /* STRUCTUREXYZ_H_ */
