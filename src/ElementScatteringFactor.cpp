/*
 * ElementScatteringFactor.cpp
 *
 *  Created on: 1 nov. 2012
 *      Author: kopp
 */

#include "ElementScatteringFactor.h"

ElementScatteringFactor::ElementScatteringFactor(std::string elname, SFCorrectionDataBaseEnum dbf1, SFDataBaseEnum dbf0)
{
	elementName = elname;

	readF0DB(dbf0);

	readF1DB(dbf1);
}

ElementScatteringFactor::~ElementScatteringFactor()
{
}

double ElementScatteringFactor::f0(double k) const
{
	static double res;

	res = c;
	for (size_t i = 0; i < NB_COEFS; ++i)
	{
		res += a[i] * exp(-b[i] * k * k);
		LOG(logDEBUG) << i << "\t" << a[i] << "\t" << b[i];
	}
	return res;
}

double ElementScatteringFactor::ref1(double en) const
{
	static std::vector<Correction>::const_iterator low,up;
	static Correction corr;

	corr.en = en;

	low = std::lower_bound(corrections.begin(), corrections.end(), corr, cmpObj);
	up = std::upper_bound(corrections.begin(), corrections.end(), corr, cmpObj);

	LOG(logDEBUG) << "boundaries:\t" << (low - 1)->en << " < " << en << " < " << up->en;
	LOG(logDEBUG) << "values:\t" << (low - 1)->ref1 << " < ref1 < " << up->ref1;

	return ((low - 1)->ref1 + up->ref1)/2;
}

double ElementScatteringFactor::imf1(double en) const
{
	static std::vector<Correction>::const_iterator low,up;
	static Correction corr;

	corr.en = en;

	low = std::lower_bound(corrections.begin(), corrections.end(), corr, cmpObj);
	up = std::upper_bound(corrections.begin(), corrections.end(), corr, cmpObj);

	LOG(logDEBUG) << "boundaries:\t" << (low - 1)->en << " < " << en << " < " << up->en;
	LOG(logDEBUG) << "values:\t" << (low - 1)->imf1 << " < imf1 < " << up->imf1;

	return ((low - 1)->imf1 + up->imf1) / 2;
}

void ElementScatteringFactor::readF0WaasmaierKirfel()
{
	std::string filename = "f0_WaasKirf.dat";
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	fin.open(filename.c_str());
	fin.exceptions( std::ifstream::failbit | std::ifstream::badbit);

	//read data
	while (true)
	{
		if(fin.eof())
		{
			LOG(logERROR) << "Cannot find " << elementName << " in " << filename << std::endl;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n\r");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(elementName.compare(elname) == 0)
			{
				//skip 2 lines
				getline(fin, line);
				getline(fin, line);

				//read atomic parameters
				getline(fin, line);
				//std::cout<<"Text(data):\t"<<text<<std::endl;

				std::istringstream is(line);
				is >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> c >> b[0] >> b[1]
						>> b[2] >> b[3] >> b[4];
				LOG(logINFO) << "Parameters for " << elementName << ":\t"
						<< c << "\t" << a[0] << "\t" << a[1] << "\t" << a[2]
						<< "\t" << a[3] << "\t" << a[4] << "\t" << b[0] << "\t"
						<< b[1] << "\t" << b[2] << "\t" << b[3] << "\t" << b[4] << std::endl;
				return;
			}
		}
	}
	fin.close();
}

void ElementScatteringFactor::readF0DB(SFDataBaseEnum db)
{
	switch(db)
	{
	case f0WAASMAIER_KIRFEL:
		readF0WaasmaierKirfel();
		break;
	default:
		//should never appear
		break;
	}
}

void ElementScatteringFactor::readF1DB(SFCorrectionDataBaseEnum db)
{
	switch(db)
	{
	case f1BRENNAN_COWAN:
		readF1BrennanCowan();
		break;
	case f1HENKE:
		readF1Henke();
		break;
	case f1SASAKI:
		readF1Sasaki();
		break;
	case f1SHADOW:
		readF1Shadow();
		break;
	case f1WINDT:
		readF1Windt();
		break;
	default:
		//should never appear
		break;
	}
}

void ElementScatteringFactor::readF1BrennanCowan()
{
	std::string filename = "f1f2_BrennanCowan.dat";
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	fin.open(filename.c_str());
	fin.exceptions( std::ifstream::failbit | std::ifstream::badbit);

	//read data
	while (true)
	{
		if(fin.eof())
		{
			LOG(logERROR) << "Cannot find " << elementName << " in " << filename;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n\r");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(elementName.compare(elname) == 0)
			{
				std::istringstream is;
				Correction corr;

				//skip 3 lines
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#')
				{
					is.str(line);
					is >> corr.en >> corr.ref1 >> corr.imf1 ;
					corrections.push_back(corr);
					is.clear();
					LOG(logDEBUG) << corrections.back().en << "\t" << corrections.back().ref1 <<"\t" << corrections.back().imf1;
					getline(fin, line);
				}
				return;
			}
		}
	}

	fin.close();
}

void ElementScatteringFactor::readF1Windt()
{
	std::string filename = "f1f2_Windt.dat";
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	fin.open(filename.c_str());
	fin.exceptions( std::ifstream::failbit | std::ifstream::badbit);

	//read data
	while (true)
	{
		if(fin.eof())
		{
			LOG(logERROR) << "Cannot find " << elementName << " in " << filename;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n\r");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(elementName.compare(elname) == 0)
			{
				std::istringstream is;
				Correction corr;
				double atNumb;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();
				LOG(logDEBUG) << "Z:\t" << atNumb << "\t" << "substr:\t" <<"\t" << substrings[1];

				//skip 11 lines
				for(size_t i = 0; i < 11; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#')
				{
					is.str(line);
					is >> corr.en >> corr.ref1 >> corr.imf1 ;
					corr.ref1 -= atNumb;
					corrections.push_back(corr);
					is.clear();
					LOG(logDEBUG) << corrections.back().en << "\t" << corrections.back().ref1 <<"\t" << corrections.back().imf1;
					getline(fin, line);
				}
				return;
			}
		}
	}

	fin.close();
}

void ElementScatteringFactor::readF1Shadow()
{
	std::string filename = "f1f2_Shadow.dat";
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	fin.open(filename.c_str());
	fin.exceptions( std::ifstream::failbit | std::ifstream::badbit);

	//read data
	while (true)
	{
		if(fin.eof())
		{
			LOG(logERROR) << "Cannot find " << elementName << " in " << filename;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n\r");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(elementName.compare(elname) == 0)
			{
				std::istringstream is;
				Correction corr;
				double atNumb;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();
				LOG(logDEBUG) << "Z:\t" << atNumb << "\t" << "substr:\t" <<"\t" << substrings[1];

				//skip 8 lines
				for(size_t i = 0; i < 8; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#')
				{
					is.str(line);
					is >> corr.en >> corr.ref1 >> corr.imf1 ;
					corr.ref1 -= atNumb;
					corrections.push_back(corr);
					is.clear();
					LOG(logDEBUG) << corrections.back().en << "\t" << corrections.back().ref1 <<"\t" << corrections.back().imf1;
					getline(fin, line);
				}
				return;
			}
		}
	}

	fin.close();
}

void ElementScatteringFactor::readF1Sasaki()
{
	std::string filename = "f1f2_Sasaki.dat";
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	fin.open(filename.c_str());
	fin.exceptions( std::ifstream::failbit | std::ifstream::badbit);

	//read data
	while (true)
	{
		if(fin.eof())
		{
			LOG(logERROR) << "Cannot find " << elementName << " in " << filename;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n\r");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(elementName.compare(elname) == 0)
			{
				std::istringstream is;
				Correction corr;
				double atNumb;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();
				LOG(logDEBUG) << "Z:\t" << atNumb << "\t" << "substr:\t" <<"\t" << substrings[1];

				//skip 8 lines
				for(size_t i = 0; i < 4; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#')
				{
					is.str(line);
					is >> corr.en >> corr.ref1 >> corr.imf1 ;
					//corr.ref1 -= atNumb;
					corrections.push_back(corr);
					is.clear();
					LOG(logDEBUG) << corrections.back().en << "\t" << corrections.back().ref1 <<"\t" << corrections.back().imf1;
					getline(fin, line);
				}
				return;
			}
		}
	}

	fin.close();
}

void ElementScatteringFactor::readF1Henke()
{
	std::string filename = "f1f2_Henke.dat";
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	fin.open(filename.c_str());
	fin.exceptions( std::ifstream::failbit | std::ifstream::badbit);

	//read data
	while (true)
	{
		if(fin.eof())
		{
			LOG(logERROR) << "Cannot find " << elementName << " in " << filename;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n\r");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(elementName.compare(elname) == 0)
			{
				std::istringstream is;
				Correction corr;
				double atNumb;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();
				LOG(logDEBUG) << "Z:\t" << atNumb << "\t" << "substr:\t" <<"\t" << substrings[1];

				//skip 8 lines
				for(size_t i = 0; i < 5; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#')
				{
					is.str(line);
					is >> corr.en >> corr.ref1 >> corr.imf1 ;
					corr.ref1 -= atNumb;
					corrections.push_back(corr);
					is.clear();
					LOG(logDEBUG) << corrections.back().en << "\t" << corrections.back().ref1 <<"\t" << corrections.back().imf1;
					getline(fin, line);
				}
				return;
			}
		}
	}

	fin.close();
}
