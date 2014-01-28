/*
 * ElementScatteringFactor.h
 *
 *  Created on: 1 nov. 2012
 *      Author: kopp
 */

#ifndef ELEMENTSCATTERINGFACTOR_H_
#define ELEMENTSCATTERINGFACTOR_H_

#include <Log.h>
#include <StringTools.h>

#include <gsl/gsl_math.h>

#include <string>
#include <fstream>
#include <algorithm>

class ElementScatteringFactor
{
public:
	enum SFDataBaseEnum {f0WAASMAIER_KIRFEL, f0UNKNOWN};
	enum SFCorrectionDataBaseEnum {f1BRENNAN_COWAN, f1HENKE, f1SASAKI, f1SHADOW, f1WINDT, f1UNKNOWN};
	ElementScatteringFactor(std::string elname, SFCorrectionDataBaseEnum dbf1 = f1BRENNAN_COWAN, SFDataBaseEnum dbf0 = f0WAASMAIER_KIRFEL);
	double f0(double k) const;
	double ref1(double en) const;
	double imf1(double en) const;
	virtual ~ElementScatteringFactor();
protected:
	enum {NB_COEFS = 5};
	/*
	 * f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]
	 *              i=1,5
	 *k = sin(theta) / lambda
	 */
	double a[NB_COEFS];
	double b[NB_COEFS];
	double c;

	struct Correction
	{
		double en;
		double ref1;
		double imf1;
	};
	struct ComparisonClass
	{
	  bool operator() (const Correction& c1, const Correction& c2)
	  {
		  return (c1.en < c2.en);
	  }
	} cmpObj;
	std::vector<Correction> corrections;

	std::string elementName;

	void readF0DB(SFDataBaseEnum db);
	void readF1DB(SFCorrectionDataBaseEnum db);
	void readF0WaasmaierKirfel();
	void readF1BrennanCowan();
	void readF1Sasaki();
	void readF1Henke();
	void readF1Shadow();
	void readF1Windt();

};

#endif /* ELEMENTSCATTERINGFACTOR_H_ */
