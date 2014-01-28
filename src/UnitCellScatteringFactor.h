/*
 * UnitCellScatteringFactor.h
 *
 *  Created on: 2 nov. 2012
 *      Author: kopp
 */

#ifndef UNITCELLSCATTERINGFACTOR_H_
#define UNITCELLSCATTERINGFACTOR_H_

#include "UnitCell.h"
#include "ElementScatteringFactor.h"
#include <map>
#include <complex>

class UnitCellScatteringFactor
{
public:
	UnitCellScatteringFactor() {}
	virtual ~UnitCellScatteringFactor() {}
	void setup(const UnitCell& cell,
			ElementScatteringFactor::SFCorrectionDataBaseEnum sfcorrdb =
					ElementScatteringFactor::f1WINDT,
			ElementScatteringFactor::SFDataBaseEnum sfdb =
					ElementScatteringFactor::f0WAASMAIER_KIRFEL);
	std::complex<double> F(double qx, double qy, double qz, double en) const;
protected:
	UnitCell unitCell;
	std::map<std::string, ElementScatteringFactor> elements;
};

#endif /* UNITCELLSCATTERINGFACTOR_H_ */
