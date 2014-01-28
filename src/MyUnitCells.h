/*
 * MyUnitCells.h
 *
 *  Created on: 25 ñ³÷. 2013
 *      Author: kopp
 */

#ifndef MYUNITCELLS_H_
#define MYUNITCELLS_H_

#include "UnitCell.h"

UnitCell getSb2Te3UnitCell(double zTe2, double zSb1,
		double pTe1 = 1.0, double pTe2 = 1.0, double pSb1 = 1.0,
		double dwTe1 = 0.0,	double dwTe2 = 0.0, double dwSb1 = 0.0);

UnitCell getSb2UnitCell(double zSb, double pSb = 1.0,	double dwSb = 0.0);

UnitCell getSb2Te2UnitCell(double zTe1, double zTe2, double zTe3,
								double zSb1z, double zSb2, double zSb3,
								double pTe1 = 1.0, double pTe2 = 1.0, double pTe3 = 1.0,
								double pSb1 = 1.0, double pSb2 = 1.0, double pSb3 = 1.0,
								double dwTe1 = 0.0,	double dwTe2 = 0.0, double dwTe3 = 0.0,
								double dwSb1 = 0.0, double dwSb2 = 0.0, double dwSb3 = 0.0);


#endif /* MYUNITCELLS_H_ */
