/*
 * MyUnitCells.cpp
 *
 *  Created on: 25 ñ³÷. 2013
 *      Author: kopp
 */

#include "MyUnitCells.h"

UnitCell getSb2Te3UnitCell(double z1, double z2, double pTe1,
		double pTe2, double pSb1, double dwTe1, double dwTe2, double dwSb1)
{
	static UnitCell cell;

	cell.clear();
	/*3a Te atoms*/
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, 0.0, pTe1, dwTe1));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, 0.3333, pTe1, dwTe1));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, 0.6667, pTe1, dwTe1));
	/*6c Te atoms*/
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, z1, pTe2, dwTe2));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, 0.3333 + z1, pTe2, dwTe2));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, 0.6667 + z1, pTe2, dwTe2));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, -z1, pTe2, dwTe2));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, 0.3333 - z1, pTe2, dwTe2));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, 0.6667 - z1, pTe2, dwTe2));
	/*6c Sb atoms*/
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, z2, pSb1, dwSb1));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, 0.3333 + z2, pSb1, dwSb1));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, 0.6667 + z2, pSb1, dwSb1));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, -z2, pSb1, dwSb1));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, 0.3333 - z2, pSb1, dwSb1));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, 0.6667 - z2, pSb1, dwSb1));

	return cell;
}

UnitCell getSb2Te2UnitCell(double zTe1, double zTe2, double zTe3,
							double zSb1, double zSb2, double zSb3,
							double pTe1, double pTe2, double pTe3,
							double pSb1, double pSb2, double pSb3,
							double dwTe1,	double dwTe2, double dwTe3,
							double dwSb1, double dwSb2, double dwSb3)
{
	static UnitCell cell;

	cell.clear();

	/*Here x and y coordinates are incorrect!*/

	/*2d Te atoms*/
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, zTe1, pTe1, dwTe1));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0,	-zTe1, pTe1, dwTe1));

	/*2d Te atoms*/
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, zTe2, pTe2, dwTe2));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0,	-zTe2, pTe2, dwTe2));

	/*2c Te atoms*/
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0, zTe3, pTe3, dwTe3));
	cell.push_back(UnitCellAtom("Te", 0.0, 0.0,	-zTe3, pTe3, dwTe3));

	/*2d Sb atoms*/
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, zSb1, pSb1, dwSb1));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0,	-zSb1, pSb1, dwSb1));

	/*2d Sb atoms*/
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, zSb2, pSb2, dwSb2));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0,	-zSb2, pSb2, dwSb2));

	/*2c Sb atoms*/
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, zSb3, pSb3, dwSb3));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0,	-zSb3, pSb3, dwSb3));

	return cell;
}

UnitCell getSb2UnitCell(double zSb, double pSb, double dwSb)
{
	static UnitCell cell;

	cell.clear();

	/*Here x and y coordinates are incorrect!*/

	/*2d Sb atoms*/
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0, zSb, pSb, dwSb));
	cell.push_back(UnitCellAtom("Sb", 0.0, 0.0,	-zSb, pSb, dwSb));

	return cell;
}



