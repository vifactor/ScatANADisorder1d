/*
 * UnitCellScatteringFactor.cpp
 *
 *  Created on: 2 nov. 2012
 *  Modified on: 22 jan. 2013
 *      Author: kopp
 */

#include "UnitCellScatteringFactor.h"

std::complex<double> UnitCellScatteringFactor::F(double qx, double qy, double qz,
		double en) const
{
	std::complex<double> res;
	std::complex<double> f;
	double Qr, k, qsq;
	std::map<std::string, ElementScatteringFactor>::const_iterator sfIt;

	//|q|^2
	qsq = qx * qx + qy * qy + qz * qz;
	//k = sin(theta) / lambda = Q / (4 pi)
	k = sqrt(qsq) / (4 * M_PI);

	res = std::complex<double>(0.0, 0.0);
	for (UnitCell::const_iterator it = unitCell.begin(); it != unitCell.end();
			++it)
	{
		sfIt = elements.find(it->_element);
		f = std::complex<double>(sfIt->second.f0(k) + sfIt->second.ref1(en),
				sfIt->second.imf1(en));

		/*scattering factor multiplied by the probability to find an atom*/
		f *= it->_occupation;
		/*scattering factor multiplied by the Debye-Waller factor*/
		f *= exp(-it->_sigmasq * qsq);

		//res = sum [f(r_i) * exp(I Q r_i), over all i]
		Qr = qx * it->_x + qy * it->_y + qz * it->_z;

		res += f * std::polar(1.0, Qr);
	}
	return res;
}

void UnitCellScatteringFactor::setup(const UnitCell& cell,
		ElementScatteringFactor::SFCorrectionDataBaseEnum sfcorrdb,
		ElementScatteringFactor::SFDataBaseEnum sfdb)
{
	unitCell = cell;

	for (UnitCell::const_iterator it = unitCell.begin(); it != unitCell.end(); ++it)
	{
		if (elements.find(it->_element) == elements.end())
		{
			elements.insert(
					std::pair<std::string, ElementScatteringFactor>(it->_element,
							ElementScatteringFactor(it->_element, sfcorrdb, sfdb)));
		}
	}
}
