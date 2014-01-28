/*
 * Calculator.h
 *
 *  Created on: 31 oct. 2012
 *      Author: kopp
 */

#ifndef CALCULATOR_H_
#define CALCULATOR_H_

#include "UnitCellScatteringFactor.h"
#include "NonlinearFitter.h"

#include <Log.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <complex>

class Calculator : public FitCalculator
{
public:
	Calculator(double wavelength = 1.542);
	void init(double l1, double l2,
			const UnitCellScatteringFactor * sf1,
			const UnitCellScatteringFactor * sf2, double L, double Q0);
	double getIntensity(double Q);
	double getTheta(double Q);
	double getPLGfactor(double theta);
	virtual ~Calculator();
	virtual void reinit(CParameterList& params);
	virtual double eval(CArgument& arg)
	{
		//LOG(logINFO) << getIntensity(arg[0]);
		return getIntensity(arg[0]);
	}
protected:
	void setMatrixP(gsl_matrix_complex * m, double p, double q);
	void setMatrixPs(gsl_matrix_complex * m, double p, double q);
	void setMatrixExp(gsl_matrix_complex * m, double q);

	gsl_complex Jns, Js, J0;

	gsl_vector_complex * F, * Fconj;
	gsl_vector_complex * W;

	gsl_matrix_complex * m_P;
	gsl_matrix_complex * T;
	gsl_matrix_complex * m_Ps;
	gsl_matrix_complex * m_Exp;

	gsl_matrix_complex * tmp_mat;
	gsl_vector_complex * tmp_vec;

	gsl_permutation * perm;

	double m_Q0;
	double m_I0, m_Ibg;

	double m_L1, m_L2, m_L, m_N;
	const UnitCellScatteringFactor * m_sf1, * m_sf2;

	/*x-ray parameters*/
	double m_lambda;
	double m_energy;

};

#endif /* CALCULATOR_H_ */
