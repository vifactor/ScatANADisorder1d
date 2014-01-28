/*
 * Calculator.cpp
 *
 *  Created on: 31 oct. 2012
 *      Author: kopp
 */

#include "Calculator.h"

Calculator::Calculator(double lambda)
{
	double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
	double eV = GSL_CONST_MKSA_ELECTRON_VOLT;
	double A = GSL_CONST_MKSA_ANGSTROM;

	F 		= gsl_vector_complex_alloc(2);
	Fconj 	= gsl_vector_complex_alloc(2);
	W 		= gsl_vector_complex_alloc(2);

	tmp_vec = gsl_vector_complex_alloc(2);

	m_P		= gsl_matrix_complex_alloc(2, 2);
	T 		= gsl_matrix_complex_alloc(2, 2);
	m_Ps		= gsl_matrix_complex_alloc(2, 2);
	m_Exp		= gsl_matrix_complex_alloc(2, 2);

	tmp_mat = gsl_matrix_complex_alloc(2, 2);

	perm = gsl_permutation_alloc(2);

	m_L1 = 0.0;
	m_L2 = 0.0;
	m_L = 0.0;
	m_N = 0.0;

	m_lambda = lambda;
	/*calculate photon's energy*/
	m_energy = (h * c) / (eV * A) * 1.0 /lambda;

	m_sf1 = NULL;
	m_sf2 = NULL;

	m_I0 = 1.0;
	m_Ibg = 0.0;
	m_Q0 = 0.0;
}

Calculator::~Calculator()
{
	gsl_vector_complex_free(F);
	gsl_vector_complex_free(Fconj);
	gsl_vector_complex_free(W);

	gsl_vector_complex_free(tmp_vec);

	gsl_matrix_complex_free(m_P);
	gsl_matrix_complex_free(T);
	gsl_matrix_complex_free(m_Ps);
	gsl_matrix_complex_free(m_Exp);

	gsl_matrix_complex_free(tmp_mat);

	gsl_permutation_free(perm);
}

void Calculator::init(double l1, double l2,
		const UnitCellScatteringFactor * sf1,
		const UnitCellScatteringFactor * sf2, double l, double Q0)
{
	m_L1 = l1;
	m_L2 = l2;
	m_L = l;

	m_sf1 = sf1;
	m_sf2 = sf2;

	/*most intensive peak coordinate*/
	m_Q0 = Q0;
}

void Calculator::setMatrixExp(gsl_matrix_complex * m, double Q)
{
	/*set matrix Exp*/
	/*
	 *			||exp(i * q * L1)	0.0		||
	 * Exp = 	||							||
	 * 			||0.0		exp(i * q * L2 )||
	 */
	gsl_matrix_complex_set(m, 0, 0, gsl_complex_polar (1.0, Q * m_L1));
	gsl_matrix_complex_set(m, 0, 1, gsl_complex_polar (0.0, 0.0));
	gsl_matrix_complex_set(m, 1, 0, gsl_complex_polar (0.0, 0.0));
	gsl_matrix_complex_set(m, 1, 1, gsl_complex_polar (1.0, Q * m_L2));
}

void Calculator::setMatrixP(gsl_matrix_complex * m, double p, double q)
{
	/*set matrix P (Markov probabilities)*/
	gsl_matrix_complex_set(m, 0, 0, gsl_complex_rect (p, 0.0));
	gsl_matrix_complex_set(m, 0, 1, gsl_complex_rect (q, 0.0));
	gsl_matrix_complex_set(m, 1, 0, gsl_complex_rect (1.0 - p, 0.0));
	gsl_matrix_complex_set(m, 1, 1, gsl_complex_rect (1.0 - q, 0.0));

	/* DEBUG
	printf("P:\n");
	gsl_matrix_complex_fprintf (stdout, P, "%g");
	*/
}

void Calculator::setMatrixPs(gsl_matrix_complex * m, double p, double q)
{
	double S1, S2;

	/*set matrix Ps (stationary probabilities)*/
	S1 = q / (1.0 - p + q);
	S2 = (1.0 - p) / (1.0 - p + q);

	gsl_matrix_complex_set(m, 0, 0, gsl_complex_rect (S1, 0.0));
	gsl_matrix_complex_set(m, 0, 1, gsl_complex_rect (0.0, 0.0));
	gsl_matrix_complex_set(m, 1, 0, gsl_complex_rect (0.0, 0.0));
	gsl_matrix_complex_set(m, 1, 1, gsl_complex_rect (S2, 0.0));


	m_N = m_L / (m_L1 * S1 + m_L2 * S2);
	/* DEBUG
	printf("Ps:\n");
	gsl_matrix_complex_fprintf (stdout, Ps, "%g");
	*/
}


double Calculator::getIntensity(double Q)
{
	std::complex<double> cppf1, cppf2;
	gsl_complex alpha, beta;
	gsl_complex gslf1, gslf2;
	int s;
	double avFactor;

	cppf1 = m_sf1->F(0.0, 0.0, Q, m_energy);
	cppf2 = m_sf2->F(0.0, 0.0, Q, m_energy);

	gslf1 = gsl_complex_rect(cppf1.real(), cppf1.imag());
	gslf2 = gsl_complex_rect(cppf2.real(), cppf2.imag());

	avFactor = exp(-2.0 / m_N);

	alpha = gsl_complex_rect (1.0, 0.0);
	beta = gsl_complex_rect (0.0, 0.0);

	/*set vector F (scattering factors)*/
	gsl_vector_complex_set(F, 0, gslf1);
	gsl_vector_complex_set(F, 1, gslf2);

	/*set vector conj(F) (scattering factors)*/
	gsl_vector_complex_set(Fconj, 0, gsl_complex_conjugate(gslf1));
	gsl_vector_complex_set(Fconj, 1, gsl_complex_conjugate(gslf2));

	/*set exp matrix*/
	setMatrixExp(m_Exp, Q);

	/*find W = P * Exp * Ps * conj(F)  vector:*/
	/* (1) W = alpha * Ps * conj(F) + beta * W */
	gsl_blas_zgemv (CblasNoTrans, alpha, m_Ps, Fconj, beta, W);

	/*printf("W(1):\n");
	gsl_vector_complex_fprintf (stdout, W, "%g");*/

	/* (2) W = alpha * Exp * tmp_vec + beta * W */
	gsl_blas_zgemv (CblasNoTrans, alpha, m_Exp, W, beta, tmp_vec);

	/*printf("W(2):\n");
	gsl_vector_complex_fprintf (stdout, tmp_vec, "%g");*/

	/* (3) W = alpha * P * tmp_vec + beta * W */
	gsl_blas_zgemv (CblasNoTrans, alpha, m_P, tmp_vec, beta, W);

	/*Find J0 = F.(Ps * conj(F)) */
	gsl_blas_zgemv (CblasNoTrans, alpha, m_Ps, Fconj, beta, tmp_vec);
	gsl_blas_zdotu (F, tmp_vec, &J0);

	/*alpha = exp(-2 / N)*/
	alpha = gsl_complex_rect (avFactor, 0.0);
	beta = gsl_complex_rect (0.0, 0.0);

	/*find T matrix: T = alpha * P * exp + beta * T*/
	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, alpha, m_P, m_Exp, beta, T);

	/*printf("T:\n");
	gsl_matrix_complex_fprintf (stdout, T, "%g");*/

	/*Find Jns = F. (G * W)  */
	/*tmp_mat = I */
	gsl_matrix_complex_set_identity (tmp_mat);
	/*tmp_mat = I - T */
	gsl_matrix_complex_sub (tmp_mat, T);
	/*LU decomposition*/
	gsl_linalg_complex_LU_decomp(tmp_mat, perm, &s);
	/*calculate product G * W = (I - T)^(-1) W directly using LU decomposition*/
	gsl_linalg_complex_LU_solve (tmp_mat, perm, W, tmp_vec);
	/*calculate F.(G * W)*/
	gsl_blas_zdotu (F, tmp_vec, &Jns);

	/*Find Js = F.(G^2 * (I - T^N) * W)  however, this term should be negligible*/

	/*result = N *(2 * Jns + J0) - Js */
	alpha = gsl_complex_mul_real (Jns, 2.0 * avFactor);
	alpha = gsl_complex_add (alpha, J0);

	return m_N * m_I0 * GSL_REAL(alpha) * getPLGfactor(getTheta(Q)) + m_Ibg;
}


double Calculator::getPLGfactor(double theta)
{
	double sintheta = sin(theta);
	double costheta = cos(theta);
	double cos2theta = cos(2 *theta);

	return (1 + cos2theta * cos2theta)/(sintheta * sintheta * costheta);
}

double Calculator::getTheta(double Q)
{
	return asin(m_lambda * Q / (4 * M_PI));
}

void Calculator::reinit(CParameterList& params)
{
	double p, q;

	m_I0 = params["I0"];
	m_Ibg = params["Ibg"];
	p = params["p"];
	q = params["q"];

	/*set matrix P (Markov probabilities)*/
	setMatrixP(m_P, p, q);

	/*set matrix Ps (stationary probabilities)*/
	setMatrixPs(m_Ps, p, q);


	/*scale coefficient is defined by the intensity of the most intensive peak*/
	/*m_I0 = 1.0 / getIntensity(m_Q0);*/

	LOG(logINFO) << "------------reinit----------" << std::endl;
	LOG(logINFO) << "p:\t" << p << std::endl;
	LOG(logINFO) << "q:\t" << q << std::endl;
	LOG(logINFO) << "I0:\t" << m_I0 << std::endl;
	LOG(logINFO) << "Ibg:\t" << m_Ibg << std::endl;
}
