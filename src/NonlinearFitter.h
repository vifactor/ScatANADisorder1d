/*
 * NonlinearFitter.h
 *
 *  Created on: 12 july 2012
 *      Author: Viktor Kopp
 *
 *      Provides a convenient interface to various nonlinear fit routines
 *      with simple boundary conditions
 *
 *      Version 0.3 from 120723
 */

#ifndef NONLINEARFITTER_H_
#define NONLINEARFITTER_H_

#include <map>
#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include <Log.h>
#include <levmar.h>
#include <port.h>

/*namespace NonlinearFit
{*/

typedef std::string ParameterNameType;
typedef double ParameterValueType;
typedef std::map<ParameterNameType, ParameterValueType> CParameterList;
typedef std::vector<double> CArgument;
typedef std::vector<CArgument> CArgumentList;

struct FParameterStruct
{
	ParameterNameType name;
	ParameterValueType value;
	ParameterValueType lbvalue;
	ParameterValueType ubvalue;
	ParameterValueType scvalue;
};

typedef std::vector<FParameterStruct> FParameterList;
typedef std::vector<double> FResidualList;

class FitCalculator
{
public:
	FitCalculator() {};
	virtual ~FitCalculator() {};
	virtual void reinit(CParameterList& params) = 0;
	virtual double eval(CArgument& arg) = 0;
protected:
};

class Fitter
{
public:
	enum FitType {ftLEVMAR, ftNL2SOL};
	enum LogMode {logSILENT, logFULL};
	Fitter(FitType t = ftLEVMAR, LogMode m = logFULL);
	void init(FitCalculator * c, const CArgumentList& alst,
			const FResidualList& rlst, const FParameterList& fplst,
			const CParameterList& cplst);
	const FParameterList& getFParameters();
	void dofit(int nbit);
	virtual ~Fitter();
protected:
	FitType ftype;
	LogMode lmode;
	/*calculates resuduals for levmar fitting proceedure; passed to dlevmar_*() */
	friend void levmarNormFunc(double *x, double *v, int m, int n, void * adata);
	/*calculates resuduals for n2fb fitting proceedure; passed to n2fb() */
	friend int n2fbNormFunc(int * m, int * n, double * x, int * nf, double * v,
			int * ui, double * ur, int (* uf)());

	void performLevmar(int nbit);
	void performNl2sol(int nbit);

	/*calculates the function value corresponding to i-th argument*/
	double calculateResidual(int i);
	/*reinitializes the calculator*/
	void resetCalculator(const double * x);
	/*resets CParameters according to x-array */
	void resetCParameters(const double * x);
	/* resets FParameters according to x-array */
	void resetFParameters(const double * x);

	/*map of calculator parameters*/
	CParameterList calcParameterList;
	/*vector of fit parameters*/
	FParameterList fitParameterList;
	/*vector of calculator arguments*/
	CArgumentList argumentsList;
	/*vector residuals*/
	std::unique_ptr<double[]> residualsList;

	/*nl2sol specific parameters*/
	friend int uf();	//dummy function

	FitCalculator * calculator;

	//information functions
	void printLevmarInfo(const double * info);
	void printNl2solInfo(const int * ninfo, const double * info);
	void printState();
};

/*}*/
#endif /* NONLINEARFITTER_H_ */
