/*
 * NonlinearFitter.cpp
 *
 *  Created on: 13 july 2012
 *      Author: kopp
 */

#include "NonlinearFitter.h"

//using namespace NonlinearFit;

int n2fbNormFunc(int * n, int * m, double * x, int * nf, double * v,
		int * ui, double * ur, int (* uf)())
{
	register int i;
	static double cresid;

	/* ui is a pointer to the fitter */
	Fitter * fitter = reinterpret_cast<Fitter *>(ui);

	/* reinitialize calculator */
	fitter->resetCalculator(x);

	/* calculate residuals */

	for (i = 0; i < *n; ++i)
	{
		//for log-fit
		cresid = fitter->calculateResidual(i);
		/*if(cresid > 0 )
		{
			v[i] = log(cresid) - fitter->residualsList[i];
		}
		else
		{
			v[i] = 0.0 - fitter->residualsList[i];
		}*/
		//for lin-fit
		v[i] = cresid - fitter->residualsList[i];
	}

	return 0;
}

void levmarNormFunc(double * x, double * v, int m, int n, void * adata)
{
	register int i;
	register double resid = 0;

	/* adata is a pointer to the fitter */
	Fitter * fitter = static_cast<Fitter *>(adata);

	/* reinitialize calculator */
	fitter->resetCalculator(x);

	/* calculate residuals */
	for (i = 0; i < n; ++i)
	{
		resid = fitter->calculateResidual(i);
		if(resid > 0)
		{
			v[i] = log(resid);
		}
		else
		{
			v[i] = 0.0;
		}
	}
}

Fitter::Fitter(Fitter::FitType t, Fitter::LogMode m)
{
	ftype = t;
	lmode = m;
	calculator = NULL;
}

Fitter::~Fitter()
{

}

void Fitter::resetCalculator(const double * x)
{
	resetCParameters(x);
	calculator->reinit(calcParameterList);
}

void Fitter::resetCParameters(const double * x)
{
	for (size_t i = 0; i < fitParameterList.size(); ++i)
	{
		calcParameterList[fitParameterList[i].name] = x[i];
		LOG(logDEBUG) << i << "\t" <<fitParameterList[i].name << "\t" << x[i];
	}
}

void Fitter::resetFParameters(const double * x)
{
	for (size_t i = 0; i < fitParameterList.size(); ++i)
	{
		fitParameterList[i].value = x[i];
	}
}

double Fitter::calculateResidual(int i)
{
	LOG(logDEBUG1)<< i << "\t" << argumentsList[i][0]
	                  << "\t" << calculator->eval(argumentsList.at(i))
	                  << "\t" << residualsList[i];
	return calculator->eval(argumentsList.at(i));
}

void Fitter::printLevmarInfo(const double * info)
{
	LOG(logINFO) << "Reason to stop iterations:";
	switch(static_cast<int>(info[6]))
	{
	case 1:
		LOG(logINFO)<<"Small gradient J^T f";
		break;
	case 2:
		LOG(logINFO)<<"Small Dp";
		break;
	case 3:
		LOG(logINFO)<<"Max nb iterations";
		break;
	case 4:
		LOG(logINFO)<<"Singular matrix";
		break;
	case 5:
		LOG(logINFO)<<"No further error reduction is possible";
		break;
	case 6:
		LOG(logINFO)<<"Small ||f||^2";
		break;
	default:
		LOG(logWARNING)<<"Invalid parameters";
		return;
		break;
	}
	LOG(logINFO)<<"In "<< info[5] <<" iterations ||f||^2 reduced from "<<sqrt(info[0])<<" to "<<sqrt(info[1]);
	LOG(logINFO)<<"Number of function evaluations:"<<info[7];
	LOG(logINFO)<<"Number of Jacobian evaluations:"<<info[8];
}

void Fitter::printNl2solInfo(const int * ninfo, const double * info)
{
	enum {ivNITER = 30,
		ivNFCALL = 5,
		ivNGCALL = 29,
		ivNGCOV = 52};
	enum {vF = 9, vF0 = 12};

	LOG(logINFO) << "Reason to stop iterations:";
	switch(ninfo[0])
	{
	case 3:
		LOG(logINFO)<<"X-convergence";
		break;
	case 4:
		LOG(logINFO)<<"Relative function convergence";
		break;
	case 5:
		LOG(logINFO)<<"Both X- and relative function convergence";
		break;
	case 6:
		LOG(logINFO)<<"Absolute function convergence";
		break;
	case 10:
		LOG(logINFO)<<"Max nb iterations";
		break;
	default:
		LOG(logWARNING)<<"Invalid parameters:\t" << ninfo[0];
		break;
	}
	LOG(logINFO)<<"In " << ninfo[ivNITER] <<" iterations ||f||^2 reduced from "<<sqrt(2 * info[vF0])<<" to "<<sqrt(2 * info[vF]);
	LOG(logINFO)<<"Number of function evaluations:\t"<<ninfo[ivNFCALL] + ninfo[ivNGCALL];
	LOG(logINFO)<<"Number of Jacobian evaluations:\t"<<ninfo[ivNGCALL] - ninfo[ivNGCOV];
}

void Fitter::dofit(int nbit)
{
	/* if number of fit parameters is not equal to zero
	 *  and
	 * max number of iterations is more than zero
	 *
	 * perform fit
	 */
	if((fitParameterList.size() > 0) && (nbit > 0))
	{
		switch(ftype)
		{
		case ftNL2SOL:
			performNl2sol(nbit);
			break;
		case ftLEVMAR:
		default:
			performLevmar(nbit);
			break;
		}
	}
	printState();
}

void Fitter::performLevmar(int nbit)
{
	int nbParams = fitParameterList.size();
	int nbResids = argumentsList.size();
	/* levmar specific parameters */
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	std::unique_ptr<double[]> x_ptr(new double[nbParams]);
	std::unique_ptr<double[]> xlb_ptr(new double[nbParams]);
	std::unique_ptr<double[]> xub_ptr(new double[nbParams]);
	std::unique_ptr<double[]> xdscl_ptr(new double[nbParams]);
	std::unique_ptr<double[]> xcovar_ptr(new double[nbParams * nbParams]);

	//initialize levmar calc options
	opts[0] = LM_INIT_MU;
	opts[1] = 1E-15;
	opts[2] = 1E-15;
	opts[3] = 1E-1;
	opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

	//initialize levmar parameters
	for(size_t i = 0; i < fitParameterList.size(); ++i)
	{
		x_ptr[i]	= fitParameterList[i].value;
		xlb_ptr[i]	= fitParameterList[i].lbvalue;
		xub_ptr[i]	= fitParameterList[i].ubvalue;
		xdscl_ptr[i]= fitParameterList[i].scvalue;
	}

	dlevmar_bc_dif(levmarNormFunc, x_ptr.get(), residualsList.get(), nbParams,
			nbResids, xlb_ptr.get(), xub_ptr.get(), xdscl_ptr.get(), nbit,
			opts, info, NULL, xcovar_ptr.get(), reinterpret_cast<void *>(this));

	resetFParameters(x_ptr.get());
	printLevmarInfo(info);
}

void Fitter::performNl2sol(int nbit)
{
	enum {vF0 = 12};
	//the value of f(xinit) for final output
	double fxinit = 0;


	int nbParams = fitParameterList.size();
	int nbResids = argumentsList.size();

	// RN2GB: "LIV...... LENGTH OF IV... LIV MUST BE AT LEAST 4*P + 82"
	int liv	= 4 * nbParams + 82;
	// RN2GB: "LV....... LENGTH OF V...  LV  MUST BE AT LEAST 105 + P*(N + 2*P + 17) + 2*N"
	int lv = 105 + nbParams * (nbResids + 2 * nbParams + 17) + 2 * nbResids;
	std::unique_ptr<double[]> x_ptr(new double[nbParams]);
	std::unique_ptr<double[]> xb_ptr(new double[2 * nbParams]);
	std::unique_ptr<int[]> iv_ptr(new int[liv]);
	std::unique_ptr<double[]> v_ptr(new double[lv]);

	//initialize nl2sol parameters
	for (size_t i = 0; i < fitParameterList.size(); ++i)
	{
		x_ptr[i] = fitParameterList[i].value;
		xb_ptr[2 * i] = fitParameterList[i].lbvalue;
		xb_ptr[2 * i + 1] = fitParameterList[i].ubvalue;
	}
	//initialize default values for iv and v optimization settings
	int kind = 1;
	iv_ptr[0] = 0;
	divset_(&kind, iv_ptr.get(), &liv, &lv, v_ptr.get());
	//turn of output
	iv_ptr[20] = 0;
	//max iterations allowed
	iv_ptr[17] = 0;

	//try to execute optimization algorithm
	dn2fb_(&nbResids, &nbParams, x_ptr.get(), xb_ptr.get(), n2fbNormFunc,
			iv_ptr.get(), &liv, &lv, v_ptr.get(),
			reinterpret_cast<int *>(this), NULL, NULL);

	//storage size for iv_ptr or v_ptr was too small
	if((iv_ptr[0] == 15) || (iv_ptr[0] == 16))
	{
		//reallocate iv_ptr and v_ptr
		liv = iv_ptr[43];
		lv = iv_ptr[44];
		iv_ptr.reset(new int[liv]);
		v_ptr.reset(new double[lv]);
		//reset defaults
		divset_(&kind, iv_ptr.get(), &liv, &lv, v_ptr.get());
		//turn of output
		iv_ptr[20] = 0;
		//max iterations allowed
		iv_ptr[17] = 0;

		//execute optimization algorithm once to get f(xinit)
		dn2fb_(&nbResids, &nbParams, x_ptr.get(), xb_ptr.get(), n2fbNormFunc,
				iv_ptr.get(), &liv, &lv, v_ptr.get(),
				reinterpret_cast<int *>(this), NULL, NULL);
		//v(vF0 = 12) is the function value of f (x) at the start of the last iteration
		fxinit = v_ptr[vF0];

		//max iterations allowed
		iv_ptr[17] = nbit;
		//execute optimization algorithm nbit-1 times to optimize
		dn2fb_(&nbResids, &nbParams, x_ptr.get(), xb_ptr.get(), n2fbNormFunc,
				iv_ptr.get(), &liv, &lv, v_ptr.get(),
				reinterpret_cast<int *>(this), NULL, NULL);
	}

	resetFParameters(x_ptr.get());
	//reset v(vF0) to fxinit
	v_ptr[vF0] = fxinit;
	printNl2solInfo(iv_ptr.get(), v_ptr.get());
}

void Fitter::init(FitCalculator * c, const CArgumentList& alst,
		const FResidualList& rlst, const FParameterList& fplst,
		const CParameterList& cplst)
{
	//init calculator
	calculator = c;

	//init arguments and residuals
	CArgument arg;
	residualsList.reset(new double[rlst.size()]);
	if(!argumentsList.empty())
	{
		argumentsList.clear();
	}
	for (size_t i = 0; i < rlst.size(); i++)
	{
		//for log-fit
		/*if(rlst.at(i) > 0)
		{
			residualsList[i] = log(rlst.at(i));
		}
		else
		{
			residualsList[i] = 0.0;
		}*/
		//for lin-fit
		residualsList[i] = rlst.at(i);

		arg = alst.at(i);
		argumentsList.push_back(arg);
	}
	//init calc and fit parameters
	fitParameterList = fplst;
	calcParameterList = cplst;

	printState();
}

const FParameterList& Fitter::getFParameters()
{
	return fitParameterList;
}

void Fitter::printState()
{
	LOG(logDEBUG1) << "nb of fit paramaters:\t" << fitParameterList.size();
	/*if number of fit parameters is not equal to zero
	print all fit parameters*/
	if(fitParameterList.size() > 0)
	{
		LOG(logDEBUG1) << "They are:";
		for (size_t i = 0; i < fitParameterList.size(); i++)
		{
			LOG(logDEBUG1) << i << ":\t" << fitParameterList[i].name << "\t"
					<< fitParameterList[i].value << "\t["
					<< fitParameterList[i].lbvalue << ":"
					<< fitParameterList[i].ubvalue << "] / "
					<< fitParameterList[i].scvalue;
		}
	}

	LOG(logDEBUG1) << "nb of calculator parameters:\t" << calcParameterList.size();
	/* list of all calculator parameters*/
	CParameterList::iterator it;
	for(it = calcParameterList.begin(); it != calcParameterList.end(); ++it)
	{
		LOG(logDEBUG1) << "\t" << it->first << "\t" << it->second;
	}

	LOG(logDEBUG1) << "Nb of residuals:\t" << argumentsList.size();

}
