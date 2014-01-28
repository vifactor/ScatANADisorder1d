/*
 * Engine.h
 *
 *  Created on: 31 oct. 2012
 *      Author: kopp
 */

#ifndef ENGINE_H_
#define ENGINE_H_

#include "DataReader.h"
#include "Calculator.h"
#include <libconfig.h++>

class Engine
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "Engine" + m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};
	Engine();
	virtual ~Engine();
	void read(std::string filename);
	void init();
	void exec();
	void save();
protected:
	struct Sample
	{
		double l1, l2;
		double p, q;
		double c1, c2, c;
		double thickness;
		std::string structfile1, structfile2;
	} sample;
	struct Settings
	{
		double lambda;
		double background;
		double scale;
		double q0;
		std::size_t nbFitSteps;
		std::string datafile;
		std::string outfile;
		double q_min, q_max;
		size_t q_samp;
		ElementScatteringFactor::SFDataBaseEnum f0db;
		ElementScatteringFactor::SFCorrectionDataBaseEnum f1db;
	} settings;
	UnitCell uc1, uc2;
	UnitCellScatteringFactor sf1, sf2;

	CParameterList m_cParameters;
	FParameterList m_fParameters;
	FResidualList m_fResiduals;
	CArgumentList m_cArguments;
	Calculator m_Calculator;
	Fitter * m_Fitter;

	ElementScatteringFactor::SFDataBaseEnum defineF0db(const std::string& str) const;
	ElementScatteringFactor::SFCorrectionDataBaseEnum defineF1db(const std::string& str) const;
	void readSample(const libconfig::Setting& smpl);
	void readSettings(const libconfig::Setting& stgs);
};

#endif /* ENGINE_H_ */
