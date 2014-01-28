/*
 * Engine.cpp
 *
 *  Created on: 31 oct. 2012
 *      Author: kopp
 */

#include "Engine.h"

Engine::Engine()
{
	m_Fitter = new Fitter(Fitter::ftNL2SOL);
}

Engine::~Engine()
{
	delete m_Fitter;
}

ElementScatteringFactor::SFDataBaseEnum Engine::defineF0db(const std::string& str) const
{
	if (str.compare("WAAS_KIRF") == 0)
	{
		LOG(logINFO) << "\tWaasmaier-Kirfel database for f0." << std::endl;
		return ElementScatteringFactor::f0WAASMAIER_KIRFEL;
	}
	return ElementScatteringFactor::f0UNKNOWN;
}

ElementScatteringFactor::SFCorrectionDataBaseEnum Engine::defineF1db(const std::string& str) const
{
	if (str.compare("HENKE") == 0)
	{
		LOG(logINFO) << "\tHenke database for f1." << std::endl;
		return ElementScatteringFactor::f1HENKE;
	}
	if (str.compare("SASAKI") == 0)
	{
		LOG(logINFO) << "\tSasaki database for f1." << std::endl;
		return ElementScatteringFactor::f1SASAKI;
	}
	if (str.compare("SHADOW") == 0)
	{
		LOG(logINFO) << "\tShadow database for f1." << std::endl;
		return ElementScatteringFactor::f1SHADOW;
	}
	if (str.compare("WINDT") == 0)
	{
		LOG(logINFO) << "\tWindt database for f1." << std::endl;
		return ElementScatteringFactor::f1WINDT;
	}
	if (str.compare("BRENNAN_COWAN") == 0)
	{
		LOG(logINFO) << "\tBrennan-Cowan database for f1."<< std::endl;
		return ElementScatteringFactor::f1BRENNAN_COWAN;
	}
	return ElementScatteringFactor::f1UNKNOWN;
}

void Engine::read(std::string filename)
{
	libconfig::Config cfg;
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(filename.c_str());

		const libconfig::Setting& root = cfg.getRoot();
		const libconfig::Setting& smpl = root["Sample"];
		const libconfig::Setting& stgs = root["Settings"];

		readSample(smpl);
		readSettings(stgs);

	} catch (const libconfig::FileIOException &fioex)
	{
		throw Engine::Exception("I/O error while reading file:\t" + filename);
	} catch (const libconfig::ParseException &pex)
	{
		throw Engine::Exception("Parse error at " +
									toString(pex.getFile()) + ":" +
									toString(pex.getLine()) + " - " +
									toString(pex.getError()));
	}catch(const libconfig::SettingNotFoundException &nfex)
	{
		throw Engine::Exception(toString(nfex.getPath()));
	}catch(libconfig::SettingTypeException& tex)
	{
		throw Engine::Exception(toString(tex.getPath()) + "(" + toString(tex.what()) + ")");
	}catch(Engine::Exception &ex)
	{
		throw ex;
	}
}

void Engine::readSample(const libconfig::Setting& smpl)
{
	sample.thickness = smpl["thickness"];
	LOG(logINFO) << "Sample thickness:\t" << sample.thickness << std::endl;

	sample.c = smpl["c"];
	sample.c1 = smpl["c1"];
	sample.c2 = smpl["c2"];
	LOG(logINFO) << "Concentrations (c, c1, c2):\t" << sample.c << "\t"
			<< sample.c1 << "\t" << sample.c2 << std::endl;

	sample.p = smpl["p"];
	sample.q = smpl["q"];
	//sample.q = (sample.c - sample.c2) / (sample.c - sample.c1) * (sample.p - 1);
	LOG(logINFO) << "Probabilities (p, q):\t" << sample.p << "\t" << sample.q << std::endl;

	sample.structfile1 = smpl["f1"].c_str();
	sample.structfile2 = smpl["f2"].c_str();
	LOG(logINFO) << "Structure files (f1, f2):\t" << sample.structfile1 << "\t" << sample.structfile2 << std::endl;

	sample.l1 = smpl["l1"];
	sample.l2 = smpl["l2"];
	LOG(logINFO) << "Lengths of units (l1, l2):\t" << sample.l1 << "\t" << sample.l2 << std::endl;
}

void Engine::readSettings(const libconfig::Setting& stgs)
{
	settings.lambda = stgs["lambda"];
	LOG(logINFO) << "Wavelength:\t" << settings.lambda << std::endl;

	settings.background = stgs["background"];
	LOG(logINFO) << "Background level:\t" << settings.background << std::endl;

	settings.scale = stgs["scale"];
	LOG(logINFO) << "Scale coefficient:\t" << settings.scale << std::endl;

	settings.q0 = stgs["q0"];
	LOG(logINFO) << "Most intensive peak:\t" << settings.q0 << std::endl;

	settings.q_max = stgs["scan"]["q_max"];
	settings.q_min = stgs["scan"]["q_min"];
	settings.q_samp = stgs["scan"]["q_samp"];
	LOG(logINFO) << "Scan settings:\t[" << settings.q_min << " : " << settings.q_max <<"] / " <<settings.q_samp << std::endl;

	settings.datafile = stgs["datafile"].c_str();
	LOG(logINFO) << "Data file:\t" << settings.datafile << std::endl;

	settings.outfile = stgs["outfile"].c_str();
	LOG(logINFO) << "Output file:\t" << settings.outfile << std::endl;

	settings.nbFitSteps = stgs["nbIterations"];
	LOG(logINFO) << "Nb of fitter iterations:\t" << settings.nbFitSteps << std::endl;

	settings.f0db = defineF0db(stgs["f0_db"]);
	settings.f1db = defineF1db(stgs["f1_db"]);
}

void Engine::init()
{
	UnitCellReader cellReader;
	FParameterStruct fParameter;
	DataReader dataReader;
	std::vector<double> arguments;
	CArgument cArgument;

	uc1 = cellReader.read(sample.structfile1);
	uc2 = cellReader.read(sample.structfile2);

	/*initialize scattering factors*/
	sf1.setup(uc1, settings.f1db, settings.f0db);
	sf2.setup(uc2, settings.f1db, settings.f0db);

	/*initialize calculator*/
	m_Calculator.init(sample.l1, sample.l2, &sf1, &sf2, sample.thickness, settings.q0);

	/*setup fit parameters*/

	/*intensity scale coefficient*/
	fParameter.name = "I0";
	fParameter.value = settings.scale;
	fParameter.scvalue = 1.0;
	fParameter.lbvalue = 1.0e-8;
	fParameter.ubvalue = 1.0e10;
	m_fParameters.push_back(fParameter);
	m_cParameters[fParameter.name] = fParameter.value;

	/*intensity background*/
	fParameter.name = "Ibg";
	fParameter.value = settings.background;
	fParameter.scvalue = 10.0;
	fParameter.lbvalue = 0.0;
	fParameter.ubvalue = 1.0e10;
	//m_fParameters.push_back(fParameter);
	m_cParameters[fParameter.name] = fParameter.value;

	/*p-parameter*/
	fParameter.name = "p";
	fParameter.value = sample.p;
	fParameter.scvalue = 0.1;
	fParameter.lbvalue = 0.0;
	fParameter.ubvalue = 1.0;
	m_fParameters.push_back(fParameter);
	m_cParameters[fParameter.name] = fParameter.value;

	/*q-parameter*/
	fParameter.name = "q";
	fParameter.value = sample.q;
	fParameter.scvalue = 0.1;
	fParameter.lbvalue = 0.0;
	fParameter.ubvalue = 1.0;
	m_fParameters.push_back(fParameter);
	m_cParameters[fParameter.name] = fParameter.value;

	m_Calculator.reinit(m_cParameters);

	dataReader.readFile(settings.datafile);
	dataReader.getColumn(m_fResiduals, "[intensity]");
	dataReader.getColumn(arguments, "[q]");

	for(std::size_t i = 0; i < arguments.size(); ++i)
	{
		cArgument.push_back(arguments[i]);
		m_cArguments.push_back(cArgument);
	}
	m_Fitter->init(&m_Calculator, m_cArguments, m_fResiduals, m_fParameters, m_cParameters);
}

void Engine::exec()
{
	double Q, qstep, intens;
	std::ofstream fout;
	//size_t n = 30;
	//double resid;

	m_Fitter->dofit(settings.nbFitSteps);

	fout.open(settings.outfile.c_str());
	qstep = (settings.q_max - settings.q_min) / (settings.q_samp-1);
	for(size_t i = 0; i < settings.q_samp; ++i)
	{
		Q = settings.q_min + i * qstep;

		intens = m_Calculator.getIntensity(Q);

		fout << Q << "\t" << intens << std::endl;
	}
	fout.close();

	/*output fitted values*/
	for(std::size_t i = 0; i < m_Fitter->getFParameters().size(); ++i)
	{
		/*merge cparams and fparams*/
		m_cParameters[m_Fitter->getFParameters()[i].name] = m_Fitter->getFParameters()[i].value;
		/*output fitted values*/
		std::cout << m_Fitter->getFParameters()[i].name << "\t" << m_Fitter->getFParameters()[i].value << std::endl;
	}

	/*output residuals*/
	/*LOG(logINFO) << "Calculating residuals...";
	fout.open("residuals.txt");
	for(double p = 0.0; p <= 1.0; p += 1.0 / n)
	{
		for(double q = 0.0; q <= 1.0; q += 1.0 / n)
		{
			m_cParameters["p"] = p;
			m_cParameters["q"] = q;
			m_Calculator.reinit(m_cParameters);

			resid = 0.0;
			for(size_t i = 0; i < m_cArguments.size(); ++i)
			{
				resid += gsl_pow_2(m_Calculator.eval(m_cArguments[i]) - m_fResiduals[i]);
			}
			fout << p << "\t" << q << "\t" << sqrt(resid) << std::endl;
		}
	}
	fout.close();*/
}
