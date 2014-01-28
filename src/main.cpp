//============================================================================
// Name        : Disorder1dAnalyticalv5.cpp
// Author      : Viktor Kopp
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Engine.h"

using namespace std;

int main(int argc, char *argv[])
{
	OutputTO::Stream()=stdout;
	Engine engine;
	try
	{
		if(argc > 1)
		{
			for(int argi = 1; argi < argc; ++argi)
			{
				engine.read(argv[argi]);
				engine.init();
				engine.exec();
			}
		}
		else
		{
			engine.read("default.cfg");
			engine.init();
			engine.exec();
		}
		LOG(logINFO) << "Done.";
	}catch(const std::exception& ex)
	{
		LOG(logERROR) << ex.what();
	}

	return 0;
}

