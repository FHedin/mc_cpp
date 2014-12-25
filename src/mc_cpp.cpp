/*
 *  mc_cpp : A Molecular Monte Carlo simulations software.
 *  Copyright (C) 2013  Florent Hedin
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <cstring>

#include <iostream>

#include <chrono> // for precise timing

#ifdef __unix__
#include "Global_include.hpp"
#else
#define PROGRAM_NAME argv[0]
#define VERSION_MAJOR 0
#define VERSION_MINOR 2
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Parser.hpp"

#include "Tools.hpp"

#include "PerConditions.hpp"

#include "Ensemble.hpp"
#include "Ens_NVT.hpp"

#ifdef NPT_EXPERIMENTAL
#include "Ens_NPT.hpp"
#endif

#include "Atom.hpp"

#include "FField.hpp"

#include "List_nonBonded.hpp"
#include "List_Moves.hpp"

#include "MC.hpp"
#include "MC_metropolis.hpp"

#ifdef SPAV_EXPERIMENTAL
#include "MC_spav.hpp"
#endif 

#ifdef CHARMM_EXPERIMENTAL
#include "IO_CHARMM.hpp"
#endif

#include "IO_MDBAS.hpp"

#include "FField_MDBAS.hpp"

void benchmark(FField* ff);

using namespace std;

int main(int argc, char* argv[])
{

	cout << "Welcome to " << PROGRAM_NAME << " version " << VERSION_MAJOR << '.' << VERSION_MINOR << "!!" << endl << endl;

	if (argc < 3)
	{
		cerr << "Error with arguments processing : please provide the input file name : " << endl;
		cerr << "Example : " << endl << argv[0] << " -i an_input_file.xml " << endl;
		exit(-1);
	}

#ifdef _OPENMP
	int nthreads;
	int maxthreads;
#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
		maxthreads = omp_get_max_threads();
	}
	cout << "OpenMP available ; using " << nthreads << " threads of a maximum of " << maxthreads << endl;
#endif

	//cmd line arguments parsing
	char* inpname = nullptr;
	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
			inpname = argv[++i];
		else
		{
			cerr << "Error : Argument '" << argv[i] << "' is unknown. " << endl;
			exit(-2);
		}
	}

	Parser_XML* xmlfp = nullptr;
	PerConditions* pbc = nullptr;
	Ensemble* ens = nullptr;
	vector<Atom> atomList;
	FField* ff = nullptr;
	List_nonBonded* exlst = nullptr;
	List_Moves* mvList = nullptr;
	MC* simulation = nullptr;

	// efficient xml parsing of parameters
	xmlfp = new Parser_XML(inpname, &pbc, &ens, atomList, &ff, &exlst, &mvList, &simulation, false);
	if (xmlfp->requiredBenchmark())
	{
		benchmark(ff);
	}
	delete xmlfp;

	// run simulation immediately as everything was parsed before
	simulation->run();

	/* freeing memory previously allocated with new */
	delete simulation;
	delete mvList;
	delete exlst;
	delete ff;
	delete ens;
	delete pbc;

	return EXIT_SUCCESS;
}

void benchmark(FField* ff)
{
	//compare standard and vectorized energy calls 
	//to see if it is really useful to use vectorized ones
	cout << endl << "benchmark standard energy call" << endl;
	auto start = chrono::system_clock::now();
	double estd = ff->getE();
	auto end = chrono::system_clock::now();
	auto elapsed_time_std = chrono::duration_cast<chrono::nanoseconds> (end - start).count();
	cout << "total energy is : " << estd << endl;
	cout << "time required for energy calculation was (nanoseconds) : " << elapsed_time_std << endl;
	cout << endl;

#ifdef VECTORIZED_ENER_EXPERIMENTAL
	cout << "benchmark vectorized energy call" << endl;
	start = chrono::system_clock::now();
	double evec = ff->getE(true);
	end = chrono::system_clock::now();
	auto elapsed_time_vect = chrono::duration_cast<chrono::nanoseconds> (end - start).count();
	cout << "total energy is : " << estd << endl;
	cout << "time required for energy calculation was (nanoseconds) : " << elapsed_time_vect << endl;
	cout << endl;
#endif
}
