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

#define _CRT_SECURE_NO_DEPRECATE

#include <cstdlib>
#include <iostream>
#include <functional>

#include "MC.hpp"
#include "Tools.hpp"
#include "FField.hpp"

using namespace std;

//MC::MC(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, List_Moves& _mvlist,
//	int _steps, int _save_freq, double _dmax_value, double _dmax_target, int _dmax_each, uint64_t _seed)
//	: at_List(_at_List), pbc(_pbc), ens(_ens), ff(_ff), mvlist(_mvlist)
//{
//	svFreq = _save_freq;
//
//	nsteps = _steps;
//	dmax = _dmax_value;
//	target = _dmax_target;
//	each = _dmax_each;
//
//	if (svFreq != 0)
//	{
//		xyz = nullptr;
//		xyz = fopen("tr.xyz", "w");
//
//		efile = nullptr;
//		efile = fopen("ener.dat", "w");
//	}
//
//	intSeed = _seed;
//	if (_seed != 0)
//		rndInit(this->intSeed);
//	else
//		rndInit();
//
//}

MC::MC(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, List_Moves& _mvlist,
	int _steps, int _save_freq, uint64_t _seed)
	: at_List(_at_List), pbc(_pbc), ens(_ens), ff(_ff), mvlist(_mvlist),
	dmax(_mvlist.getMoveLimitsList()), target(_mvlist.getTargetAcceptanceList()), each(_mvlist.getMoveUpdateFreqList())
{
	svFreq = _save_freq;

	nsteps = _steps;

	if (svFreq != 0)
	{
		xyz = nullptr;
		xyz = fopen("tr.xyz", "w");

		efile = nullptr;
		efile = fopen("ener.dat", "w");
	}

	intSeed = _seed;
	if (_seed != 0)
		rndInit(this->intSeed);
	else
		rndInit();
}

MC::~MC()
{
	if (svFreq != 0)
	{
		fclose(xyz);
		fclose(efile);
	}
}

void MC::write_traj(int st) const
{
	double crd[3] = { 0. };
	int n = ens.getN();
	const char *symb = nullptr;

	fprintf(xyz, "%d\n", n);
	fprintf(xyz, "#At step %d\n", st);

	for (int it = 0; it < n ;it++)
	{
		at_List[it].getCoords(crd);
        symb = at_List[it].getSymbol().c_str();
		fprintf(xyz, "%s\t%10.5lf\t%10.5lf\t%10.5lf\n", symb, crd[0], crd[1], crd[2]);
	}
}

void MC::rndInit()
{
	unsigned int lseed = seed();
	cout << "SEED : " << lseed << endl;
	generator.seed(lseed);
	distributionAlpha = uniform_real_distribution<double>(0.0, 1.0);
	distributionMove = uniform_real_distribution<double>(-1.0, 1.0);
}

void MC::rndInit(uint64_t _seed)
{
	cout << "SEED : " << _seed << endl;
	generator.seed(_seed);
	distributionAlpha = uniform_real_distribution<double>(0.0, 1.0);
	distributionMove = uniform_real_distribution<double>(-1.0, 1.0);
}

double MC::rndUnifMove()
{
	return distributionMove(generator);
}

double MC::rndUnifAlpha()
{
	return distributionAlpha(generator);
}

int MC::rndIntCandidate(int _n)
{
	return (int)_n * (rndUnifAlpha());
}

void MC::rndSphere(double rnd[3])
{
	double rx, ry, rz;

	do
	{
		rx = rndUnifMove();
		ry = rndUnifMove();
		rz = rndUnifMove();
	} while ((rx * rx + ry * ry + rz * rz) > 1.0);

	rnd[0] = rx;
	rnd[1] = ry;
	rnd[2] = rz;
}

void MC::scaleVec(double r[3], double dmax)
{
	r[0] *= dmax;
	r[1] *= dmax;
	r[2] *= dmax;
}

void MC::adjust_dmax(double& l_dmax, const double l_target, const int l_each, const int acc) const
{
		double ratio = (double)acc / (double)l_each;

		if (ratio > l_target / 100)
			l_dmax *= 1.10;
		else
			l_dmax *= 0.9;
}

bool MC::initial_checks_before_running()
{
	if (nsteps == 0)
	{
		cout << "Dummy MC simulation because nsteps = 0 ; returning ..." << endl;
		return false;
	}

	int nmvtyp = mvlist.getNMoveTypes();
	if (nmvtyp == 0)
	{
		cerr << "Error : Trying to perform a MC Metropolis simulation without any move selection !" << endl;
		exit(-30);
	}

	return true;
}


