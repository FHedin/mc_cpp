/*
 *  mc_cpp : A basic Monte Carlo simulations software.
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

#include "Ensemble.h"

Ensemble::Ensemble()
{
}

Ensemble::~Ensemble()
{
}

ens_type Ensemble::getType() const
{
    return type;
}

int Ensemble::getN() const
{
    return N;
}

double Ensemble::getPress() const
{
    return P;
}

double Ensemble::getVol() const
{
    return V;
}

double Ensemble::getTemp() const
{
    return T;
}

double Ensemble::getE() const
{
    return E;
}

double Ensemble::get_mu() const
{
    return mu;
}

void Ensemble::setE(double _E)
{
    E = _E;
}

void Ensemble::addE(double _dE)
{
    E += _dE;
}

void Ensemble::setP(double _P)
{
    P = _P;
}

void Ensemble::set_mu(double _mu)
{
    mu = _mu;
}

