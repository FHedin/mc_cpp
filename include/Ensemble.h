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

#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "Global_include.h"

#include <string>

enum ens_type
{
    NVT = 1/*, NPT = 2, NVE = 3, muVT = 4*/
};

class Ensemble
{
public:
    Ensemble();
    virtual ~Ensemble();

    virtual std::string whoami() = 0;

    ens_type getType() const;
    int getN() const;
    double getPress() const;
    double getVol() const;
    double getVolT() const;
    double getTemp() const;
    //    double getE() const;
    double get_mu() const;

    //    void setE(double _E);
    //    void addE(double _dE);
    void setP(double _P);
    void set_mu(double _mu);

protected:
    ens_type type;
    int N;
    double P;
    double V;
    double T;
    //    double E;
    double mu;

};

#endif // ENSEMBLE_H
