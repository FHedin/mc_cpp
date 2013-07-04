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

#include <cmath>
#include <iostream>

#include "FField.h"

const double FField::elemchg = 1.602176565e-19;
const double FField::angstr = 1.e-10;
const double FField::calory = 4.184;
const double FField::kcaltoiu = 418.4;
const double FField::clight = 299792458.;
const double FField::NA = 6.02214129e+23;
const double FField::bartoiu = 6.02214129e-3;
const double FField::kboltz = 1.3806488e-23;
const double FField::rboltz = 8.3144621;
const double FField::rboltzui = 0.83144621;
const double FField::mu0 = 1.e-7;
const double FField::chgcharmm = 332.0716;
const double FField::chgnamd = 332.0636;
const double FField::chgdlpolyiu = 138935.4835;
const double FField::sq6rt2 = 1.122462048309373;
const double FField::PI = 3.141592653589793;
const double FField::TWOPI = 6.283185307179586;
const double FField::SQRTPI = 1.772453850905516;
const double FField::watercomp = 0.007372;

FField::FField(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens,
               double _ctoff, double _dcut)
: at_List(_at_List), pbc(_pbc), ens(_ens), cutoff(_ctoff), deltacut(_dcut)
{
    nBond = nConst = nUb = nAngle = nDihedral = nImproper = 0;
    tot = pot = kin = elec = vdw = bond = ang = ub = dihe = impr = 0.0;
    excl = nullptr;
    mcmvlst = nullptr;
}

int FField::getNImproper() const
{
    return nImproper;
}

int FField::getNDihedral() const
{
    return nDihedral;
}

int FField::getNAngle() const
{
    return nAngle;
}

int FField::getNUb() const
{
    return nUb;
}

int FField::getNConst() const
{
    return nConst;
}

int FField::getNBond() const
{
    return nBond;
}

const std::vector<Dihedral_improper>& FField::getImprList() const
{
    return imprList;
}

const std::vector<Dihedral>& FField::getDiheList() const
{
    return diheList;
}

const std::vector<Angle>& FField::getAngList() const
{
    return angList;
}

const std::vector<Bond_UB>& FField::getUbList() const
{
    return ubList;
}

const std::vector<Bond>& FField::getBndList() const
{
    return bndList;
}

double FField::getDeltacut() const
{
    return deltacut;
}

double FField::getCutoff() const
{
    return cutoff;
}

FField::~FField() { }

void FField::setNImproper(int nImproper)
{
    this->nImproper = nImproper;
}

void FField::setNDihedral(int nDihedral)
{
    this->nDihedral = nDihedral;
}

void FField::setNAngle(int nAngle)
{
    this->nAngle = nAngle;
}

void FField::setNUb(int nUb)
{
    this->nUb = nUb;
}

void FField::setNConst(int nConst)
{
    this->nConst = nConst;
}

void FField::setNBond(int nBond)
{
    this->nBond = nBond;
}

void FField::setImprList(std::vector<Dihedral_improper> imprList)
{
    this->imprList = imprList;
}

void FField::setDiheList(std::vector<Dihedral> diheList)
{
    this->diheList = diheList;
}

void FField::setAngList(std::vector<Angle> angList)
{
    this->angList = angList;
}

void FField::setUbList(std::vector<Bond_UB> ubList)
{
    this->ubList = ubList;
}

void FField::setBndList(std::vector<Bond> bndList)
{
    this->bndList = bndList;
}

void FField::setExcl(List_nonBonded& _excl)
{
    this->excl = &_excl;
}

void FField::setMcmvlst(List_Moves& _mcmvlst)
{
    this->mcmvlst = &_mcmvlst;
}

std::ostream& operator<<(std::ostream& overloadStream, const FField& forf)
{
    forf.toString(overloadStream);

    return overloadStream;
}

void FField::toString(std::ostream& stream) const
{
    int i;

    for ( i = 0; i < ens.getN(); i++ )
        stream << at_List.at(i) << std::endl;

    for ( i = 0; i < nBond; i++ )
        stream << bndList.at(i) << std::endl;

    for ( i = 0; i < nUb; i++ )
        stream << ubList.at(i) << std::endl;

    for ( i = 0; i < nAngle; i++ )
        stream << angList.at(i) << std::endl;

    for ( i = 0; i < nDihedral; i++ )
        stream << diheList.at(i) << std::endl;

    for ( i = 0; i < nImproper; i++ )
        stream << imprList.at(i) << std::endl;
}

