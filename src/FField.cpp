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
#include <chrono>

#include "FField.hpp"

using namespace std;

FField::FField(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens,
               string _cutMode, double _ctoff, double _cuton, double _dcut)
    : at_List(_at_List), pbc(_pbc), ens(_ens), cutoff(_ctoff), cuton(_cuton), deltacut(_dcut)
{
    nBond = nConst = nUb = nAngle = nDihedral = nImproper = 0;
    tot = pot = kin = elec = vdw = bond = ang = ub = dihe = impr = 0.0;
    excl = nullptr;
    mcmvlst = nullptr;

    if(!_cutMode.compare("full"))
    {
        this->cutMode=FULL;
        cout << "Using full evaluation for non-bonded terms." << endl;
    }
    else if(!_cutMode.compare("switch"))
    {
        this->cutMode = SWITCH;
        cout << "Using switching function for non-bonded terms with [ cutoff | cuton | deltacut ] of [ " << cutoff << " | " << cuton << " | " << deltacut << " ] " << endl;
    }
//     else if(!_cutMode.compare("shift"))
//     {
//         this->cutMode = SHIFT;
//         cout << "Using shifting function for non-bonded terms." << endl;
//     }
}

FField::~FField()
{
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

const vector<Dihedral_improper>& FField::getImprList() const
{
    return imprList;
}

const vector<Dihedral>& FField::getDiheList() const
{
    return diheList;
}

const vector<Angle>& FField::getAngList() const
{
    return angList;
}

const vector<Bond_UB>& FField::getUbList() const
{
    return ubList;
}

const vector<Bond>& FField::getBndList() const
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

CUT_TYPE FField::getCutMode() const
{
    return cutMode;
}

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

void FField::setImprList(vector<Dihedral_improper>& _imprList)
{
    this->imprList = vector<Dihedral_improper>(_imprList) ;
}

void FField::setDiheList(vector<Dihedral>& _diheList)
{
  this->diheList = vector<Dihedral>(_diheList);
}

void FField::setAngList(vector<Angle>& _angList)
{
  this->angList = vector<Angle>(_angList);
}

void FField::setUbList(vector<Bond_UB>& _ubList)
{
  this->ubList = vector<Bond_UB>(_ubList);
}

void FField::setBndList(vector<Bond>& _bndList)
{
  this->bndList = vector<Bond>(_bndList);
}

void FField::setExcl(List_nonBonded& _excl)
{
    this->excl = &_excl;
}

void FField::setMcmvlst(List_Moves& _mcmvlst)
{
    this->mcmvlst = &_mcmvlst;
}

ostream& operator<<(ostream& overloadStream, const FField& forf)
{
    forf.toString(overloadStream);

    return overloadStream;
}

void FField::toString(ostream& stream) const
{
    int i;

    stream << at_List << endl;

    for ( i = 0; i < nBond; i++ )
        stream << bndList.at(i) << endl;

    for ( i = 0; i < nUb; i++ )
        stream << ubList.at(i) << endl;

    for ( i = 0; i < nAngle; i++ )
        stream << angList.at(i) << endl;

    for ( i = 0; i < nDihedral; i++ )
        stream << diheList.at(i) << endl;

    for ( i = 0; i < nImproper; i++ )
        stream << imprList.at(i) << endl;
}

// ask the FF for updating the non bonded list excl
void FField::askListUpdate(int st)
{
    if (cutMode != FULL)
    {
        cout << "Verlet list updated at step " << st;
        auto start = chrono::system_clock::now();

        excl->update_verlet_list();

        auto end = chrono::system_clock::now();
        auto elapsed_time = chrono::duration_cast<chrono::milliseconds> (end - start).count();
        cout << " : time required (milliseconds) : " << elapsed_time << endl;
    }
}
