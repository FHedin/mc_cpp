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

#ifndef FFIELD_H
#define FFIELD_H

#include <iostream>
#include <vector>

#include "Atom.h"
#include "Angle.h"
#include "Bond.h"
#include "Bond_UB.h"
#include "Dihedral.h"
#include "Dihedral_improper.h"

#include "Ensemble.h"
#include "PerConditions.h"

class FField
{
    friend std::ostream& operator<<(std::ostream& overloadStream, const FField& forf);

public:
    static const double elemchg, angstr, calory, kcaltoiu, clight;
    static const double NA, bartoiu, kboltz, rboltz, rboltzui;
    static const double mu0, chgcharmm, chgnamd, chgdlpolyiu;
    static const double sq6rt2, PI, TWOPI, SQRTPI, watercomp;
    
    FField(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
    virtual ~FField();
    
    void setNImproper(int nImproper);
    void setNDihedral(int nDihedral);
    void setNAngle(int nAngle);
    void setNUb(int nUb);
    void setNConst(int nConst);
    void setNBond(int nBond);
    void setImprList(std::vector<Dihedral_improper> imprList);
    void setDiheList(std::vector<Dihedral> diheList);
    void setAngList(std::vector<Angle> angList);
    void setUbList(std::vector<Bond_UB> ubList);
    void setBndList(std::vector<Bond> bndList);

    // The Lennard-Jones potential
    virtual double getLJ(bool dV) = 0; //all atoms
    virtual double getLJ(std::vector<Atom>& candidateVec, bool dV) = 0; //all atoms
    virtual double getLJ(Atom const& newAt, int candidate, bool dV) = 0; //one atom
    
    int getNImproper() const;
    int getNDihedral() const;
    int getNAngle() const;
    int getNUb() const;
    int getNConst() const;
    int getNBond() const;
    const std::vector<Dihedral_improper>& getImprList() const;
    const std::vector<Dihedral>& getDiheList() const;
    const std::vector<Angle>& getAngList() const;
    const std::vector<Bond_UB>& getUbList() const;
    const std::vector<Bond>& getBndList() const;

protected:
    std::vector<Atom>& at_List;
    PerConditions& pbc;
    Ensemble& ens;

    // number of bonds, angles, etc ...
    int nBond;
    int nConst;
    int nUb;
    int nAngle;
    int nDihedral;
    int nImproper;

    // data structures (vectors of objects)) storing list of bonds, angles ...
    std::vector<Bond> bndList;
    std::vector<Bond_UB> ubList;
    std::vector<Angle> angList;
    std::vector<Dihedral> diheList;
    std::vector<Dihedral_improper> imprList;
    
    // components of the energy
    double tot,pot,kin;
    double elec,vdw;
    double bond,ang,ub,dihe,impr;

    virtual void toString(std::ostream& stream) const;

};

#endif // FFIELD_H
