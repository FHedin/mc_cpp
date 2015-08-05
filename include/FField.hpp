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

#ifndef FFIELD_H
#define FFIELD_H

class List_nonBonded;
class List_Moves;

#include <iostream>
#include <vector>
#include <string>

//#include "Global_include.hpp"

#include "Ensemble.hpp"
#include "PerConditions.hpp"

#include "Angle.hpp"
#include "Bond.hpp"
#include "Bond_UB.hpp"
#include "Dihedral.hpp"
#include "Dihedral_improper.hpp"

#include "AtomList.hpp"

#include "List_nonBonded.hpp"
#include "List_Moves.hpp"

enum CUT_TYPE
{
    FULL=0,
    SWITCH=1,
    SHIFT=2
};

class FField
{
    friend std::ostream& operator<<( std::ostream& overloadStream, const FField& forf );

public:

    FField(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens,
           std::string _cutMode="switch", double _ctoff=12.0, double _cuton=10.0, double _dcut=2.0);

    virtual ~FField();

    // setters and getters
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
    void setExcl(List_nonBonded& excl);
    void setMcmvlst(List_Moves& mcmvlst);

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

    double getDeltacut() const;
    double getCutoff() const;
    CUT_TYPE getCutMode() const;

    virtual double getE() = 0;
    virtual void getE(double ener[10]) = 0;
    
    virtual void askListUpdate(int st);

protected:

    AtomList& at_List;
    PerConditions& pbc;
    Ensemble& ens;

    //exclude list info stored in an object of type List_nonBonded
    List_nonBonded* excl;
    List_Moves* mcmvlst;

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
    double tot, pot, kin;
    double elec, vdw;
    double bond, ang, ub, dihe, impr;

    // cutoff, switching, shifting, etc ...
    CUT_TYPE cutMode;
    double cutoff;
    double cuton;
    double deltacut;

    virtual double getEtot() = 0;
    virtual double getEswitch() = 0;

    virtual void computeNonBonded_full() = 0;
    virtual void computeNonBonded14() = 0;

    virtual void computeNonBonded_switch() = 0;
    virtual void computeNonBonded14_switch() = 0;

    virtual void computeEbond() = 0;
    virtual void computeEang() = 0;
    virtual void computeEub() = 0;
    virtual void computeEdihe() = 0;
    virtual void computeEimpr() = 0;

    virtual void toString(std::ostream& stream) const;

};

#endif // FFIELD_H
