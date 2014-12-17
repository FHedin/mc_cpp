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

#ifndef FFIELD_MDBAS_H
#define	FFIELD_MDBAS_H

#include "Global_include.hpp"

#include "FField.hpp"

class FField_MDBAS : public FField
{
public:
    FField_MDBAS(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens,
                 std::string _cutMode="switch", double _ctoff=12.0, double _cton=10.0, double _dcut=2.0);
    ~FField_MDBAS();

    virtual double getE();
    virtual double getE(bool useVect);

private:
    
#ifdef VECTORIZED_ENER_EXPERIMENTAL
    //vectorizable buffers
    double* vect_vdw_6;
    double* vect_vdw_12;
    double* crds;
    double* q;
    double* e;
    double* s;
    
    double* rt;
    double* qij;
    double* eij;
    double* sij;
#endif

    virtual double getEtot(bool useVect);
    virtual double getEswitch(bool useVect);
        
    virtual void computeNonBonded_full();
    virtual void computeNonBonded14();
    
    virtual void computeNonBonded_switch();
    virtual void computeNonBonded14_switch();
    
    inline double computeEelec(const double qi, const double qj, const double rt);
    inline double computeEvdw(const double epsi, const double epsj, const double sigi,
                               const double sigj, const double rt);
    
#ifdef VECTORIZED_ENER_EXPERIMENTAL
    // vectorised versions
    virtual void computeNonBonded_full_VECT();
    virtual void computeNonBonded_switch_VECT();
    inline double computeEelec_VECT(double qij[], const double rt[], size_t len);
    inline double computeEvdw_VECT(double epsij[], double sigij[], const double rt[], size_t len, size_t offset);
    inline void computeEelec_VECT_SWITCH(double qij[], const double rt[], size_t len);
    inline void computeEvdw_VECT_SWITCH(double epsij[], double sigij[], const double rt[], size_t len, size_t offset);
#endif
    
#ifdef RANGED_E_EXPERIMENTAL
    virtual double computeNonBonded_full_range(int first, int last);
    virtual double computeNonBonded14_full_range(int first, int last);
#endif
    
    virtual void computeEbond();
    virtual void computeEang();
    virtual void computeEub();
    virtual void computeEdihe();
    virtual void computeEimpr();
    
    virtual double E_moving_set(int moveAtomList[], int moveBondList[]);

};

#endif	/* FFIELD_MDBAS_H */

