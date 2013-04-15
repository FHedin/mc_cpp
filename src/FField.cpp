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

#include <cmath>
#include "FField.h"

// White (1999)
//const double FField::epsilon = 0.2497878566;
//const double FField::sigma = 3.3345;

// from  Frenkel & Smit book 
const double FField::epsilon = 0.238;
const double FField::sigma = 3.41;

const double FField::sigma3 = sigma*sigma*sigma;
const double FField::sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
const double FField::sigma12 = sigma6*sigma6;
const double FField::kb_ch = 1.98719e-03;
const double FField::kb_si = 1.3806504e-23;
const double FField::NA = 6.02214129e23;

const double FField::rc=2.5*sigma;

FField::FField(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : at_List(_at_List),pbc(_pbc),ens(_ens)
{
    e=w=0.;
}

FField::~FField()
{
    //dtor
}

void FField::resetE()
{
    e=0.;
}

void FField::resetW()
{
    w=0.;
}

// get total potential of the system
void FField::getLJV(bool dV)
{
    double crd1[3],crd2[3];
    double d2,dx,dy,dz;
    double d6,d12;
    double e6,e12;
    double LT,L,dL,dL6=1.,dL12=1.;

    int n = ens.getN();

    for (int i=0 ; i<(n-1) ; i++)
    {
        at_List.at(i).getCoords(crd1);

        for (int j=i+1 ; j<n ; j++)
        {
            at_List.at(j).getCoords(crd2);

            dx = crd2[0]-crd1[0];
            dy = crd2[1]-crd1[1];
            dz = crd2[2]-crd1[2];

            pbc.applyPBC(dx,dy,dz);
            dx *= dx ; dy *= dy ; dz *= dz ;
            d2 = dx+dy+dz;
            
            if (d2 > rc*rc)
                continue;
            
            d6 = d2*d2*d2;
            d12=d6*d6;

            e6  = -1.0*sigma6/d6;
            e12 = sigma12/d12;

//            if (dV)
//            {
//                switch (ens.getType())
//                {
//                    case NPT:
//                        LT = ens.getVolT();
//                        L  = ens.getVol();
//                        dL = pow(LT,1.0/3.0)/pow(L,1.0/3.0);
//                        dL6 = dL*dL*dL*dL*dL*dL;
//                        dL12 = dL6 *dL6 ;
//                        break;
//                    default:
//                        break;
//                }
//            }

//            e += 4.0*epsilon*( sigma12/d12 - sigma6/d6 ) ;
//            w += 24.0*epsilon*( 2.0*sigma12/d12 - sigma6/d6 )  ;
            e += 4.0*epsilon*( e12*dL12 + e6*dL6);
            w += 24.0*epsilon*( 2.0*e12*dL12 + e6*dL6 );
        }
    }
}

//get potential for a given atom and all the others atoms from the box
void FField::getLJV(Atom const& newAt, int candidate, bool dV)
{
    double crd[3],crdw[3];
    double d2,dx,dy,dz;
    double d6,d12;
    double e6,e12;
    double LT,L,dL,dL6=1.,dL12=1.;

    newAt.getCoords(crd);

    for (std::vector<Atom>::iterator it = at_List.begin() ; it != at_List.end() ; ++it)
    {
        if (it->getID() == candidate)
            continue;

        it->getCoords(crdw);

        dx = crd[0]-crdw[0];
        dy = crd[1]-crdw[1];
        dz = crd[2]-crdw[2];

        pbc.applyPBC(dx,dy,dz);
        dx *= dx ; dy *= dy ; dz *= dz ;
        d2 = dx+dy+dz;
        
        if (d2 > rc*rc)
                continue;
        
        d6 = d2*d2*d2;
        d12=d6*d6;

        e6  = -1.0*sigma6/d6;
        e12 = sigma12/d12;

//        if (dV)
//        {
//            switch (ens.getType())
//            {
//                case NPT:
//                    LT = ens.getVolT();
//                    L  = ens.getVol();
//                    dL = pow(LT,1.0/3.0)/pow(L,1.0/3.0);
//                    dL6 = dL*dL*dL*dL*dL*dL;
//                    dL12 = dL6 *dL6 ;
//                    break;
//                default:
//                    break;
//            }
//        }

//        e += 4.0*epsilon*( sigma12/d12 - sigma6/d6 ) ;
//        w += 24.0*epsilon*( 2.0*sigma12/d12 - sigma6/d6 )  ;
        e += 4.0*epsilon*( e12*dL12 + e6*dL6);
        w += 24.0*epsilon*( 2.0*e12*dL12 + e6*dL6 );
    }
}

//estimate pressure from virial ; result returned in pascal
double FField::PressFromVirial(int each)
{
    double V = ens.getVol();
    V *= 1.0e-30; //for m^3

    int n = ens.getN();

    w /= each;

    w *= 1.0/(3.0*V);
    w *= 4184.0/NA; //unit conversion

    //Estimation of pressure
    return ( (n/V) * kb_si * ens.getTemp() + w ) ;
}

double FField::get_E(){return e;}
double FField::get_W() {return w;}

double FField::tail_energy()
{
    double pi=3.14159265359;
    double rho=ens.getN()/ens.getVol();
    
    double u_tail=(8./3.)*pi*rho*epsilon*sigma3;
    u_tail *= (1./3.) * ((sigma6*sigma3)/pow(rc,9.0)) - (sigma3/pow(rc,3.0));
    
    return u_tail;
}

double FField::tail_pressure()
{
    double pi=3.14159265359;
    double rho=ens.getN()/ens.getVol();
    
    double p_tail=(16./3.)*pi*rho*rho*epsilon*sigma3;
    p_tail *= (2./3.) * ((sigma6*sigma3)/pow(rc,9.0)) - (sigma3/pow(rc,3.0));
    
    return p_tail;
}

