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
//#include <iostream>

#include "FField.h"

// White (1999)
//const double FField::epsilon = 0.2497878566;
//const double FField::sigma = 3.3345;

// from  Frenkel & Smit book 
//const double FField::epsilon = 0.238;
//const double FField::sigma = 3.41;

// reduced units
//const double FField::epsilon = 1.0;
//const double FField::sigma = 1.0;
//const double FField::kb_ch = 1.0;
//const double FField::kb_si = 1.0;
//const double FField::NA = 1.0;

//const double FField::sigma3 = sigma*sigma*sigma;
//const double FField::sigma6 = sigma3*sigma3;
//const double FField::sigma12 = sigma6*sigma6;

//const double FField::kb_ch = 1.98719e-03;
//const double FField::kb_si = 1.3806504e-23;
//const double FField::NA = 6.02214129e23;

//const double FField::rc = 2.5*sigma;

//const double FField::rconstr = 4.0*sigma;
//const double FField::rrconstrsq = 1.0/(rconstr*rconstr);


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

FField::FField(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : at_List(_at_List), pbc(_pbc), ens(_ens)
{
    //    w=0.;
    //    extraEnergy.resize(_at_List.size(),0.0);
}

FField::~FField()
{
    //dtor
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

//void FField::resetE()
//{
//    e=0.;
//}

//void FField::resetW()
//{
//    w=0.;
//}

/*
// get total potential of the system
double FField::getLJ(bool dV)
{
    double crd1[3],crd2[3],cm[3];
    
    double dx1,dy1,dz1;
    double dx,dy,dz;
    double d2,d6,d12;
    
    double e=0.;
    double e6,e12;
    double LT,L,dL,dL6=1.,dL12=1.;
    
    double epsii,epsij,epsilon;
    double sigmi,sigmj,sigma,sigma6,sigma12;
//    double rc;
    
//    double dconstr;
//    double Vconstr=0.0;

    int n = ens.getN();
    
//    Atom::getCentreOfMass(at_List,cm,n);

//    for (int i=0 ; i<(n-1) ; i++)
    for (int i=0 ; i<n ; i++)
    {
        at_List.at(i).getCoords(crd1);
        epsii = at_List.at(i).getEpsilon();
        sigmi = at_List.at(i).getSigma();
        
        dx1 = crd1[0];
        dy1 = crd1[1];
        dz1 = crd1[2];
        
//        dconstr = (dx1-cm[0])*(dx1-cm[0]) + (dy1-cm[1])*(dy1-cm[1]) + (dz1-cm[2])*(dz1-cm[2]);
//        Vconstr += Vconstraint(dconstr);
//        extraEnergy.at(i) = Vconstraint(dconstr);
        
        for (int j=i+1 ; j<n ; j++)
        {
            at_List.at(j).getCoords(crd2);
            epsij = at_List.at(j).getEpsilon();
            sigmj = at_List.at(j).getSigma();
            
            epsilon = sqrt(epsii*epsij);
            sigma   = 0.5*(sigmi+sigmj);
            sigma6  = sigma*sigma*sigma*sigma*sigma*sigma;
            sigma12 = sigma6*sigma6;
//            std::cerr << epsilon <<" "<< sigma <<" "<< sigma6 <<" "<< sigma12 << std::endl;
        
            dx = crd2[0]-dx1;
            dy = crd2[1]-dy1;
            dz = crd2[2]-dz1;

            pbc.applyPBC(dx,dy,dz);
            
            dx *= dx ; dy *= dy ; dz *= dz ;
            d2 = dx+dy+dz;
            
//            if (d2 > rc*rc)
//                continue;
            
            d6 = d2*d2*d2;
            d12= d6*d6;

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
//            w += 24.0*epsilon*( 2.0*e12*dL12 + e6*dL6 );
        }
    }
    
//    std::cerr << "Vconstr (void FField::getLJV(bool dV)) : " << Vconstr << std::endl;
//    e += Vconstr;
    
    return e;
}

// get total potential of an alternate (i.e trial) system
double FField::getLJ(std::vector<Atom>& candidateVec, bool dV)
{
    double crd1[3],crd2[3],cm[3];
    
    double dx1,dy1,dz1;
    double dx,dy,dz;
    double d2,d6,d12;
    
    double e=0.;
    double e6,e12;
    double LT,L,dL,dL6=1.,dL12=1.;
    
    double epsii,epsij,epsilon;
    double sigmi,sigmj,sigma,sigma6,sigma12;
//    double rc;
    
//    double dconstr;
//    double Vconstr=0.0;

    int n = ens.getN();
    
//    Atom::getCentreOfMass(candidateVec,cm,n);

//    for (int i=0 ; i<(n-1) ; i++)
    for (int i=0 ; i<n ; i++)
    {
        candidateVec.at(i).getCoords(crd1);
        epsii = candidateVec.at(i).getEpsilon();
        sigmi = candidateVec.at(i).getSigma();
            
        dx1 = crd1[0];
        dy1 = crd1[1];
        dz1 = crd1[2];
        
//        dconstr = (dx1-cm[0])*(dx1-cm[0]) + (dy1-cm[1])*(dy1-cm[1]) + (dz1-cm[2])*(dz1-cm[2]);
//        Vconstr += Vconstraint(dconstr);
//        extraEnergy.at(i) = Vconstraint(dconstr);

        for (int j=i+1 ; j<n ; j++)
        {
            candidateVec.at(j).getCoords(crd2);
            epsij = candidateVec.at(j).getEpsilon();
            sigmj = candidateVec.at(j).getSigma();
            
            epsilon = sqrt(epsii*epsij);
            sigma   = 0.5*(sigmi+sigmj);
            sigma6  = sigma*sigma*sigma*sigma*sigma*sigma;
            sigma12 = sigma6*sigma6;
//            std::cerr << epsilon <<" "<< sigma <<" "<< sigma6 <<" "<< sigma12 << std::endl;
            
            dx = crd2[0]-dx1;
            dy = crd2[1]-dy1;
            dz = crd2[2]-dz1;

            pbc.applyPBC(dx,dy,dz);
            
            dx *= dx ; dy *= dy ; dz *= dz ;
            d2 = dx+dy+dz;
            
//            if (d2 > rc*rc)
//                continue;
            
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
//            w += 24.0*epsilon*( 2.0*e12*dL12 + e6*dL6 );
        }
    }
    
//    std::cerr << "Vconstr (void FField::getLJV(std::vector<Atom>& candidateVec, bool dV)) : " << Vconstr << std::endl;
//    e += Vconstr;
    
    return e;
}

//get potential for a given atom and all the others atoms from the box
double FField::getLJ(Atom const& newAt, int candidate, bool dV)
{
    double crd1[3],crd2[3],cm[3];
    
    double dx1,dy1,dz1;
    double dx,dy,dz;
    double d2,d6,d12;
    
    double e=0.;
    double e6,e12;
    double LT,L,dL,dL6=1.,dL12=1.;
    
    double epsii,epsij,epsilon;
    double sigmi,sigmj,sigma,sigma6,sigma12;
//    double rc;
    
//    double dconstr;
//    double Vconstr=0.0;
    
    Atom::getCentreOfMass(at_List,cm,ens.getN());
    newAt.getCoords(crd1);
    
    epsii = newAt.getEpsilon();
    sigmi = newAt.getSigma();
    
    dx1 = crd1[0];
    dy1 = crd1[1];
    dz1 = crd1[2];
    
//    std::cerr << "x & xcm : " << dx1 << "\t" << cm[0] << std::endl;
//    std::cerr << "y & ycm : " << dy1 << "\t" << cm[1] << std::endl;
//    std::cerr << "z & zcm : " << dz1 << "\t" << cm[2] << std::endl;
//    std::cerr << std::endl;
    
//    dconstr = (dx1-cm[0])*(dx1-cm[0]) + (dy1-cm[1])*(dy1-cm[1]) + (dz1-cm[2])*(dz1-cm[2]);
//    Vconstr = Vconstraint(dconstr);
//    extraEnergy.at(candidate) = Vconstraint(dconstr);

    for (std::vector<Atom>::iterator it = at_List.begin() ; it != at_List.end() ; ++it)
    {
        it->getCoords(crd2);
        
        if (it->getID() == candidate)
            continue;
        
        epsij = it->getEpsilon();
        sigmj = it->getSigma();
            
        epsilon = sqrt(epsii*epsij);
        sigma   = 0.5*(sigmi+sigmj);
        sigma6  = sigma*sigma*sigma*sigma*sigma*sigma;
        sigma12 = sigma6*sigma6;
//        std::cerr << epsilon <<" "<< sigma <<" "<< sigma6 <<" "<< sigma12 << std::endl;

        dx = crd2[0]-dx1;
        dy = crd2[1]-dy1;
        dz = crd2[2]-dz1;

        pbc.applyPBC(dx,dy,dz);
        dx *= dx ; dy *= dy ; dz *= dz ;
        d2 = dx+dy+dz;
        
//        if (d2 > rc*rc)
//                continue;
        
        d6 = d2*d2*d2;
        d12= d6*d6;

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
//        w += 24.0*epsilon*( 2.0*e12*dL12 + e6*dL6 );
    }
    
//    std::cerr << "Vconstr (void FField::getLJV(Atom const& newAt, int candidate, bool dV)) : " << Vconstr << std::endl;
//    e += Vconstr;
    
    return e;
}

// a constraint (confining) potential for avoiding "clusters evaporation"
//double FField::Vconstraint(double distToCM)
//{
//    double Vc = distToCM*rrconstrsq;
//    Vc = pow(Vc,10);
//    Vc *= epsilon;
//    
////    std::cerr << "d = " << distToCM << "\t v = " << Vc << std::endl;
//    
//    return Vc;
//}

//estimate pressure from virial ; unit is Pascal
//double FField::PressFromVirial(int each)
//{
//    double V = ens.getVol();
//    V *= 1.0e-30; //for m^3
//    
////    double wl = w;
////    resetW();
//    
//    int n = ens.getN();
//
//    wl /= each;
//
//    wl *= 1.0/(3.0*V);
//    wl *= 4184.0/NA; //unit conversion
//
//    //Estimation of pressure
//    return ( (n/V) * kb_si * ens.getTemp() + wl ) ;
//}

//double FField::get_E()
//{
//    return e;
//}
//
//double FField::get_W()
//{
//    return w;
//}

// LJ tail energy in units of epsilon
//double FField::tail_energy()
//{
//    const double pi=3.14159265359;
//    double rho=ens.getN()/ens.getVol();
//    
//    double u_tail=(8./3.)*pi*rho*epsilon*sigma3;
//    u_tail *= (1./3.) * ((sigma6*sigma3)/pow(rc,9.0)) - (sigma3/pow(rc,3.0));
//    
//    return u_tail;
//}

// LJ tail pressure in Pascal units.
//double FField::tail_pressure()
//{
//    const double pi=3.14159265359;
//    double rho=ens.getN()/ens.getVol();
//    
//    double p_tail=(16./3.)*pi*rho*rho*epsilon*sigma3;
//    p_tail *= (2./3.) * ((sigma6*sigma3)/pow(rc,9.0)) - (sigma3/pow(rc,3.0));
//    
//    return p_tail;
//}

//double FField::getExtraE(int candidate) const
//{
//    return extraEnergy.at(candidate);
//}
//
//double FField::getExtraE() const
//{
//    double vtot=0.;
//    
//    for (std::vector<double>::const_iterator it = extraEnergy.begin() ; it != extraEnergy.end() ; ++it)
//    {
//        vtot += *it;
//    }
//    
//    return vtot;
//}

 */

std::ostream& operator<<(std::ostream& overloadStream, const FField& forf)
{
    forf.toString(overloadStream);
    
    return overloadStream;
}

void FField::toString(std::ostream& stream) const
{
    int i;
    
    for(i=0 ; i<ens.getN(); i++ )
         stream << at_List.at(i) << std::endl;
        
    for(i=0 ; i<nBond; i++ )
        stream << bndList.at(i) << std::endl;
    
    for(i=0 ; i<nUb; i++ )
        stream << ubList.at(i) << std::endl;
    
    for(i=0 ; i<nAngle; i++ )
        stream << angList.at(i) << std::endl;
    
    for(i=0 ; i<nDihedral; i++ )
        stream << diheList.at(i) << std::endl;
    
    for(i=0 ; i<nImproper; i++ )
        stream << imprList.at(i) << std::endl;   
}

