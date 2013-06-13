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

#include <iomanip>

#include <cmath>

#include "FField_MDBAS.h"
#include "Tools.h"

using namespace std;

FField_MDBAS::FField_MDBAS(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : FField(_at_List, _pbc, _ens)
{
}

FField_MDBAS::~FField_MDBAS()
{
}

double FField_MDBAS::getEtot()
{
    /* --- Preliminar work should come here --- */
    cout << std::fixed << std::setprecision(15);
    /* --- */

    // electrostatic and vdw
    computeNonBonded_full();
    computeNonBonded14_full();
    cout << "Electrostatic Full (kcal/mol) : " << this->elec / FField::kcaltoiu << endl;
    cout << "Van der Waals Full (kcal/mol) : " << this->vdw / FField::kcaltoiu << endl;

    /* --- Other types of energies here --- */

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

    return tot;
}

void FField_MDBAS::computeNonBonded_full()
{
    int i, j, k, exclude;
    double lelec = 0., pelec; // delec;
    double levdw = 0., pvdw; // dvdw;
    double r, r2, rt; // fxi, fyi, fzi, fxj, fyj, fzj;
    double di[3], dj[3]; //,delta[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;

    int nAtom = ens.getN();

    const vector<int>& exclPair = excl->getExclPair();
    const vector < vector<int >> &exclList = excl->getExclList();

    for (i = 0; i < nAtom - 1; i++)
    {
        //        fxi = 0.;
        //        fyi = 0.;
        //        fzi = 0.;

        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon();
        sigi = at_List[i].getSigma();

        for (j = i + 1; j < nAtom; j++)
        {

            exclude = 0;
            for (k = 0; k < exclPair[i]; k++)
            {
                if (exclList[i][k] == j)
                {
                    exclude = 1;
                    break;
                }
            }

            at_List[j].getCoords(dj);
            qj = at_List[j].getCharge();
            epsj = at_List[j].getEpsilon();
            sigj = at_List[j].getSigma();

            if (!exclude)
            {
                /*
                delta[0] = x[j] - x[i];
                delta[1] = y[j] - y[i];
                delta[2] = z[j] - z[i];

                r2 = dist(box, delta);*/

                r2 = Atom::distance2(di, dj, pbc);

                r = sqrt(r2);
                rt = 1. / r;

                pelec = computeEelec(qi, qj, rt);
                pvdw = computeEvdw(epsi, epsj, sigi, sigj, r);

                lelec += pelec;
                levdw += pvdw;

                //                delec = -pelec*rt;

                //                fxj = delec * delta[0] * rt;
                //                fyj = delec * delta[1] * rt;
                //                fzj = delec * delta[2] * rt;
                //
                //                fxi += fxj;
                //                fyi += fyj;
                //                fzi += fzj;
                //
                //                fx[j] += -fxj;
                //                fy[j] += -fyj;
                //                fz[j] += -fzj;

            }

        }

        //        fx[i] += fxi;
        //        fy[i] += fyi;
        //        fz[i] += fzi;
    }

    this->elec += lelec;
    this->vdw += levdw;
}

void FField_MDBAS::computeNonBonded14_full()
{
    int i, j, k;
    double lelec = 0., pelec;
    double levdw = 0., pvdw;
    double r, r2, rt;
    double di[3], dj[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;

    int nPair14 = excl->getNPair14();

    const vector<int>& neighList14 = excl->getNeighList14();

    for (k = 0; k < nPair14; k++)
    {
        i = neighList14[2 * k];
        j = neighList14[2 * k + 1];

        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon14();
        sigi = at_List[i].getSigma14();

        at_List[j].getCoords(dj);
        qj = at_List[j].getCharge();
        epsj = at_List[j].getEpsilon14();
        sigj = at_List[j].getSigma14();

        r2 = Atom::distance2(di, dj, pbc);

        r = sqrt(r2);
        rt = 1. / r;

        pelec = computeEelec(qi, qj, rt);
        pvdw = computeEvdw(epsi, epsj, sigi, sigj, r);

        lelec += pelec;
        levdw += pvdw;
    }

    this->elec += lelec;
    this->vdw += levdw;
}

double FField_MDBAS::computeEelec(const double qi, const double qj, const double rt)
{
    return FField::chgcharmm * FField::kcaltoiu * qi * qj * rt;
}

double FField_MDBAS::computeEvdw(const double epsi, const double epsj, const double sigi,
        const double sigj, const double r)
{
    return 4. * epsi * epsj * (Tools::X12((sigi + sigj) / r) - Tools::X6((sigi + sigj) / r));
}

void FField_MDBAS::computeEbond()
{
}

void FField_MDBAS::computeEang()
{
}

void FField_MDBAS::computeEub()
{
}

void FField_MDBAS::computeEdihe()
{
}

void FField_MDBAS::computeEimpr()
{
}

