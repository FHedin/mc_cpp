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
#include <algorithm> // for std::max
#include <limits> // for std::numeric_limits<double>
#include <chrono> // for precise timing

#include <cmath>
#include <string>

#include "FField_MDBAS.h"

#include "Tools.h"

using namespace std;

FField_MDBAS::FField_MDBAS(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : FField(_at_List, _pbc, _ens) {
}

FField_MDBAS::~FField_MDBAS() {
}

double FField_MDBAS::getEtot() {
  
    /* --- Preliminar work should come here --- */
    cout << std::fixed << std::setprecision(15);
    /* --- */

    // electrostatic and vdw are performed together for minimising computations
    auto start = chrono::system_clock::now();
    computeNonBonded_full();
    computeNonBonded14_full();
    auto end = chrono::system_clock::now();
    auto elapsed_time =  chrono::duration_cast<chrono::milliseconds> (end-start).count();
    cout << "Electrostatic Full (kcal/mol) : " << this->elec / FField::kcaltoiu << endl;
    cout << "Van der Waals Full (kcal/mol) : " << this->vdw / FField::kcaltoiu << endl;
    cout << "Time required for NonBonded was (milliseconds) : " << elapsed_time << endl;
    
    // all the components of potential energy
    if (nBond > 0)
        computeEbond();
    cout << "Bonds energy (kcal/mol) : " << this->bond / FField::kcaltoiu << endl;

    if (nAngle > 0)
        computeEang();
    cout << "Angles energy (kcal/mol) : " << this->ang / FField::kcaltoiu << endl;

    if (nUb > 0)
        computeEub();
    cout << "Urey Bradley energy (kcal/mol) : " << this->ub / FField::kcaltoiu << endl;

    if (nDihedral > 0)
        computeEdihe();
    cout << "Dihedrals Energy (kcal/mol) : " << this->dihe / FField::kcaltoiu << endl;

    if (nImproper > 0)
        computeEimpr();
    cout << "Impropers energy (kcal/mol) : " << this->impr / FField::kcaltoiu << endl;
    /* --- Other types of energies here --- */

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

    cout << "Total energy (kcal/mol) : " << this->tot / FField::kcaltoiu << endl;

    return tot;
}

void FField_MDBAS::computeNonBonded_full() {
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

    for (i = 0; i < nAtom - 1; i++) {
        //        fxi = 0.;
        //        fyi = 0.;
        //        fzi = 0.;

        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon();
        sigi = at_List[i].getSigma();

        for (j = i + 1; j < nAtom; j++) {

            exclude = 0;
            for (k = 0; k < exclPair[i]; k++) {
                if (exclList[i][k] == j) {
                    exclude = 1;
                    break;
                }
            }

            at_List[j].getCoords(dj);
            qj = at_List[j].getCharge();
            epsj = at_List[j].getEpsilon();
            sigj = at_List[j].getSigma();

            if (!exclude) {
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

    this->elec = lelec;
    this->vdw = levdw;
}

void FField_MDBAS::computeNonBonded14_full() {
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

    for (k = 0; k < nPair14; k++) {
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

double FField_MDBAS::computeEelec(const double qi, const double qj, const double rt) {
    return FField::chgcharmm * FField::kcaltoiu * qi * qj * rt;
}

double FField_MDBAS::computeEvdw(const double epsi, const double epsj, const double sigi,
        const double sigj, const double r) {
//    return 4. * epsi * epsj * (Tools::X12((sigi + sigj) / r) - Tools::X6((sigi + sigj) / r));
	return 4. * epsi * epsj * ( pow( ((sigi + sigj) / r),12) - pow( ((sigi + sigj) / r),6) );
}

void FField_MDBAS::computeEbond() {
    int i, j, ll;
    int type;
    double di[3], dj[3];
    double r0, k;
    double d;
    double ebond = 0.0;

    for (ll = 0; ll < nBond; ll++) {
        i = bndList[ll].getAt1();
        j = bndList[ll].getAt2();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        d = Atom::distance2(di, dj, pbc);
        d = sqrt(d);

        r0 = bndList[ll].getR0();
        k = bndList[ll].getK();
        type = bndList[ll].getType();

        switch (type) {
            case BHARM:
                ebond += 0.5 * k * Tools::X2(d - r0);
                break;

            case BMORSE:
            {
                double beta, morsea, morseb;
                beta = bndList[ll].getBeta();
                morsea = exp(-beta * (d - r0));
                morseb = Tools::X2(morsea);
                ebond += k * (morseb - 2. * morsea) + k;
            }
                break;

            default:
                ebond += 0.5 * k * Tools::X2(d - r0);
                break;
        }
    }
    this->bond = ebond;
}

void FField_MDBAS::computeEang() {
    int i, j, k, ll;
    double di[3], dj[3], dk[3], dab[3], dbc[3];
    double rab, rbc, rabt, rbct, cost, sint, theta;
    double kst, theta0;
    double eang = 0.0;

    const double dbl_epsilon = numeric_limits<double>::epsilon();

    for (ll = 0; ll < nAngle; ll++) {
        i = angList[ll].getAt1();
        j = angList[ll].getAt2();
        k = angList[ll].getAt3();
        kst = angList[ll].getK();
        theta0 = angList[ll].getTheta0();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        at_List[k].getCoords(dk);

        rab = Atom::distance2(di, dj, pbc, dab);
        rab = sqrt(rab);
        rabt = 1. / rab;

        rbc = Atom::distance2(dk, dj, pbc, dbc);
        rbc = sqrt(rbc);
        rbct = 1. / rbc;

        cost = (dab[0] * dbc[0] + dab[1] * dbc[1] + dab[2] * dbc[2]) / (rab * rbc);
        sint = max(dbl_epsilon, sqrt(1.0 - (cost * cost)));
        theta = acos(cost);

        eang += 0.5 * kst * Tools::X2(theta - theta0);
    }
    this->ang = eang;
}

void FField_MDBAS::computeEub() {
    int i, j, ll;
    double di[3], dj[3];
    double r0, k;
    double d;
    double ebond = 0.0;

    for (ll = 0; ll < nUb; ll++) {
        i = ubList[ll].getAt1();
        j = ubList[ll].getAt2();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        d = Atom::distance2(di, dj, pbc);
        d = sqrt(d);

        r0 = ubList[ll].getR0();
        k = ubList[ll].getK();

        ebond += 0.5 * k * Tools::X2(d - r0);
    }
    this->ub = ebond;
}

void FField_MDBAS::computeEdihe() {
    int i, j, k, l, ll;
    double di[3], dj[3], dk[3], dl[3];
    double dab[3], dbc[3], dcd[3];
    double pb[3], pc[3];
    double rbc, rpb, rpc, r2pb, r2pc;
    double pbpc, cosp, sinp, phi;
    double edihe = 0.;
    double kst, phi0, mult;
    int order, type;

    const double twopi = FField::PI;
    const double dbl_epsilon = numeric_limits<double>::epsilon();

    for (ll = 0; ll < nDihedral; ll++) {
        i = diheList[ll].getAt1();
        j = diheList[ll].getAt2();
        k = diheList[ll].getAt3();
        l = diheList[ll].getAt4();
        kst = diheList[ll].getK();
        phi0 = diheList[ll].getPhi0();
        mult = diheList[ll].getMult();
        order = diheList[ll].getOrder();
        type = diheList[ll].getType();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        at_List[k].getCoords(dk);
        at_List[l].getCoords(dl);

        Tools::vec_substract(dj, di, dab);
        pbc.applyPBC(dab);

        rbc = sqrt(Atom::distance2(dk, dj, pbc, dbc));

        Tools::vec_substract(dl, dk, dcd);
        pbc.applyPBC(dcd);

        // construct first dihedral vector
        pb[0] = dab[1] * dbc[2] - dab[2] * dbc[1];
        pb[1] = dab[2] * dbc[0] - dab[0] * dbc[2];
        pb[2] = dab[0] * dbc[1] - dab[1] * dbc[0];
        r2pb = Tools::X2(pb[0]) + Tools::X2(pb[1]) + Tools::X2(pb[2]);
        rpb = sqrt(r2pb);

        // construct second dihedral vector
        pc[0] = dbc[1] * dcd[2] - dbc[2] * dcd[1];
        pc[1] = dbc[2] * dcd[0] - dbc[0] * dcd[2];
        pc[2] = dbc[0] * dcd[1] - dbc[1] * dcd[0];
        r2pc = Tools::X2(pc[0]) + Tools::X2(pc[1]) + Tools::X2(pc[2]);
        rpc = sqrt(r2pc);

        // determine dihedral angle 
        pbpc = pb[0] * pc[0] + pb[1] * pc[1] + pb[2] * pc[2];
        cosp = pbpc / (rpb * rpc);
        sinp = (dbc[0]*(pc[1] * pb[2] - pc[2] * pb[1]) + dbc[1]*(pb[0] * pc[2] - pb[2] * pc[0]) +
                dbc[2]*(pc[0] * pb[1] - pc[1] * pb[0])) / (rpb * rpc * rbc);
        phi = atan2(sinp, cosp);

        // avoid singularity in sinp
        if (sinp >= 0.) {
            sinp = max(dbl_epsilon, fabs(sinp));
        } else {
            sinp = -(max(dbl_epsilon, fabs(sinp)));
        }

        // calculate potential energy
        switch (type) {
            case DCOS: // cosine dihedral
                edihe += kst * (1. + cos(mult * phi - phi0));
                break;

            case DHARM: // harmonic dihedral
                phi = phi - phi0;
                phi = phi - PerConditions::rint(phi / twopi) * twopi;
                edihe += 0.5 * kst * (phi * phi);
                break;

            default:
                edihe += kst * (1. + cos(mult * phi - phi0));
                break;
        }

    } // end of for loop on dihedrals
    this->dihe = edihe;
}

void FField_MDBAS::computeEimpr() {
    int i, j, k, l, ll;
    double di[3], dj[3], dk[3], dl[3];
    double dab[3], dbc[3], dcd[3];
    double pb[3], pc[3];
    double rbc, rpb, rpc, r2pb, r2pc;
    double pbpc, cosp, sinp, phi;
    double eimpr = 0.;
    double kst, phi0, mult;
    int order, type;

    const double twopi = FField::PI;
    const double dbl_epsilon = numeric_limits<double>::epsilon();

    for (ll = 0; ll < nImproper; ll++) {
        i = imprList[ll].getAt1();
        j = imprList[ll].getAt2();
        k = imprList[ll].getAt3();
        l = imprList[ll].getAt4();
        kst = imprList[ll].getK();
        phi0 = imprList[ll].getPhi0();
        mult = imprList[ll].getMult();
        order = imprList[ll].getOrder();
        type = imprList[ll].getType();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        at_List[k].getCoords(dk);
        at_List[l].getCoords(dl);

        Tools::vec_substract(dj, di, dab);
        pbc.applyPBC(dab);

        rbc = sqrt(Atom::distance2(dk, dj, pbc, dbc));

        Tools::vec_substract(dl, dk, dcd);
        pbc.applyPBC(dcd);

        // construct first dihedral vector
        pb[0] = dab[1] * dbc[2] - dab[2] * dbc[1];
        pb[1] = dab[2] * dbc[0] - dab[0] * dbc[2];
        pb[2] = dab[0] * dbc[1] - dab[1] * dbc[0];
        r2pb = Tools::X2(pb[0]) + Tools::X2(pb[1]) + Tools::X2(pb[2]);
        rpb = sqrt(r2pb);

        // construct second dihedral vector
        pc[0] = dbc[1] * dcd[2] - dbc[2] * dcd[1];
        pc[1] = dbc[2] * dcd[0] - dbc[0] * dcd[2];
        pc[2] = dbc[0] * dcd[1] - dbc[1] * dcd[0];
        r2pc = Tools::X2(pc[0]) + Tools::X2(pc[1]) + Tools::X2(pc[2]);
        rpc = sqrt(r2pc);

        // determine dihedral angle 
        pbpc = pb[0] * pc[0] + pb[1] * pc[1] + pb[2] * pc[2];
        cosp = pbpc / (rpb * rpc);
        sinp = (dbc[0]*(pc[1] * pb[2] - pc[2] * pb[1]) + dbc[1]*(pb[0] * pc[2] - pb[2] * pc[0]) +
                dbc[2]*(pc[0] * pb[1] - pc[1] * pb[0])) / (rpb * rpc * rbc);
        phi = atan2(sinp, cosp);

        // avoid singularity in sinp
        if (sinp >= 0.) {
            sinp = max(dbl_epsilon, fabs(sinp));
        } else {
            sinp = -(max(dbl_epsilon, fabs(sinp)));
        }

        // calculate potential energy
        switch (type) {
            case DCOS: // cosine dihedral
                eimpr += kst * (1. + cos(mult * phi - phi0));
                break;

            case DHARM: // harmonic dihedral
                phi = phi - phi0;
				phi = phi - PerConditions::rint(phi / twopi) * twopi;
                eimpr += 0.5 * kst * (phi * phi);
                break;

            default:
                eimpr += kst * (1. + cos(mult * phi - phi0));
                break;
        }

    } // end of for loop on dihedrals
    this->impr = eimpr;
}

