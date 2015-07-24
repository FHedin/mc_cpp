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

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

#include "List_nonBonded.hpp"

#include "AtomList.hpp"
#include "Angle.hpp"
#include "Bond.hpp"
#include "Bond_UB.hpp"
#include "Dihedral.hpp"
#include "Dihedral_improper.hpp"
#include "FField.hpp"
#include "Tools.hpp"

#ifdef VECTORCLASS_EXPERIMENTAL
#include "vectorclass.h"
#endif

const double List_nonBonded::TOLLIST = 0.02;

#ifdef VECTORCLASS_EXPERIMENTAL
const VERLET_ALGORITHM verlet_type = BASIC_VECT;
#else
const VERLET_ALGORITHM verlet_type = BASIC;
#endif

using namespace std;

List_nonBonded::List_nonBonded(AtomList& _at_List, FField& _ff, PerConditions& _pbc,
                               Ensemble& _ens) : at_List(_at_List), ff(_ff), pbc(_pbc), ens(_ens)
{
    nAtom = ens.getN();

    cout << "Building exclude list ..." << std::endl;
    auto start = chrono::system_clock::now();
    build_exclude_list();
    auto end = chrono::system_clock::now();
    auto elapsed_time = chrono::duration_cast<chrono::milliseconds> (end - start).count();
    cout << "Building of exclude list done. ";
    cout << "Time required (milliseconds) : " << elapsed_time << endl;

    if (ff.getCutMode() != FULL)
    {
        cout << "Building verlet list";
        start = chrono::system_clock::now();

        switch(verlet_type)
        {
        case BASIC:
            cout << " using the standard list method." << endl;
            init_verlet_list();
            update_verlet_list();
            break;
#ifdef VECTORCLASS_EXPERIMENTAL
        case BASIC_VECT:
            cout << " using a VECTORIZED version of the standard list method." << endl;
            init_verlet_list();
            update_verlet_list_VECT();
            break;
#endif
#ifdef BALDRICH_EXPERIMENTAL
        case BALDRICH:
            cout << " using the B-Aldrich optimized method." << endl;
            init_verlet_list_BAldrich();
            update_verlet_list_BAldrich();
            break;
#endif
        default:
            init_verlet_list();
            update_verlet_list();
            break;
        }

        end = chrono::system_clock::now();
        elapsed_time = chrono::duration_cast<chrono::milliseconds> (end - start).count();
        cout << "Building of verlet list done. ";
        cout << "Time required (milliseconds) : " << elapsed_time << endl;
//         cout << *this << endl;
    }

}

List_nonBonded::~List_nonBonded()
{
}

void List_nonBonded::resize_tempAtom(int ii, int jj)
{
    if ( tmpPair[ii] >= nAlloc || tmpPair[jj] >= nAlloc )
    {
        nAlloc += nIncr;
        for ( int j = 0; j < nAtom; j++ )
            tempAtom[j].resize(nAlloc);
    }
}

void List_nonBonded::resize_tempConnect(int ii, int jj)
{
    if ( tempConnectNum[ii] >= nConnect || tempConnectNum[jj] >= nConnect )
    {
        nConnect += nIncr;
        for ( int j = 0; j < nAtom; j++ )
            tempConnect[j].resize(nConnect);
    }
}

void List_nonBonded::resize_exclList(int idx)
{
    exclList[idx].resize(tmpPair[idx]);
}

void List_nonBonded::delete_all_temp()
{
    vector<int>().swap(tmpPair);
    vector<int>().swap(tempConnectNum);

    vector < vector<int >> ().swap(tempAtom);
    vector < vector<int >> ().swap(tempConnect);
    vector < vector<int >> ().swap(tempVer14);
}

void List_nonBonded::build_exclude_list()
{
    nAlloc = 16;
    nIncr = 16;
    nConnect = 16;

    tmpPair = vector<int>(nAtom, 0);
    tempAtom = vector < vector<int >> (nAtom, vector<int>(nAlloc));

    tempConnectNum = vector<int>(nAtom, 0);
    tempConnect = vector < vector<int >> (nAtom, vector<int>(nConnect));

    // step 1 : bond connectivity
    excl_bonds();

    // step 2 : angles
    excl_angles();

    // step 3 : dihedrals
    excl_dihedrals();

    // step 4 : impropers
    excl_impropers();

    // step 5 : connectivity
    excl_connectivity();

    exclPair = vector<int>(nAtom);
    exclList = vector < vector<int >> (nAtom, vector<int>());
    neighList14 = vector<int>(nPair14 * 2);

    // step 6 : build the real list from tmp arrays
    excl_final_Lists();

    // step 7 : delete temporary arrays
    delete_all_temp();

    //    cout << nAtom << '\t' << nAlloc << '\t' << nIncr << '\t' << nConnect << endl;

}

void List_nonBonded::excl_bonds()
{
    int i, ia, ib, ii, jj;

    int nBond = ff.getNBond();
    //    cout << "From List_nonBonded::excl_bonds() nBond is : " << nBond << endl;

    const vector<Bond>& bond = ff.getBndList();
    for ( i = 0; i < nBond; i++ )
    {
        ia = bond[i].getAt1();
        ib = bond[i].getAt2();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii, jj);
        resize_tempConnect(ii, jj);

        tempConnect[ii][tempConnectNum[ii]] = jj;
        tempConnect[jj][tempConnectNum[jj]] = ii;
        tempConnectNum[ii]++;
        tempConnectNum[jj]++;

        tempAtom[ii][tmpPair[ii]] = jj;
        tempAtom[jj][tmpPair[jj]] = ii;
        tmpPair[ii]++;
        tmpPair[jj]++;
    }// end of bonds job
}

void List_nonBonded::excl_angles()
{
    int i, j, ia, ib, ic, ii, jj;
    int exclude;

    int nAngle = ff.getNAngle();
    //    cout << "From List_nonBonded::excl_angles() nAngle is : " << nAngle << endl;

    const vector<Angle>& angle = ff.getAngList();
    for ( i = 0; i < nAngle; i++ )
    {
        ia = angle[i].getAt1();
        ib = angle[i].getAt2();
        ic = angle[i].getAt3();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ia;
        jj = ic;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ib;
        jj = ic;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }
    }// end of angles job
}

void List_nonBonded::excl_dihedrals()
{
    int i, j, ia, ib, ic, id, ii, jj;
    int exclude;

    tempVer14 = vector < vector<int >> (5 * ff.getNDihedral(), vector<int>(2));
    nPair14 = 0;

    int nDihe = ff.getNDihedral();
    //    cout << "From List_nonBonded::excl_dihedrals() nDihe is : " << nDihe << endl;

    const vector<Dihedral>& dihe = ff.getDiheList();
    for ( i = 0; i < nDihe; i++ )
    {
        ia = dihe[i].getAt1();
        ib = dihe[i].getAt2();
        ic = dihe[i].getAt3();
        id = dihe[i].getAt4();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ia;
        jj = ic;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ia;
        jj = id;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempVer14[nPair14][0] = ia;
            tempVer14[nPair14][1] = id;
            nPair14++;

            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ib;
        jj = ic;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ib;
        jj = id;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ic;
        jj = id;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }
    } // end of dihedrals job
}

void List_nonBonded::excl_impropers()
{
    int i, j, ia, ib, ic, id, ii, jj;
    int exclude;

    int nImproper = ff.getNImproper();
    //    cout << "From  List_nonBonded::excl_impropers() nImproper is : " << nImproper << endl;

    const vector<Dihedral_improper>& impr = ff.getImprList();
    for ( i = 0; i < nImproper; i++ )
    {
        ia = impr[i].getAt1();
        ib = impr[i].getAt2();
        ic = impr[i].getAt3();
        id = impr[i].getAt4();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ia;
        jj = ic;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ia;
        jj = id;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ib;
        jj = ic;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ib;
        jj = id;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

        ii = ic;
        jj = id;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for ( j = 0; j < tmpPair[ii]; j++ )
        {
            if ( tempAtom[ii][j] == jj )
            {
                exclude = 0;
                break;
            }
        }

        if ( exclude )
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

    } // end of impropers job
}

void List_nonBonded::excl_connectivity()
{
    int j, ia, ib, ic, id, ii, jj, k, kk, l;
    int exclude;

    for ( ia = 0; ia < nAtom; ia++ )
    {
        for ( j = 0; j < tempConnectNum[ia]; j++ )
        {
            ib = tempConnect[ia][j];
            for ( k = 0; k < tempConnectNum[ib]; k++ )
            {
                ic = tempConnect[ib][k];

                if ( ic == ia )
                    continue;

                ii = ia;
                jj = ic;

                resize_tempAtom(ii, jj);

                exclude = 1;
                for ( kk = 0; kk < tmpPair[ii]; kk++ )
                {
                    if ( tempAtom[ii][kk] == jj )
                    {
                        exclude = 0;
                        break;
                    }
                }

                if ( exclude )
                {
                    tempAtom[ii][tmpPair[ii]] = jj;
                    tempAtom[jj][tmpPair[jj]] = ii;
                    tmpPair[ii]++;
                    tmpPair[jj]++;
                }

                for ( l = 0; l < tempConnectNum[ic]; l++ )
                {
                    id = tempConnect[ic][l];

                    if ( id == ia || id == ib )
                        continue;

                    ii = ia;
                    jj = id;

                    resize_tempAtom(ii, jj);

                    //                    if (neigh->nPair14 >= 5 * param->nDihedral)
                    //                        my_error(DIHE_NONPARAM_ERROR, __FILE__, __LINE__, 0);

                    exclude = 1;
                    for ( kk = 0; kk < tmpPair[ii]; kk++ )
                    {
                        if ( tempAtom[ii][kk] == jj )
                        {
                            exclude = 0;
                            break;
                        }
                    }

                    if ( exclude )
                    {
                        tempVer14[nPair14][0] = ia;
                        tempVer14[nPair14][1] = id;
                        nPair14++;

                        tempAtom[ii][tmpPair[ii]] = jj;
                        tempAtom[jj][tmpPair[jj]] = ii;
                        tmpPair[ii]++;
                        tmpPair[jj]++;
                    }

                }
            }
        }
    } // end of loop on nAtom

}

void List_nonBonded::excl_final_Lists()
{
    int i = 0, ii = 0;
    int j = 0, jj = 0;
    int k = 0, kk = 0;

    int exclAtom;
    int hnAtom, hm1nAtom;

    hnAtom = nAtom / 2;
    hm1nAtom = (nAtom - 1) / 2;

    for ( i = 0; i < nAtom; i++ )
    {
        resize_exclList(i);

        jj = 0;
        for ( j = 0; j < tmpPair[i]; j++ )
        {
            exclAtom = tempAtom[i][j];

            if ( ((exclAtom > i) && ((exclAtom - i) < hnAtom)) || ((exclAtom < i) && ((exclAtom - i + nAtom) < hm1nAtom)) )
            {
                exclList[ii][jj] = exclAtom;
                if ( jj > 0 )
                {
                    for ( k = jj; k > 0; k-- )
                    {
                        if ( exclList[ii][k] < exclList[ii][k - 1] )
                        {
                            exclAtom = exclList[ii][k];
                            exclList[ii][k] = exclList[ii][k - 1];
                            exclList[ii][k - 1] = exclAtom;
                        }
                    }
                }
                jj++;
            }
        } // end of for on tmppair
        exclPair[ii] = jj;
        ii++;
    } // end of for on natom

    ii = 0;
    for ( i = 0; i < nAtom; i++ )
    {
        for ( jj = 0; jj < exclPair[ii]; jj++ )
        {
            if ( exclList[ii][0] < i )
            {
                exclAtom = exclList[ii][0];

                for ( kk = 0; kk < exclPair[ii] - 1; kk++ )
                {
                    exclList[ii][kk] = exclList[ii][kk + 1];
                }

                exclList[ii][exclPair[ii]] = exclAtom;
            }
        }
        ii++;
    } // end of second for on natom

    //1-4 neighbours list
    for ( i = 0; i < nPair14; i++ )
    {
        neighList14[2 * i] = tempVer14[i][0];
        neighList14[2 * i + 1] = tempVer14[i][1];
    }

}

void List_nonBonded::init_verlet_list()
{
    int i, j, k, exclude;
    double r2, cutnb2, di[3], dj[3];

    cutnb2 = ff.getCutoff() + ff.getDeltacut();
    cutnb2 *= cutnb2;

    sizeList=0;
    neighPair = vector<int>(nAtom, 0);
    neighOrder = vector<int>(nAtom, 0);

    for (i = 0; i < nAtom - 1; i++)
    {
        at_List.getCoords(i,di);

        for (j = i + 1; j < nAtom; j++)
        {
            at_List.getCoords(j,dj);

            r2 = Tools::distance2(di, dj, pbc);

            if (r2 <= cutnb2)
            {
                exclude = 0;

                if ( at_List.Is_frozen(i) && at_List.Is_frozen(j) )
                {
                    exclude = 1;
                }
                else
                {
                    for (k = 0; k < exclPair[i]; k++)
                    {
                        if (exclList[i][k] == j)
                        {
                            exclude = 1;
                            break;
                        }
                    }
                }

                if(!exclude)
                    neighPair[i]++;

            } // if r2 <= cutnb2
        } // second for loop
        neighOrder[i]=i;
    } // first for loop

    sizeList = *max_element(neighPair.begin(), neighPair.end());

    sizeList = (int) (sizeList * (1. + 2. * TOLLIST)) + 1;
    neighList = vector < vector<int >> (nAtom, vector<int>(sizeList));
    
    cout << "sizeList is " << sizeList << " for " << nAtom << " atoms " << endl;

}

void List_nonBonded::update_verlet_list()
{
    int i, j, k, exclude=0;
    double r2, cutnb2, di[3], dj[3];

    cutnb2 = ff.getCutoff() + ff.getDeltacut();
    cutnb2 *= cutnb2;

    neighPair = vector<int>(nAtom, 0);
    neighOrder = vector<int>(nAtom, 0);

#ifdef _OPENMP
    #pragma omp parallel default(shared) private(i,j,k,r2,di,dj,exclude)
    {
        #pragma omp for schedule(dynamic) nowait
#endif
        for (i = 0; i < nAtom; i++)
        {
            at_List.getCoords(i,di);

            for (j = i + 1; j < nAtom; j++)
            {

                at_List.getCoords(j,dj);

                r2 = Tools::distance2(di, dj, pbc);

                if (r2 <= cutnb2)
                {
                    exclude = 0;

                    if (at_List.Is_frozen(i) && at_List.Is_frozen(j))
                    {
                        exclude = 1;
                    }
                    else
                    {
                        for (k = 0; k < exclPair[i]; k++)
                        {
                            if (exclList[i][k] == j)
                            {
                                exclude = 1;
                                break;
                            }
                        }
                    }

                    if (!exclude)
                    {
//                         if (neighPair[i] >= sizeList)
//                         {
// #ifdef _OPENMP
//                             #pragma omp critical
//                             {
// #endif
//                                 cout << "Warning : List larger than estimated. Size increased from " << sizeList;
//                                 sizeList = (int)(sizeList * (1. + TOLLIST)) + 1;
//                                 cout << " to " << sizeList << endl;
//
//                                 for (int l = 0; l < nAtom; l++)
//                                 {
//                                     neighList[l].resize(sizeList, 0);
//                                 }
// #ifdef _OPENMP
//                             }
// #endif
//                         }

                        neighList[i][neighPair[i]] = j;
                        neighPair[i]++;
                    }

                } // if r2
            } // second for
            neighOrder[i] = i;
        } // first for
#ifdef _OPENMP
    } // END OF parallel zone
#endif
}

#ifdef VECTORCLASS_EXPERIMENTAL

// void List_nonBonded::init_verlet_list_VECT()
// {
//
//     Vec4d r2,xi,yi,zi,xj,yj,zj,dx,dy,dz;
//     Vec4db test1,test2,test3,exclude;
//     Vec4db frozi,frozj;
//
//     const Vec4d cutnb2 = square(Vec4d(ff.getCutoff()+ff.getDeltacut())) ;
//     const Vec4d ones(1.0);
//     const Vec4d zeroes(0.0);
//
//     const vector<double>& x = at_List.getXvect();
//     const vector<double>& y = at_List.getYvect();
//     const vector<double>& z = at_List.getZvect();
//     const vector<bool>& frozList = at_List.getFrozenList();
//
//     const size_t psize = 4;
//     size_t remaining,end;
//
//     sizeList=0;
//     neighPair = vector<int>(nAtom, 0);
//     neighOrder = vector<int>(nAtom, 0);
//
//     for (size_t i = 0; i < nAtom - 1; i++)
//     {
//         xi = Vec4d(x[i]);
//         yi = Vec4d(y[i]);
//         zi = Vec4d(z[i]);
//         frozi = Vec4db(frozList[i]);
//
//         remaining = (nAtom-(i+1))%psize;
//         end = nAtom - remaining;
//
//         for (size_t j = i + 1; j < end; j+=psize)
//         {
//
//             xj.load(x.data()+j);
//             yj.load(y.data()+j);
//             zj.load(z.data()+j);
//             frozj = Vec4db(frozList[j],frozList[j+1],frozList[j+2],frozList[j+3]);
//
//             dx = xi - xj;
//             dy = yi - yj;
//             dz = zi - zj;
//             pbc.applyPBC(dx,dy,dz);
//             r2 = square(dx) + square(dy) + square(dz);
//
//             test1 = (r2 <= cutnb2);
//             if (horizontal_or(test1))
//             {
//                 exclude = false;
//
//                 test2 = frozi && frozj;
//
//                 test3 = Vec4db( binary_search(exclList[i].begin(),exclList[i].end(),j),
//                                 binary_search(exclList[i].begin(),exclList[i].end(),j+1),
//                                 binary_search(exclList[i].begin(),exclList[i].end(),j+2),
//                                 binary_search(exclList[i].begin(),exclList[i].end(),j+3)
//                               );
//
//                 exclude = test2 || test3;
//
//                 const Vec4d addToNPair = select(exclude,ones,zeroes);
//                 for(size_t p=0; p<psize; p++)
//                   neighPair[i] += round(addToNPair[p]);
//
//             } // if r2 <= cutnb2
//
//         } // second for loop on j
//
//         neighOrder[i]=i;
//
//     } // first for loop on i
//
//     sizeList = *max_element(neighPair.begin(), neighPair.end());
//
//     sizeList = (int) (sizeList * (1. + 2. * TOLLIST)) + 1;
//     neighList = vector < vector<int >> (nAtom, vector<int>(sizeList));
//
// }

void List_nonBonded::update_verlet_list_VECT()
{

    Vec4d r2,xi,yi,zi,xj,yj,zj,dx,dy,dz;
    Vec4db test1,test2,test3,exclude;
    Vec4db frozi,frozj;

    const Vec4d cutnb2 = square(Vec4d(ff.getCutoff()+ff.getDeltacut())) ;
    const Vec4d ones(1.0);
    const Vec4d zeroes(0.0);
    const Vec4d inf(std::numeric_limits<double>::infinity());

    const vector<double>& x = at_List.getXvect();
    const vector<double>& y = at_List.getYvect();
    const vector<double>& z = at_List.getZvect();
    const vector<bool>& frozList = at_List.getFrozenList();

    const size_t psize = 4;
    size_t remaining,end;

    neighPair = vector<int>(nAtom, 0);
    neighOrder = vector<int>(nAtom, 0);

    for (size_t i = 0; i < nAtom; i++)
    {
        xi = Vec4d(x[i]);
        yi = Vec4d(y[i]);
        zi = Vec4d(z[i]);
        frozi = Vec4db(frozList[i]);

        remaining = (nAtom-(i+1))%psize;
        end = nAtom - remaining;
        
        sort(exclList[i].begin(),exclList[i].end());
        // removes all elements with the value 0
        exclList[i].erase( std::remove( exclList[i].begin(), exclList[i].end(), 0 ), exclList[i].end() ); 

        for (size_t j = i + 1; j < end; j+=psize)
        {

            xj.load(x.data()+j);
            yj.load(y.data()+j);
            zj.load(z.data()+j);
            frozj = Vec4db(frozList[j],frozList[j+1],frozList[j+2],frozList[j+3]);

            dx = xi - xj;
            dy = yi - yj;
            dz = zi - zj;
            pbc.applyPBC(dx,dy,dz);
            r2 = square(dx) + square(dy) + square(dz);

            test1 = (r2 <= cutnb2);
            if (horizontal_or(test1))
            {
                exclude = false;

                test2 = frozi && frozj;

                // TODO : possibility to skip test3 if test2 to true everywhere

                test3 = Vec4db( binary_search(exclList[i].begin(),exclList[i].end(),j),
                                binary_search(exclList[i].begin(),exclList[i].end(),j+1),
                                binary_search(exclList[i].begin(),exclList[i].end(),j+2),
                                binary_search(exclList[i].begin(),exclList[i].end(),j+3)
                );

                exclude = (!test1) || (test2 || test3);

                const Vec4d addToNList = select(exclude,zeroes,Vec4d(j,j+1,j+2,j+3));
                const Vec4d addToNPair = select(exclude,zeroes,ones);

                for(size_t m=0; m<psize; m++)
                {
                  neighList[i][neighPair[i]] = round(addToNList[m]);
                  neighPair[i] += round(addToNPair[m]);
                }


            } // if r2 <= cutnb2

        } // second for loop on j

        // TODO : the same from end to remaining
        
//         remaining=0;
        if(remaining>0)
        {
          
          dx = zeroes;
          dy = zeroes;
          dz = zeroes;
          r2 = zeroes;
          frozj = false;
          
          size_t j=end;
          for(size_t m=0; m<remaining; m++)
          {
            dx.insert(m, x[i]-x[j+m]);
            dy.insert(m, y[i]-y[j+m]);
            dz.insert(m, z[i]-z[j+m]);
            frozj.insert(m, frozList[j+m]);
          }

          pbc.applyPBC(dx,dy,dz);
          r2 = square(dx) + square(dy) + square(dz);
          r2 = select(r2==0,inf,r2);
          
          test1 = (r2 <= cutnb2);
          if (horizontal_or(test1))
          {
            exclude = false;
            
            test2 = frozi && frozj;
            
            // TODO : possibility to skip test3 if test2 to true everywhere
            test3 = Vec4db( binary_search(exclList[i].begin(),exclList[i].end(),j),
                            binary_search(exclList[i].begin(),exclList[i].end(),j+1),
                            binary_search(exclList[i].begin(),exclList[i].end(),j+2),
                            binary_search(exclList[i].begin(),exclList[i].end(),j+3)
            );
            
            exclude = (!test1) || (test2 || test3);
            
            const Vec4d addToNList = select(exclude,zeroes,Vec4d(j,j+1,j+2,j+3));
            const Vec4d addToNPair = select(exclude,zeroes,ones);
            
            for(size_t m=0; m<psize; m++)
            {
              neighList[i][neighPair[i]] = round(addToNList[m]);
              neighPair[i] += round(addToNPair[m]);
            }
            
            
          } // if r2 <= cutnb2
          
        } // second loop on remainings

        neighOrder[i]=i;

    } // first for loop on i

}
#endif

#ifdef BALDRICH_EXPERIMENTAL

void List_nonBonded::init_verlet_list_BAldrich()
{
    int i, ii, j, m, latm;
    int exclude;
    double r2, cutnb2, delta[3], di[3], dj[3];

    cutnb2 = ff.getCutoff() + ff.getDeltacut();
    cutnb2 *= cutnb2;

    latm = nAtom;
    int hnAtom = nAtom / 2;
    int hm1nAtom = (nAtom - 1) / 2;

    counter.resize(nAtom, 0);
    neighPair.resize(nAtom, 0);

    for ( m = 0; m < hnAtom; m++ )
    {
        if ( m >= hm1nAtom )
            latm = hnAtom;

        ii = 0;

        for ( i = 0; i < latm; i++ )
        {
            j = i + m + 1;

            if ( j >= nAtom )
                j = j - nAtom;

            cout << "i : " << i << "\t ii : " << ii << "\t exclPair[ii] : " << exclPair.at(ii);
            cout << "\t size of exclist[ii] : " << exclList.at(ii).size();
            cout << "\t counter[ii] : " << counter.at(ii) << endl;
            cout.flush();

            if ( (exclPair.at(ii) > 0) && (exclList.at(ii).at(counter.at(ii)) == j) )
//             if ( (exclPair[ii] > 0) && (exclList[ii][counter[ii]] == j) )
            {
                counter[ii]++;
            }
            else
            {
                exclude = 0;
                if ( at_List[i].Is_frozen() && at_List[j].Is_frozen() )
                    exclude = 1;

                if ( !exclude )
                {
                    at_List[i].getCoords(di);
                    at_List[j].getCoords(dj);

                    r2 = Tools::distance2(di, dj, pbc, delta);

                    //                    cout << i << '\t' << j << '\t' << r2 << '\t' << cutnb2 << endl;
                    if ( r2 <= cutnb2 )
                    {
                        neighPair[ii]++;
                    }
                }
            }
            ii++;
        }
    } // end main loop

    sizeList = *max_element(neighPair.begin(), neighPair.end());

    sizeList = (int) (sizeList * (1. + 2. * TOLLIST)) + 1;
    neighList = vector < vector<int >> (nAtom, vector<int>(sizeList, 0));

}

void List_nonBonded::update_verlet_list_BAldrich()
{
    int i, ii, j, k, l, m, latm;
    int exclude;
    double r2, cutnb2, delta[3], di[3], dj[3];

    cutnb2 = ff.getCutoff() + ff.getDeltacut();
    cutnb2 *= cutnb2;

    latm = nAtom;
    int hnAtom = nAtom / 2;
    int hm1nAtom = (nAtom - 1) / 2;

    counter = vector<int>(nAtom, 0);
    neighPair = vector<int>(nAtom, 0);

    for ( m = 0; m < hnAtom; m++ )
    {
        if ( m >= hm1nAtom )
            latm = hnAtom;

        ii = 0;

        for ( i = 0; i < latm; i++ )
        {
            j = i + m + 1;

            if ( j >= nAtom )
                j = j - nAtom;

            if ( (exclPair.at(ii) > 0) && (exclList.at(ii).at(counter[ii]) == j) )
//             if ( (exclPair[ii] > 0) && (exclList[ii][counter[ii]] == j) )
            {
                counter[ii]++;
            }
            else
            {
                exclude = 0;
                if ( at_List[i].Is_frozen() && at_List[j].Is_frozen() )
                    exclude = 1;

                if ( !exclude )
                {
                    at_List[i].getCoords(di);
                    at_List[j].getCoords(dj);

                    r2 = Tools::distance2(di, dj, pbc, delta);

                    //                    cout << i << '\t' << j << '\t' << r2 << '\t' << cutnb2 << endl;
                    if ( r2 <= cutnb2 )
                    {
                        if ( neighPair[ii] >= sizeList )
                        {
                            cout << "Warning : List larger than estimated. Size increased from " << sizeList;
                            sizeList = (int) (sizeList * (1. + TOLLIST)) + 1;
                            cout << " to " << sizeList << endl;

                            for ( l = 0; l < nAtom; l++ )
                            {
                                neighList[l].resize(sizeList,0);
                            }
                        }

                        neighList[ii][neighPair[ii]] = j;
                        neighPair[ii]++;

                    }
                }
            }
            ii++;
        } // end of loop on natom
    } // end of main loop
}

#endif

const std::vector<std::vector<int >> &List_nonBonded::getExclList() const
{
    return exclList;
}

const std::vector<int>& List_nonBonded::getExclPair() const
{
    return exclPair;
}

const std::vector<int>& List_nonBonded::getNeighList14() const
{
    return neighList14;
}

int List_nonBonded::getNPair14() const
{
    return nPair14;
}

const std::vector<int>& List_nonBonded::getNeighPair() const
{
    return neighPair;
}

const std::vector<int>& List_nonBonded::getNeighOrder() const
{
    return neighOrder;
}

const std::vector<std::vector<int>>& List_nonBonded::getNeighList() const
{
    return neighList;
}

std::ostream& operator<<(std::ostream& overloadStream, const List_nonBonded& exlst)
{
    exlst.toString(overloadStream);

    return overloadStream;
}

void List_nonBonded::toString(std::ostream& stream) const
{
    for ( int i = 0; i < nAtom; i++ )
    {
        stream << "exclPair[" << i << "] : " << exclPair[i] << endl;
        stream << "exclList " << endl;
        for ( int j = 0; j < exclPair[i]; j++ )
        {
            stream << exclList[i][j] << '\t';
        }
        stream << endl << endl;
    }
    stream << endl << endl;

    for ( int i = 0; i < nAtom; i++ )
    {
        stream << "neighPair[" << i << "] : " << neighPair[i] << endl;
        stream << "neighOrder[" << i << "] : " << neighOrder[i] << endl;
        stream << "neighList " << endl;
        for ( int j = 0; j < neighPair[i]; j++ )
        {
            stream << neighList[i][j] << '\t';
        }
        stream << endl << endl;
    }
    stream << endl << endl;
}
