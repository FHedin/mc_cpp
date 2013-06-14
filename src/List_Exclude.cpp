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

#include <iostream>
#include <vector>
#include <chrono>

#include "List_Exclude.h"

#include "Atom.h"
#include "Angle.h"
#include "Bond.h"
#include "Bond_UB.h"
#include "Dihedral.h"
#include "Dihedral_improper.h"

using namespace std;

List_exclude::List_exclude(FField& _ff, Ensemble& _ens) : ff(_ff), ens(_ens)
{
    
    cout << "Building exclude list ..." << std::endl;
    
    auto start = chrono::system_clock::now();
    build_exclude_list();
    auto end = chrono::system_clock::now();
    auto elapsed_time =  chrono::duration_cast<chrono::milliseconds> (end-start).count();
    cout << "Time required for Exclude List was (milliseconds) : " << elapsed_time << endl;
    cerr << *this << endl;
    
    cout << "Building of exclude list done" << std::endl;
}

List_exclude::~List_exclude()
{
}

void List_exclude::resize_tempAtom(int ii, int jj)
{
    if (tmpPair[ii] >= nAlloc || tmpPair[jj] >= nAlloc)
    {
        nAlloc += nIncr;
        for (int j = 0; j < nAtom; j++)
            tempAtom[j].resize(nAlloc);
    }
}

void List_exclude::resize_tempConnect(int ii, int jj)
{
    if (tempConnectNum[ii] >= nConnect || tempConnectNum[jj] >= nConnect)
    {
        nConnect += nIncr;
        for (int j = 0; j < nAtom; j++)
            tempConnect[j].resize(nConnect);
    }
}

void List_exclude::resize_exclList(int idx)
{
    exclList[idx].resize(tmpPair[idx]);
}

void List_exclude::delete_all_temp()
{
    vector<int>().swap(tmpPair);
    vector<int>().swap(tempConnectNum);
    
    vector <vector<int>> ().swap(tempAtom);
    vector <vector<int>> ().swap(tempConnect);
    vector <vector<int>> ().swap(tempVer14);
}

void List_exclude::build_exclude_list()
{
    nAtom = ens.getN();
    nAlloc = 16;
    nIncr = 16;
    nConnect = 16;

    tmpPair = vector<int>(nAtom, 0);
    tempAtom = vector <vector<int>> (nAtom, vector<int>(nAlloc));

    tempConnectNum = vector<int>(nAtom, 0);
    tempConnect = vector <vector<int>> (nAtom, vector<int>(nConnect));

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
    exclList = vector <vector<int>> (nAtom, vector<int>());
    neighList14 = vector<int>(nPair14*2);

    // step 6 : build the real list from tmp arrays
    excl_final_Lists();

    // step 7 : delete temporary arrays
    delete_all_temp();
    
//    cout << nAtom << '\t' << nAlloc << '\t' << nIncr << '\t' << nConnect << endl;

}

void List_exclude::excl_bonds()
{
    int i, ia, ib, ii, jj;
    
    int nBond = ff.getNBond();
//    cout << "From List_exclude::excl_bonds() nBond is : " << nBond << endl;
    
    const vector<Bond>& bond = ff.getBndList();
    for (i = 0; i < nBond; i++)
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

void List_exclude::excl_angles()
{
    int i, j, ia, ib, ic, ii, jj;
    int exclude;
    
    int nAngle = ff.getNAngle();
//    cout << "From List_exclude::excl_angles() nAngle is : " << nAngle << endl;

    const vector<Angle>& angle = ff.getAngList();
    for (i = 0; i < nAngle; i++)
    {
        ia = angle[i].getAt1();
        ib = angle[i].getAt2();
        ic = angle[i].getAt3();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }
    }// end of angles job
}

void List_exclude::excl_dihedrals()
{
    int i, j, ia, ib, ic, id, ii, jj;
    int exclude;

    tempVer14 = vector < vector<int >> (5 * ff.getNDihedral(), vector<int>(2));
    nPair14 = 0;
    
    int nDihe = ff.getNDihedral();
//    cout << "From List_exclude::excl_dihedrals() nDihe is : " << nDihe << endl;

    const vector<Dihedral>& dihe = ff.getDiheList();
    for (i = 0; i < nDihe; i++)
    {
        ia = dihe[i].getAt1();
        ib = dihe[i].getAt2();
        ic = dihe[i].getAt3();
        id = dihe[i].getAt4();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }
    } // end of dihedrals job
}

void List_exclude::excl_impropers()
{
    int i, j, ia, ib, ic, id, ii, jj;
    int exclude;
    
    int nImproper = ff.getNImproper();
//    cout << "From  List_exclude::excl_impropers() nImproper is : " << nImproper << endl;

    const vector<Dihedral_improper>& impr = ff.getImprList();
    for (i = 0; i < nImproper; i++)
    {
        ia = impr[i].getAt1();
        ib = impr[i].getAt2();
        ic = impr[i].getAt3();
        id = impr[i].getAt4();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii, jj);

        exclude = 1;
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
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
        for (j = 0; j < tmpPair[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][tmpPair[ii]] = jj;
            tempAtom[jj][tmpPair[jj]] = ii;
            tmpPair[ii]++;
            tmpPair[jj]++;
        }

    } // end of impropers job
}

void List_exclude::excl_connectivity()
{
    int j, ia, ib, ic, id, ii, jj, k, kk, l;
    int exclude;

    for (ia = 0; ia < nAtom; ia++)
    {
        for (j = 0; j < tempConnectNum[ia]; j++)
        {
            ib = tempConnect[ia][j];
            for (k = 0; k < tempConnectNum[ib]; k++)
            {
                ic = tempConnect[ib][k];

                if (ic == ia)
                    continue;

                ii = ia;
                jj = ic;

                resize_tempAtom(ii, jj);

                exclude = 1;
                for (kk = 0; kk < tmpPair[ii]; kk++)
                {
                    if (tempAtom[ii][kk] == jj)
                    {
                        exclude = 0;
                        break;
                    }
                }

                if (exclude)
                {
                    tempAtom[ii][tmpPair[ii]] = jj;
                    tempAtom[jj][tmpPair[jj]] = ii;
                    tmpPair[ii]++;
                    tmpPair[jj]++;
                }

                for (l = 0; l < tempConnectNum[ic]; l++)
                {
                    id = tempConnect[ic][l];

                    if (id == ia || id == ib)
                        continue;

                    ii = ia;
                    jj = id;

                    resize_tempAtom(ii, jj);

                    //                    if (neigh->nPair14 >= 5 * param->nDihedral)
                    //                        my_error(DIHE_NONPARAM_ERROR, __FILE__, __LINE__, 0);

                    exclude = 1;
                    for (kk = 0; kk < tmpPair[ii]; kk++)
                    {
                        if (tempAtom[ii][kk] == jj)
                        {
                            exclude = 0;
                            break;
                        }
                    }

                    if (exclude)
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

void List_exclude::excl_final_Lists()
{
    int i = 0, ii = 0;
    int j = 0, jj = 0;
    int k = 0, kk = 0;

    int exclAtom;
    int hnAtom, hm1nAtom;

    hnAtom = nAtom / 2;
    hm1nAtom = (nAtom - 1) / 2;

    for (i = 0; i < nAtom; i++)
    {
        resize_exclList(i);

        jj = 0;
        for (j = 0; j < tmpPair[i]; j++)
        {
            exclAtom = tempAtom[i][j];

            if (((exclAtom > i) && ((exclAtom - i) < hnAtom)) || ((exclAtom < i) && ((exclAtom - i + nAtom) < hm1nAtom)))
            {
                exclList[ii][jj] = exclAtom;
                if (jj > 0)
                {
                    for (k = jj; k > 0; k--)
                    {
                        if (exclList[ii][k]<exclList[ii][k - 1])
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
    for (i = 0; i < nAtom; i++)
    {
        for (jj = 0; jj < exclPair[ii]; jj++)
        {
            if (exclList[ii][0] < i)
            {
                exclAtom = exclList[ii][0];

                for (kk = 0; kk < exclPair[ii] - 1; kk++)
                {
                    exclList[ii][kk] = exclList[ii][kk + 1];
                }

                exclList[ii][exclPair[ii]] = exclAtom;
            }
        }
        ii++;
    } // end of second for on natom
    
    //1-4 neighbours list
    for (i = 0; i < nPair14; i++)
    {
        neighList14[2 * i] = tempVer14[i][0];
        neighList14[2 * i + 1] = tempVer14[i][1];
    }

}

const std::vector<std::vector<int>>& List_exclude::getExclList() const
{
    return exclList;
}

const std::vector<int>& List_exclude::getExclPair() const
{
    return exclPair;
}

const std::vector<int>& List_exclude::getNeighList14() const
{
    return neighList14;
}

int List_exclude::getNPair14() const
{
    return nPair14;
}

std::ostream& operator<<(std::ostream& overloadStream, const List_exclude& exlst)
{
    exlst.toString(overloadStream);
    
    return overloadStream;
}

void List_exclude::toString(std::ostream& stream) const
{
    for(int i=0; i < nAtom; i++)
    {
        stream << "exclPair[" << i << "] : " << exclPair[i] << endl;
        stream << "exclList " << endl;
        for(int j=0; j < exclPair[i]; j++)
        {
            stream << exclList[i][j] << '\t';
        }
        stream << endl << endl;
    }
    stream << endl << endl;
}
