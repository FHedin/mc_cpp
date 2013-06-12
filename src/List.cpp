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

#include <vector>

#include "List.h"

#include "Atom.h"
#include "Angle.h"
#include "Bond.h"
#include "Bond_UB.h"
#include "Dihedral.h"
#include "Dihedral_improper.h"

using namespace std;

List::List(FField& _ff, Ensemble& _ens) : ff(_ff), ens(_ens)
{
}

List::~List()
{
}

void List::resize_tempAtom(int ii, int jj)
{
    if (tmpPair[ii] >= nAlloc || tmpPair[jj] >= nAlloc)
    {
        nAlloc += nIncr;
        for (int j = 0; j < nAtom; j++)
            tempAtom[j].resize(nAlloc);
    }  
}

void List::resize_tempConnect(int ii, int jj)
{
    if (tempConnectNum[ii] >= nConnect || tempConnectNum[jj] >= nConnect)
    {
        nConnect += nIncr;
        for (int j = 0; j < nAtom; j++)
            tempConnect[j].resize(nConnect);
    }
}

void List::exclude_list()
{
    nAtom = ens.getN();
    nAlloc = 16;
    nIncr = 16;
    nConnect = 16;

    tmpPair  = vector<int>(nAtom, 0);
    tempAtom = vector<vector<int>>(nAtom, vector<int>(nAlloc) );
    
    tempConnectNum  = vector<int>(nAtom, 0);
    tempConnect     = vector<vector<int>>(nAtom, vector<int>(nConnect) );
    
    // step 1 : bond connectivity
    excl_bonds();

    // step 2 : angles 
    excl_angles();
    
    //step 3 : dihedrals
    excl_dihedrals();
    
} // end of function

void List::excl_bonds()
{
    int i,ia,ib,ii,jj;
    
    const vector<Bond>& bond = ff.getBndList();
    for (i = 0; i < ff.getNBond(); i++)
    {
        ia = bond[i].getAt1();
        ib = bond[i].getAt2();
        
        ii = ia;
        jj = ib;

        resize_tempAtom(ii,jj);
        resize_tempConnect(ii,jj);
        
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

void List::excl_angles()
{
    int i,j,ia,ib,ic,ii,jj;
    int exclude;
    
    const vector<Angle>& angle = ff.getAngList();
    for (i = 0; i < ff.getNAngle(); i++)
    {
        ia = angle[i].getAt1();
        ib = angle[i].getAt2();
        ic = angle[i].getAt3();

        ii = ia;
        jj = ib;
        
        resize_tempAtom(ii,jj);
        
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
        
        resize_tempAtom(ii,jj);
        
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
        
        resize_tempAtom(ii,jj);
        
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

void List::excl_dihedrals()
{
    int i,j,ia,ib,ic,id,ii,jj;
    int exclude;
    
    tempVer14 = vector<vector<int>>(5*ff.getNDihedral(), vector<int>(2) );
    nPair14 = 0;

    const vector<Dihedral>& dihe = ff.getDiheList();
    for (i = 0; i < ff.getNDihedral(); i++)
    {
        ia = dihe[i].getAt1();
        ib = dihe[i].getAt2();
        ic = dihe[i].getAt3();
        id = dihe[i].getAt4();

        ii = ia;
        jj = ib;

        resize_tempAtom(ii,jj);

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

        resize_tempAtom(ii,jj);

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

        resize_tempAtom(ii,jj);

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

        resize_tempAtom(ii,jj);

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

        resize_tempAtom(ii,jj);

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

        resize_tempAtom(ii,jj);

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
