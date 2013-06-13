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

#ifndef LIST_H
#define	LIST_H

#include <vector>

#include "FField.h"

class List_exclude
{
    friend std::ostream& operator<<(std::ostream& overloadStream, const List_exclude& exlst);
    
public:
    List_exclude(FField& _ff, Ensemble& _ens);
    ~List_exclude();
    
private:
    FField& ff;
    Ensemble& ens;
    
    int nAtom, nAlloc, nIncr, nConnect;
    
    std::vector<int> tmpPair;               // equivalent to tmpPair[nAtom]
    std::vector<std::vector<int>> tempAtom; // equivalent to tempAtom[nAtom][nAlloc]

    // tempConnectNum : connectivity for each bond (i.e. who's connected  whith who))
    // tempConnect : for each member of tempConnectNum, list of directly bonded atom
    std::vector<int> tempConnectNum;
    std::vector<std::vector<int>> tempConnect;
    
    std::vector<std::vector<int>> tempVer14;
    int nPair14;
    
    std::vector<int> exclPair; 
    std::vector<std::vector<int>> exclList;
    
    void resize_tempAtom(int ii, int ij);
    void resize_tempConnect(int ii, int jj);
    void resize_exclList(int idx);
    void delete_all_temp();
    
    void build_exclude_list();
    
    void excl_bonds();
    void excl_angles();
    void excl_dihedrals();
    void excl_impropers();
    void excl_connectivity();
    void excl_final_Lists();
    
    void toString(std::ostream& stream) const;
    
};

#endif	/* LIST_H */

