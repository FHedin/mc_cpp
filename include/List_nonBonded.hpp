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

#ifndef LIST_NONBONDED_H
#define	LIST_NONBONDED_H

class FField;

#include <vector>

//#include "Global_include.hpp"

#include "AtomList.hpp"
#include "FField.hpp"

enum LIST_ALGORITHM
{
    BASIC=0,
#ifdef VECTORCLASS_EXPERIMENTAL
    BASIC_VECT=1,
#endif
#ifdef BALDRICH_EXPERIMENTAL
    BALDRICH=2
#endif
};

class List_nonBonded
{
    friend std::ostream& operator<<( std::ostream& overloadStream, const List_nonBonded& exlst );

public:
    List_nonBonded(AtomList& _at_List, FField& _ff, PerConditions& _pbc, Ensemble& _ens);
    ~List_nonBonded();

    const std::vector<std::vector<int >> &getExclList() const;
    const std::vector<int>& getExclPair() const;
    const std::vector<int>& getNeighList14() const;
    int getNPair14() const;

    const std::vector<int>& getNeighPair() const;
    const std::vector<int>& getNeighOrder() const;
    const std::vector<std::vector<int>>& getNeighList() const;

    void update_verlet_list();

#ifdef VECTORCLASS_EXPERIMENTAL
    void update_verlet_list_VECT();
#endif

#ifdef BALDRICH_EXPERIMENTAL
    void update_verlet_list_BAldrich();
#endif

private:
    AtomList& at_List;
    FField& ff;
    PerConditions& pbc;
    Ensemble& ens;

    int nAtom, nAlloc, nIncr, nConnect;

    std::vector<int> tmpPair; // equivalent to tmpPair[nAtom]
    std::vector<std::vector<int>> tempAtom; // equivalent to tempAtom[nAtom][nAlloc]

    // tempConnectNum : connectivity for each bond (i.e. who's connected  with who))
    // tempConnect : for each member of tempConnectNum, list of directly bonded atom
    std::vector<int> tempConnectNum;
    std::vector<std::vector<int>> tempConnect;

    std::vector<std::vector<int>> tempVer14;
    int nPair14;

    // exclude list for
    std::vector<int> exclPair;
    std::vector<std::vector<int>> exclList;
    std::vector<int> neighList14;

    // neighbours list
    static const double TOLLIST;
    int sizeList;
    std::vector<int> neighPair;
    std::vector<int> neighOrder;
    std::vector<std::vector<int>> neighList;

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

    void init_verlet_list();

#ifdef VECTORCLASS_EXPERIMENTAL
    void init_verlet_list_VECT();
#endif

#ifdef BALDRICH_EXPERIMENTAL
    std::vector<int> counter;
    void init_verlet_list_BAldrich();
#endif

    void toString(std::ostream& stream) const;

};

#endif	/* LIST_NONBONDED_H */

