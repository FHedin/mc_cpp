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

#ifndef MOVE_TRN_H
#define	MOVE_TRN_H

#include <vector>

#include "Global_include.h"

#include "Move.h"
#include "Atom.h"
#include "PerConditions.h"

class Move_TRN : public Move
{
public:
    Move_TRN();
    virtual ~Move_TRN() = 0;

    static void translate_set(std::vector<Atom>& at_List, PerConditions& pbc, int moveAtomList[],
                                   double rndx, double rndy, double rndz);

private:
    static void tr_a_b(std::vector<Atom>& at_List, PerConditions& pbc, int atFirst, int atLast,
                       double rndx, double rndy, double rndz);

};

#endif	/* MOVE_TRN_H */

