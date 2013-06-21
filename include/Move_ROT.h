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

#ifndef MOVE_ROT_H
#define	MOVE_ROT_H

#include <vector>

#include "Move.h"
#include "Atom.h"

class Move_ROT : public Move
{
public:
    Move_ROT();
    virtual ~Move_ROT() = 0;

    static void rotate_set(std::vector<Atom>& at_List, int moveAtomList[],
                           int Pivot, double rndAngle, double rndvec[3]);
    
private:
    static void find_rot_mat(double mat[3][3], double rvec[3], double angle);
    static void apply_rotation(std::vector<Atom>& at_List, double pivot[3], double rotmat[3][3],
                               int first, int last);
};

#endif	/* MOVE_ROT_H */

