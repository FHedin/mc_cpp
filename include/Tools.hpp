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

#ifndef TOOLS_H
#define	TOOLS_H

// class PerConditions;

#include <vector>
#include <string>
#include <algorithm> //for std::remove_if
#include <functional> //for std::ptr_fun

#include <cstring>

// #include <scorep/SCOREP_User.h>

//#include "Global_include.hpp"

#include "AtomList.hpp"
#include "PerConditions.hpp"

namespace Tools
{
    inline void str_rm_blank_spaces(std::string& str)
    {
        str.erase(std::remove_if(str.begin(), str.end(), std::ptr_fun(isspace)), str.end());
    }
    
    inline void str_to_lower_case(std::string& str)
    {
        std::transform(str.begin(), str.end(), str.begin(), std::ptr_fun(tolower));
    }

    inline void str_to_lower_case(char* str)
    {
        size_t n = strlen(str);

        for ( size_t i = 0; i < n; i++ )
        {
            str[i] = (char) tolower(str[i]);
        }
    }

    // stores b-a in c
    // 3 FLOP
    inline void vec_substract(const double a[3], const double b[3], double c[3])
    {
//         SCOREP_USER_FUNC_BEGIN();
        c[0] = b[0] - a[0];
        c[1] = b[1] - a[1];
        c[2] = b[2] - a[2];
//         SCOREP_USER_FUNC_END();
    }

//     inline void getCentreOfMass(AtomList& at_List, double cmass[3])
//     {
//         double crd[3];
//         double mass, mtot = 0.0;
//         cmass[0] = cmass[1] = cmass[2] = 0.0;
// 
//         for ( auto& it : at_List )
//         {
//             it.getCoords(crd);
//             mass = it.getMass();
//             mtot += mass;
//             cmass[0] += mass * crd[0];
//             cmass[1] += mass * crd[1];
//             cmass[2] += mass * crd[2];
//         }
//         cmass[0] /= mtot;
//         cmass[1] /= mtot;
//         cmass[2] /= mtot;
//     }

    inline void getCentreOfMass(AtomList& at_List, int first, int last, double cmass[3])
    {
        double crd[3];
        double mass, mtot = 0.0;
        cmass[0] = cmass[1] = cmass[2] = 0.0;

        for ( int i = first; i <= last; i++ )
        {
            at_List[i].getCoords(crd);
            mass = at_List[i].getMass();
            mtot += mass;
            cmass[0] += mass * crd[0];
            cmass[1] += mass * crd[1];
            cmass[2] += mass * crd[2];
        }
        cmass[0] /= mtot;
        cmass[1] /= mtot;
        cmass[2] /= mtot;
    }
    
    inline void getCentreOfMass(AtomList& at_List, int moveAtmList[], double cmass[3])
    {
        double crd[3];
        double mass, mtot = 0.0;
        cmass[0] = cmass[1] = cmass[2] = 0.0;

        int ng = moveAtmList[0];
        int endng = ng + 2;
        int nn;
        int iaf, ial;
        for ( int it1 = 1; it1 <= ng; it1++ )
        {
            nn = moveAtmList[it1];
            for ( int it2 = endng; it2 <= nn; it2 += 2 )
            {
                iaf = moveAtmList[it2 - 1];
                ial = moveAtmList[it2];
                for ( int it3 = iaf; it3 <= ial; it3++ )
                {
                    at_List[it3].getCoords(crd);
                    mass = at_List[it3].getMass();
                    mtot += mass;
                    cmass[0] += mass * crd[0];
                    cmass[1] += mass * crd[1];
                    cmass[2] += mass * crd[2];
                }
            }
            endng = nn + 2;
        }

        cmass[0] /= mtot;
        cmass[1] /= mtot;
        cmass[2] /= mtot;
    }

    // 23 FLOP
    inline double distance2(const double a1[3], const double a2[3],
                            const PerConditions& pbc)
    {
//         SCOREP_USER_FUNC_BEGIN();
        double r2;
        double delta[3];

        // 3FLOP
        Tools::vec_substract(a1, a2, delta);

        // 15 FLOP
        pbc.applyPBC(delta);

        // 5 FLOP
        r2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

//         SCOREP_USER_FUNC_END();
        
        return r2;
    }
    
    // 23 FLOP
    inline double distance2(const double a1[3], const double a2[3],
                            const PerConditions& pbc, double delta[3])
    {
        double r2;

        // 3 FLOP
        Tools::vec_substract(a1, a2, delta);

        // 15 FLOP
        pbc.applyPBC(delta);

        // 5 FLOP
        r2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

        return r2;
    }

    template <typename T>
    inline T X2(T x)
    {
        return (( x )*( x ) );
    }
    
    template <typename T>
    inline T X3(T x)
    {
        return (( x )*( x )*( x ) );
    }

    template <typename T>
    inline T X4(T x)
    {
        return (( x )*( x )*( x )*( x ) );
    }

    template <typename T>
    inline T X6(T x)
    {
        return (( x )*( x )*( x )*( x )*( x )*( x ) );
    }

    template <typename T>
    inline T X12(T x)
    {
        return (( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x ) );
    }

    // 2 FLOP
    inline double rint(double x)
    {
        return (double) (x >= 0. ? (int) (x + 0.5) : (int) (x - 0.5));
    }

}

// #ifdef VECTORIZED_ENER_EXPERIMENTAL
// 
// #ifdef __GNUC__ || __clang__
// #define RESTRICT __restrict__
// #elif _MSC_VER
// #define RESTRICT __restrict
// #endif
// 
// // A collection of vectorizable and inlined functions for fast floating point operations
// namespace Vectorized_Tools
// {
//     // adds values of 'from' to 'to'
//     // condition : arrays are distinct, and not overlapping
// 	inline void fast_double_add(double* RESTRICT to, double* RESTRICT from, const size_t len)
//     {
// #ifdef __GNUC__ || __clang__
//         double *x = (double*)__builtin_assume_aligned(to, 16);
//         double *y = (double*)__builtin_assume_aligned(from, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = to;
// 		__declspec(align(16)) double *y = from;
// #endif
//         for (size_t i = 0; i < len; i++)
//         {
//             x[i] += y[i];
//         } 
//     }
//     
//     // substracts values of 'from' to 'to'
//     // condition : arrays are distinct, and not overlapping
// 	inline void fast_double_sub(double* RESTRICT to, double* RESTRICT from, const size_t len)
//     {  
// #ifdef __GNUC__ || __clang__
// 		double *x = (double*)__builtin_assume_aligned(to, 16);
// 		double *y = (double*)__builtin_assume_aligned(from, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = to;
// 		__declspec(align(16)) double *y = from;
// #endif
//         for (size_t i = 0; i < len; i++)
//         {
//             x[i] -= y[i];
//         }
//     }
//     
//     // computes c = b - a
//     // condition : arrays are distinct, and not overlapping
// 	inline void fast_double_sub(double* RESTRICT a, double* RESTRICT b,
// 		double* RESTRICT c, const size_t len)
//     { 
// #ifdef __GNUC__ || __clang__
//         double *x = (double*)__builtin_assume_aligned(a, 16);
//         double *y = (double*)__builtin_assume_aligned(b, 16);
//         double *z = (double*)__builtin_assume_aligned(c, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = a;
// 		__declspec(align(16)) double *y = b;
// 		__declspec(align(16)) double *z = c;
// #endif
//         for (size_t i = 0; i < len; i++)
//         {
//             z[i] = y[i] - x[i];
//         }
//     }
//     
//     // multiplies values of 'to' by 'from'
//     // condition : arrays are distinct, and not overlapping
// 	inline void fast_double_mul(double* RESTRICT to, const double* RESTRICT from, const size_t len)
//     {  
// #ifdef __GNUC__ || __clang__
// 		double *x = (double*)__builtin_assume_aligned(to, 16);
// 		double *y = (double*)__builtin_assume_aligned(from, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = to;
// 		__declspec(align(16)) const double *y = from;
// #endif
//         for (size_t i = 0; i < len; i++)
//         {
//             x[i] *= y[i];
//         }
//     }
//     
//     // divides values of 'to' by 'from'
//     // condition : arrays are distinct, and not overlapping
// 	inline void fast_double_div(double* RESTRICT to, double* RESTRICT from, const size_t len)
//     {    
// #ifdef __GNUC__ || __clang__
// 		double *x = (double*)__builtin_assume_aligned(to, 16);
// 		double *y = (double*)__builtin_assume_aligned(from, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = to;
// 		__declspec(align(16)) double *y = from;
// #endif
//         for (size_t i = 0; i < len; i++)
//         {
//             x[i] /= y[i];
//         }
//     }
//     
//     // stores in 'to' the sqrt of 'from'
//     // condition : arrays are distinct, and not overlapping
// 	inline void fast_double_sqrt(double* RESTRICT to, double* RESTRICT from, const size_t len)
//     {
// #ifdef __GNUC__ || __clang__
// 		double *x = (double*)__builtin_assume_aligned(to, 16);
// 		double *y = (double*)__builtin_assume_aligned(from, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = to;
// 		__declspec(align(16)) double *y = from;
// #endif
// 
//         for (size_t i = 0; i < len; i++)
//         {
//             x[i] = sqrt(y[i]);
//         }
//     }
//     
//     // computes sqrt of array dat and overwrites it with new computed values
//     inline void fast_double_sqrt(double* RESTRICT dat, const size_t len)
//     {
// #ifdef __GNUC__ || __clang__
// 		double *x = (double*)__builtin_assume_aligned(dat, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = dat;
// #endif
//         for (size_t i=0;i<len;i++)
//         {
//             x[i] = sqrt(x[i]);
//         }
//     }
//     
//     // inverts an array, i.e. for an array dat, returns 1.0/dat
// 	inline void fast_double_invert_array(double* RESTRICT dat, const size_t len)
//     {
// #ifdef __GNUC__ || __clang__
// 		double *x = (double*)__builtin_assume_aligned(dat, 16);
// #elif _MSC_VER
// 		__declspec(align(16)) double *x = dat;
// #endif
//         for (size_t i=0;i<len;i++)
//         {
//             x[i] = 1.0/x[i];
//         }
//     }
//     
// }// end of namespace
// 
// #endif /* VECTORIZED_ENER_EXPERIMENTAL */

#endif	/* TOOLS_H */

