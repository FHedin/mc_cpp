/* 
 * File:   IO.cpp
 * Author: hedin
 * 
 * Created on June 7, 2013, 5:17 PM
 */

#include "IO.h"

IO::IO(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : at_List(_at_List),pbc(_pbc),ens(_ens)
{
}

IO::~IO()
{
}

