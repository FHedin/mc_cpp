/* 
 * File:   IO.h
 * Author: hedin
 *
 * Created on June 7, 2013, 5:17 PM
 */

#ifndef IO_H
#define	IO_H

#include "Atom.h"
#include "Ensemble.h"
#include "PerConditions.h"

class IO
{
public:
    IO(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
    virtual ~IO();
    
protected:
    std::vector<Atom>& at_List;
    PerConditions& pbc;
    Ensemble& ens;
    
    virtual void read_CONF()=0;

};

#endif	/* IO_H */

