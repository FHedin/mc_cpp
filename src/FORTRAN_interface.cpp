#include <iostream>

#include "PerConditions.hpp"

#include "Ens_NVT.hpp"

#include "AtomList.hpp"

#include "FField_MDBAS.hpp"

#include "IO_MDBAS.hpp"

#include "List_nonBonded.hpp"

static PerConditions* pbc = nullptr;
static Ensemble* ens = nullptr;
static AtomList atomList;
static FField* ff = nullptr;
static List_nonBonded* exlst = nullptr;

extern "C" {

void set_pbc_(char* type, double* a, double* b,  double* c, double* alpha, double* beta, double *gamma)
{
  std::string pbtype(type);
  Tools::str_rm_blank_spaces(pbtype);
  Tools::str_to_lower_case(pbtype);
  pbc = new PerConditions(pbtype, *a, *b, *c, *alpha, *beta, *gamma);
}


void set_ens_(int* natom, double* Temp)
{
  ens = new Ens_NVT(*natom,pbc->computeVol(),*Temp);
}

void set_ff_and_cor_(char* cutoffmode, double* ctoff, double* cton, double* dcut, char* FF_file, char* COR_file)
{

  std::string cutMode(cutoffmode);
  Tools::str_rm_blank_spaces(cutMode);
  Tools::str_to_lower_case(cutMode);
  ff = new FField_MDBAS(atomList, *pbc, *ens, cutMode, *ctoff, *cton, *dcut);

  std::string corname(COR_file);
  std::string fffile(FF_file);
  Tools::str_rm_blank_spaces(corname);
  Tools::str_rm_blank_spaces(fffile);

  IO* io = new IO_MDBAS(corname, fffile, *ff, atomList, *pbc, *ens);
  exlst = new List_nonBonded(atomList, *ff, *pbc, *ens);
  ff->setExcl(*exlst);

  delete io;
  
}

// all energies returned in kcal/mol
// ener[0]=total energy
// ener[1]=potential energy
// ener[2]=kinetic energy
// ener[3]=electrostatic energy
// ener[4]=vdw energy
// ener[5]=bond energy
// ener[6]=angles energy
// ener[7]=urey bradley energy
// ener[8]=dihedrals energy
// ener[9]=impropers energy
void get_energy_(double ener[10])
{
  ff->getE(ener);
  for(size_t i=0; i<10; i++)
    ener[i] /= CONSTANTS::kcaltoiu;
}



void clean_()
{
  if(pbc != nullptr)
    delete pbc;

  if(ens != nullptr)
    delete ens;

  if(ff != nullptr)
    delete ff;
  
  if(exlst != nullptr)
    delete exlst;

}

}//extern C