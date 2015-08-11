#include <iostream>
#include <random>

#include "PerConditions.hpp"

#include "Ens_NVT.hpp"

#include "AtomList.hpp"

#include "FField_MDBAS.hpp"

#include "IO_MDBAS.hpp"

#include "List_nonBonded.hpp"

#include "List_Moves.hpp"

static PerConditions* pbc = nullptr;
static Ensemble* ens = nullptr;
static AtomList atomList;
static FField* ff = nullptr;
static List_nonBonded* exlst = nullptr;
static List_Moves* mvlist = nullptr;

static std::random_device seed;
static std::mt19937 generator;
static std::uniform_real_distribution<double> distributionAlpha;
static std::uniform_real_distribution<double> distributionMove;


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

    /*
      * move_type is either "rot" for rotations or "trn" for translation
      * move_mode is "all", "residue" or "atom" : when applying a move it decides if we apply this to all
        atoms of the selection at once, or only to the residue containing the atom selection, or only to the atom selection.

      * note that move_type="rot" and  move_mode="atom" is impossible as we can't rotate just one atom

      * Then the atom selection is given using this selection node
      * sel_mode are :
        "all" or "none" which are explicit enough...
        "residue_id" which uses the column 2 of cor files, "residue_name" using column 3
        "atom_idx" which uses the column 1 of the cor file

      * sele="..." is simply the string containing the selection (unused for sel_mode "all" or "none")
    */
    void add_move_type_(char* move_type, char* move_mode, char* sel_mode, char* sele)
    {
        int natom = ens->getN();
        if(mvlist == nullptr)
        {
            mvlist = new List_Moves(atomList,*ff,natom);
            unsigned int lseed = seed();
            generator.seed(lseed);
            distributionAlpha = std::uniform_real_distribution<double>(0.0, 1.0);
            distributionMove = std::uniform_real_distribution<double>(-1.0, 1.0);
        }

        std::string mvtyp(move_type);
        std::string mvmode(move_mode);
        std::string smode(sel_mode);
        std::string selec(sele);
        mvlist->addNewMoveType(mvtyp,mvmode,smode,selec);

    }
    
    
    static double rndUnifMove_()
    {
      return distributionMove(generator);
    }
    
    static double rndUnifAlpha_()
    {
      return distributionAlpha(generator);
    }
    
    
    static double rndIntCandidate_(int _n)
    {
      return (int)_n * (rndUnifAlpha_());
    }

    static void rndSphere_(double rnd[3])
    {
      double rx, ry, rz;
      
      do
      {
        rx = rndUnifMove_();
        ry = rndUnifMove_();
        rz = rndUnifMove_();
      } while ((rx * rx + ry * ry + rz * rz) > 1.0);
      
      rnd[0] = rx;
      rnd[1] = ry;
      rnd[2] = rz;
    }
    
    static void scaleVec(double r[3], double dmax)
    {
      r[0] *= dmax;
      r[1] *= dmax;
      r[2] *= dmax;
    }

    
    void random_move_(int* type, double* maxMove)
    {
//       int natom = ens->getN();
      int imvtyp = 0;
      int imvatm = 0;
      int nmvtyp = mvlist->getNMoveTypes();
      
      const std::vector<int>& nMoveAt = mvlist->getNMoveAtm();
      const std::vector<MOVETYPE>& movetypList = mvlist->getMoveTypeList();
      const std::vector<int**>& moveAtomList = mvlist->getMoveAtomList();
      const std::vector<int**>& movePivotList = mvlist->getMovePivotList();
      
      // get a random movetype and moveinstance from this movetype
      //imvtyp = rndIntCandidate_(nmvtyp); // get an int between 0 and (nmvtyp-1)

      if(*type > nmvtyp)
      {
        std::cerr << "Error, you have provided a move index larger than the current number of registered allowed move which is : " << nmvtyp << std::endl;
      }
      
      imvtyp = ( *type - 1 );
      imvatm = rndIntCandidate_(nMoveAt[imvtyp]); // get an int between 0 and (nMoveAt[imvtyp]-1)

      double r[3] = {0.0, 0.0, 0.0};
      double rang = 0.0;
      double dmax = *maxMove;
      
      switch ( movetypList[imvtyp] )
      {
        case TRN:
        {
          rndSphere_(r);
          scaleVec(r,dmax);
          MOVE_TRN::translate_set(atomList, *pbc, moveAtomList[imvtyp][imvatm],
                                  r[0], r[1], r[2]);
          break;
        }
        case ROT:
        {
          rndSphere_(r);
          rang = dmax * rndUnifMove_();
          MOVE_ROT::rotate_set(atomList, *pbc, moveAtomList[imvtyp][imvatm],
                               movePivotList[imvtyp][imvatm][0], rang, r);
          break;
        }
        default:
          break;
      }
      
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

// sets new coordinates
    void set_coords_(double x[], double y[], double z[])
    {
        int natom = ens->getN();
        for (int i=0; i<natom; i++)
            atomList.setCoords(i,x[i],y[i],z[i]);
    }
    
    // get coordinates
    void get_coords_(double x[], double y[], double z[])
    {
      int natom = ens->getN();
      
      const double* X = atomList.getXvect().data();
      const double* Y = atomList.getYvect().data();
      const double* Z = atomList.getZvect().data();
      
      memcpy(x,X,natom * sizeof(double));
      memcpy(y,Y,natom * sizeof(double));
      memcpy(z,Z,natom * sizeof(double));
      
//       for (int i=0; i<natom; i++)
//         atomList.get
//         atomList.setCoords(i,x[i],y[i],z[i]);
    }

// stores atomic masses to provided mass array
    void get_mass_(double mass[])
    {
        int natom = ens->getN();
        for (int i=0; i<natom; i++)
            mass[i] = atomList.getMass(i);
    }
    
    void write_xyz_(char* filename)
    {
      double crd[3] = { 0. };
      int n = ens->getN();
      const char *symb = nullptr;
      
      FILE* xyz = NULL;
      xyz = fopen(filename,"at");
      if(xyz==NULL)
      {
        std::cerr << "Error opening file : " << filename << std::endl;
        exit(-98);
      }
      
      fprintf(xyz, "%d\n", n);
      fprintf(xyz, "# Dump of coordinates in XYZ format compatible with VMD\n");
      
      for (int it = 0; it < n ;it++)
      {
        atomList.getCoords(it,crd);
        symb = atomList.getSymbol(it).c_str();
        fprintf(xyz, "%s\t%10.5lf\t%10.5lf\t%10.5lf\n", symb, crd[0], crd[1], crd[2]);
      }
      
      fclose(xyz);
      
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

        if(mvlist != nullptr)
            delete mvlist;
    }

}//extern C
