    /*************************************************************************************

    This is a test code for computing staggered correlators using the
    Grid library. At the moment this is indended as a test code, to make
    sure we the phases correct.

    The Grid library doesn't currently have routines to fatten the links,
    because Grid routines are intended to be called by application code.

    This uses naive staggered fermions.

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_unprec.cc
    Starting point:    Test_staggered_cg_unprec.cc 


    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

#include"staggeredInclude.h"


//////////////////////////////////////////////
// Load Fermion vector into propagator 
//////////////////////////////////////////////

void FermToProp_s(LatticeStaggeredPropagator & Qprop, LatticeStaggeredFermion & psi , const int c) ;


//////////////////////////////////////////////////
//  Apply anti-peroidic boundary conditions in time
////////////////////////////////////////////////////

void anti_peroidic( LatticeGaugeField & Umu , int nt)
{
  int mu = 3 ;  // time directiom
  GridBase *grid = Umu._grid;
  // conformable(grid,g._grid);
  LatticeColourMatrix U(grid);

  U= PeekIndex<LorentzIndex>(Umu,mu);

  // code hacked from Grid/Grid/qcd/action/fermion/FermionOperatorImpl.h`
  Lattice<iScalar<vInteger> > t(grid); LatticeCoordinate(t,3);
  LatticeComplex phases(grid); phases=1.0;

  --nt ;
  phases = where( t ==(Integer)nt, phases,-phases);
  U = phases * U ;

  PokeIndex<LorentzIndex>(Umu,U,mu);

  cout << "Anti-peroidic boundary conditions in time applied" << endl ; 
}



int main (int argc, char ** argv)
{

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  for(int mu=0;mu<Nd;mu++) 
    {
      cout << "Lattice[ " << mu << " ]= " << latt_size[mu] << endl ;
    }

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  FermionField result(&Grid); result=zero;
  
   LatticeGaugeField Umu(&Grid); 

   // https://stackoverflow.com/questions/12183008/how-to-use-enums-in-c
   enum cfgStart  { Cold, Hot, Load } ;
   cfgStart cfg_start = Load ;


   if ( cfg_start  == Hot ) {
     std::cout << "Hot configuration created\n";
     SU3::HotConfiguration(pRNG,Umu);
   }
   else if ( cfg_start  == Cold ) {
     std::cout << "Cold configuration created\n";
    SU3::ColdConfiguration(Umu); // Umu = 1  
   }
   else if ( cfg_start  == Load ) {
     std::cout << "Loading the configuration from disk\n";
     // read configuration in (Test_ildg_io.cc)
     FieldMetaData header;
     
     std::cout <<GridLogMessage<<"**************************************"<<std::endl;
     std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
     std::cout <<GridLogMessage<<"**************************************"<<std::endl;
     emptyUserRecord record;
     std::string file("./ckpoint_ildg.4000");
     
     IldgReader _IldgReader;
     _IldgReader.open(file);
     _IldgReader.readConfiguration(Umu,header);
     _IldgReader.close();
   }
   else
     {
       cout << "ERROR cfg_start = " << cfg_start << "\n" ;
       exit(1) ;
     }


  int t_dir = 3;
  int nt =latt_size[t_dir];

  anti_peroidic( Umu , nt ) ;

  enum gtrans_options  { do_trans , no_trans } ;
  gtrans_options g_trans  = no_trans ;

  if( g_trans == do_trans )
    {
      LatticeColourMatrix   g(&Grid); // Gauge xform
      SU3::RandomGaugeTransform(pRNG,Umu,g); // Unit gauge
      cout << "Random Gauge Transform applied "  << endl ; 
    }
  else
    {
   cout << "NO Gauge Transform applied "  << endl ; 
    }

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  
  
  RealD mass=0.1;

  // This uses the Milc conventions. See gridStaggInvert.cc in the MILC code.
  // Naive staggered action
  RealD u0=1.0;
  RealD c1= 2.0 ;
  RealD c2= 0.0 ;

  // Naik coefficients (See reference in Grid library)
  //  RealD c1=9.0/8.0;
  // RealD c2=-1.0/24.0;

  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,2.0*mass,c1,c2,u0);


  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);
  ConjugateGradient<FermionField> CG(1.0e-8,10000);


  // compute the meson spectrum
  compute_local_mesons(Grid, HermOp, CG, Ds, nt, Tp );  

  /* ----   start of peturbative QED analysis  ---- **/

  std::cout << "Start of QED analysis\n";

  QED_mesons(Grid, HermOp, CG, Ds, nt, Tp );  

  // End of the Grid
  Grid_finalize();
}

