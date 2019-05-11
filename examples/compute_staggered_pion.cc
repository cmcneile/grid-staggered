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


//////////////////////////////////////////////
// Load Fermion vector into propagator 
//////////////////////////////////////////////

void FermToProp_s(LatticeStaggeredPropagator & Qprop, LatticeStaggeredFermion & psi , const int c)
{
  const int nc = 3 ;  // later use the dimension 

  for(int i = 0; i < nc ; ++i)
    {
      pokeColour(Qprop, peekColour(psi, i), i, c);
    }

}

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
  typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermionR::ComplexField ComplexField; 
  typename ImprovedStaggeredFermionR::ImplParams params; 

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

   //  create a hot su3 configuration
  //  SU3::HotConfiguration(pRNG,Umu);
  SU3::ColdConfiguration(Umu); // Umu = 1  

  int t_dir = 3;
  int nt =latt_size[t_dir];

  anti_peroidic( Umu , nt ) ;

  const int g_trans =  0;
  if( g_trans == 1)
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


  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop(&Grid)  ;

  LatticeStaggeredFermion D_out(&Grid) ;
  LatticeStaggeredFermion res(&Grid) ;
  LatticeStaggeredFermion tmp(&Grid) ;

  Qprop = zero ;

  cout << "Inversion for mass "   << mass << endl ; 

  // Compute the staggered quark propator
  // Solve for x 
  // M^dagger * M * x = M^dagger * src

  for(int ic = 0 ; ic < 3 ; ++ic)
    {
      cout << "---------------------------------------------" << endl ; 
      cout << "Inversion for colour " << ic << endl ; 
      cout << "---------------------------------------------" << endl ; 

      // create point source

      // ideas from tests/core/Test_staggered5Dvec.cc
      std::vector<int> site({0,0,0,0});
      ColourVector cv = zero;
      cv()()(ic)=1.0;
      local_src = zero;
      pokeSite(cv,local_src,site);

      // apply Mdagg
      Ds.Mdag(local_src, out) ;
      local_src = out ;

      // invert 
       out = zero ;  // intial guess

      CG(HermOp,local_src,out);

      // add solution to propagator structure
      FermToProp_s(Qprop, out , ic  ) ; 

      // compute the residual
       Ds.M(out, D_out) ;
       Ds.Mdag(D_out, tmp) ;
       D_out = tmp ; 

       res = D_out - local_src ;
       RealD nrm = norm2(res); 
       double xxx = (double) nrm*1.0 ;
       cout << "Residual = " <<  std::sqrt(xxx)  << endl ; 

    }

  // 
  //  -----  Use the quark propagator to compute the pion correlator
  //

  // pion correlator
  std::vector<TComplex> corr(nt)  ;

  // contract the quark propagators
  LatticeComplex  c(&Grid)  ;

   c = trace(Qprop * adj(Qprop)) ; 

  //  The correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  cout << "Tp = " << Tp  << endl; 
  sliceSum(c, corr, Tp);

  // output the correlators
  cout << "Pseuodscalar Goldstone pion \n" ;
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      double ttt = real(corr[tt]) ;
      cout << tt << " "  <<  ttt  << endl ;
    }

  // End of the Grid
  Grid_finalize();
}
