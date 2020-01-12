    /*************************************************************************************


    Perhaps a better starting point:
    typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 


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

/*template<class d>
struct scal {
  d internal;
};

Gamma::Algebra Gmu [] = {
  Gamma::Algebra::GammaX,
  Gamma::Algebra::GammaY,
  Gamma::Algebra::GammaZ,
  Gamma::Algebra::GammaT
  };
*/



//////////////////////////////////////////////
// Load Fermion into propagator 
//////////////////////////////////////////////

void FermToProp_s(LatticeStaggeredPropagator &Qprop, LatticeStaggeredFermion &psi , const int c)
{
  const int nc = 3 ;  // later use the dimension 

  for(int i = 0; i < nc ; ++i)
    {
      pokeColour(Qprop, peekColour(psi, i), i, c);
    }

}

//////////////////////////////////////////////////
//  Apply anti-periodic boundary conditions in time
////////////////////////////////////////////////////

void anti_periodic( LatticeGaugeField & Umu , int nt)
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

  cout << "Anti-periodic boundary conditions in time applied" << endl ; 
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

///////////////////////////////////////////////////////////////
//		Staggered Phases
///////////////////////////////////////////////////////////////

  LatticeComplex phases[3] = {&Grid, &Grid, &Grid}; 
  LatticeComplex one(&Grid), minusOne(&Grid); one = 1; minusOne = -1;
  LatticeInteger coor(&Grid);
  for(int i=0; i<3; i++) {
    LatticeCoordinate(coor,i);	// fills coor with value of coord in i dir.
    phases[i] = where((mod(coor,2)==1), minusOne, one);
  }

//////////////////////////////////////////////////////////////

  LatticeGaugeField Umu(&Grid); 
  //SU3::HotConfiguration(pRNG,Umu); 
  SU3::ColdConfiguration(Umu); // Umu = 1

  int t_dir = 3;
  int nt = latt_size[t_dir];

  anti_periodic(Umu, nt);

  const int g_trans = 1;
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

  // This uses the Milc conventions. See gridStaggInvert.cc in the MILC code.
  // Naive staggered action
  RealD u0=1.0;
  RealD c1= 2.0 ;
  RealD c2= 0.0 ;

  // Naik coefficients (See reference in Grid library)
  // RealD c1=9.0/8.0; RealD c2=-1.0/24.0;

  RealD mass=0.1;

  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,2.0*mass,c1,c2,u0);

  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);
  ConjugateGradient<FermionField> CG(1.0e-8,10000);
  //LinearOperatorBase<FermionField> LinOp(Ds);

  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop(&Grid)  ;

  Qprop = zero ;

  cout << "Inversion for mass "   << mass << endl ; 

  // Compute the staggered quark propator
  for(int ic = 0 ; ic < 3 ; ++ic)
    {
      cout << "---------------------------------------------" << endl ; 
      cout << "Inversion for colour " << ic << endl ; 
      cout << "---------------------------------------------" << endl ; 

      // create point source
      // tests/core/Test_staggered5Dvec.cc
      
      std::vector<int> site({0,0,0,0});
      ColourVector cv = zero;
      cv()()(ic)=1.0;  
      local_src = zero;
      pokeSite(cv,local_src,site);

      Ds.Mdag(local_src, out) ; // apply Mdagger
      local_src = out;

      // invert 
      out = zero ;  // intial guess

      CG(HermOp,local_src,out);

      // add solution to propagator structure
      FermToProp_s(Qprop, out , ic  ) ; 
    }

  //  -----  Use the quark propagator to compute the pion correlator

  // rho correlator
  std::vector<vector<TComplex>> corr(3, vector<TComplex>(nt)) ;

  // contract the quark propagators
  LatticeComplex  c[3] = {&Grid, &Grid, &Grid};

  for(int j=0; j<3; j++) {
    c[j] = trace(Qprop * adj(Qprop)) ; 
    c[j] = c[j] * phases[j];	// goto line 104 to see how phases are generated.
  }

  //  this correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  cout << "Tp = " << Tp  << "\n"; 
  for(int k=0; k<3; k++) {
    sliceSum(c[k], corr[k], Tp);
  }

  // output the correlators
  cout << "Vector rho meson \n\n";

  for(int m=0; m<3; m++) {
    cout << "\nCorrelator in spin component " << m << endl; 
    for(int tt = 0 ; tt < nt ; ++tt) {
      
        double ttt = real(corr[m][tt]) ;
        cout << tt << " "  <<  ttt  << endl ;
    }
  }

Grid_finalize();

}
