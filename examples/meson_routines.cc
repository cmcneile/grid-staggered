#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

#include"staggeredInclude.h"

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


void symm_shift(LatticeGaugeField &Umu, LatticeStaggeredFermion &q, int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);

  q =  0.5 * U * Cshift(q, dir, shift) + adj(Cshift(U, dir, -shift))*Cshift(q, dir, -shift) ;

}


void compute_local_mesons(GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
			  int nt, int Tp)
{


  // This workspace should be not specific to this routine
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop(&Grid)  ;

  LatticeStaggeredFermion D_out(&Grid) ;
  LatticeStaggeredFermion res(&Grid) ;
  LatticeStaggeredFermion tmp(&Grid) ;



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




}



void compute_onelink_rho(LatticeGaugeField & Umu, GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
			  int nt, int Tp)
{



///////////////////////////////////////////////////////////////
//		Staggered Phases
///////////////////////////////////////////////////////////////

  LatticeComplex phases(&Grid); 
  LatticeComplex one(&Grid), minusOne(&Grid); one = 1; minusOne = -1;
  LatticeInteger coor[4] = {&Grid, &Grid, &Grid, &Grid}; LatticeInteger r(&Grid);

  for(int m=0; m<4; m++) {
    LatticeCoordinate(coor[m], m);  
  }

  r = coor[0] + coor[1] + coor[2] + coor[3];

    phases = where((mod(r,2)== (Integer) 0), minusOne, one);
  

//////////////////////////////////////////////////////////////



  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop[2] = {&Grid, &Grid}  ;

  Qprop[1] = Qprop[0] = zero ;

  // Compute the staggered quark propagators
  for(int k=0; k<2; k++) {

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
      
        // shift the source
        if(k) symm_shift(Umu, local_src, 0, 1); // do the unshifted qprop first

        Ds.Mdag(local_src, out) ; // apply Mdagger
        local_src = out;

        // invert 
        out = zero ;  // intial guess
        CG(HermOp,local_src,out);

        // apply the shifted operator to the sink
	if(k) symm_shift(Umu, out, 0, 1); 

        // add solution to propagator structure
        FermToProp_s(Qprop[k], out , ic  ) ; 
      }
  }
  //  -----  Use the quark propagator to compute the correlator

  // 1-link rho correlator
  std::vector<vector<TComplex>> corr(3, vector<TComplex>(nt)) ;

  // contract the quark propagators
  LatticeComplex  c[3] = {&Grid, &Grid, &Grid};

  for(int j=0; j<3; j++) {
    c[j] = trace(adj(Qprop[0]) * Qprop[1]) ; 
    c[j] = c[j] * phases;	
  }

  //  this correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  cout << "\nTp = " << Tp  << "\n"; 
  for(int k=0; k<3; k++) {
    sliceSum(c[k], corr[k], Tp);
  }

  // output the correlators
  cout << "\n\nSHIFTED rho meson \n\n";

  for(int m=0; m<3; m++) {
    cout << "\nCorrelator in spin component " << m << endl; 
    for(int tt = 0 ; tt < nt ; ++tt) {
      
        double ttt = real(corr[m][tt]) ;
        cout << tt << " "  <<  ttt  << endl ;
    }
  }


}
