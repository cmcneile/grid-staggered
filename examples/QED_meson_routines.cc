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


/**
   Remove the trace from a colour matrix.
 **/
  void make_traceless(LatticeColourMatrix & B) ;



void QED_mesons(GridCartesian & Grid,   
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


///////////////////////////////////////////////////////////////
//		compute quark propagator
///////////////////////////////////////////////////////////////



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
  std::vector<TComplex> pion_corr(nt)  ;
  std::vector<TComplex> rho_corr(nt)  ;
  std::vector<TComplex> rho_corr_av(nt)  ;

  std::vector<TComplex> a1_corr(nt)  ;
  std::vector<TComplex> a1_corr_av(nt)  ;

  // contract the quark propagators
  LatticeComplex  c(&Grid)  ;

   c = trace(Qprop * adj(Qprop)) ; 


  //  The correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  cout << "Tp = " << Tp  << endl; 
  sliceSum(c, pion_corr, Tp);

  // output the correlators
  cout << "Pseuodscalar Goldstone pion \n" ;
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      double ttt = real(pion_corr[tt]) ;
      cout << "PION " << tt << " "  <<  ttt  << endl ;
    }






}

