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


void PropToFerm_s( LatticeStaggeredFermion & psi ,  LatticeStaggeredPropagator & Qprop, const int c) ;


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

  LatticeStaggeredFermion out_seq(&Grid) ;


  LatticeStaggeredFermion D_out(&Grid) ;
  LatticeStaggeredFermion res(&Grid) ;
  LatticeStaggeredFermion tmp(&Grid) ;


  LatticeStaggeredPropagator Qprop(&Grid)  ;

  LatticeStaggeredPropagator Seq1_Qprop(&Grid)  ;
  LatticeStaggeredPropagator Seq2_Qprop(&Grid)  ;

  //  QED stuff  
  //  tests//core/Test_qed.ccb
  PhotonR photon(&Grid, PhotonR::Gauge::coulomb, PhotonR::ZmScheme::qedL);



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
       cout << "Residual (first inversion) = " <<  std::sqrt(xxx)  << endl ; 

    }

  //
  //  compute the first sequential propagator
  //



  for(int ic = 0 ; ic < 3 ; ++ic)
    {
      cout << "---------------------------------------------" << endl ; 
      cout << "Sequential Inversion for colour " << ic << endl ; 
      cout << "---------------------------------------------" << endl ; 

      PropToFerm_s(local_src, Qprop, ic  ) ; 

      // apply Mdagg
      Ds.Mdag(local_src, out) ;
      local_src = out ;

      // invert 
       out = zero ;  // intial guess

      CG(HermOp,local_src,out);

      // add solution to propagator structure
      FermToProp_s(Seq1_Qprop, out , ic  ) ; 

      // compute the residual
       Ds.M(out, D_out) ;
       Ds.Mdag(D_out, tmp) ;
       D_out = tmp ; 

       res = D_out - local_src ;
       RealD nrm = norm2(res); 
       double xxx = (double) nrm*1.0 ;
       cout << "Residual (first inversion) = " <<  std::sqrt(xxx)  << endl ; 


       //
       //  compute the first sequential propagator
       //

      CG(HermOp,out,out_seq);

      // add solution to propagator structure
      FermToProp_s(Seq1_Qprop, out_seq , ic  ) ; 


      // compute the residual
       Ds.M(out_seq, D_out) ;
       Ds.Mdag(D_out, tmp) ;
       D_out = tmp ; 

       res = D_out - local_src ;
       nrm = norm2(res); 
       xxx = (double) nrm*1.0 ;
       cout << "Residual (sequential inversion) = " <<  std::sqrt(xxx)  << endl ; 
       // probaly need norm of source.


    }






  // 
  //  -----  Use the quark propagator to compute the pion correlator
  //

  // pion correlator
  std::vector<TComplex> pion_corr(nt)  ;
  std::vector<TComplex> pion_corr_QED(nt)  ;

  //  std::vector<TComplex> rho_corr(nt)  ;
  //std::vector<TComplex> rho_corr_av(nt)  ;

  //std::vector<TComplex> a1_corr(nt)  ;
  //std::vector<TComplex> a1_corr_av(nt)  ;

  // contract the quark propagators
  LatticeComplex  c(&Grid)  ;


   c = trace(Qprop * adj(Qprop)) ; 
  //  The correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  sliceSum(c, pion_corr, Tp);

  // output the correlators
  cout << "Pseuodscalar Goldstone pion \n" ;
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      double ttt = real(pion_corr[tt]) ;
      cout << "PION " << tt << " "  <<  ttt  << endl ;
    }

  /**
   **     ***  QED  *****
   **/


   c = trace(Seq1_Qprop * adj(Seq1_Qprop)) ; 
  //  The correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  sliceSum(c, pion_corr_QED, Tp);

  // output the correlators
  cout << "Pseuodscalar (QED correction) Goldstone pion \n" ;
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      double ttt = real(pion_corr_QED[tt]) ;
      cout << "PION-QED0 " << tt << " "  <<  ttt  << endl ;
    }


}

