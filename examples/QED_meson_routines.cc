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
//		Rho  Staggered Phases
///////////////////////////////////////////////////////////////

  LatticeComplex rho_phases[3] = {&Grid, &Grid, &Grid}; 
  LatticeComplex one(&Grid), minusOne(&Grid); 
  one = 1; 
  minusOne = -1;

  LatticeInteger coor(&Grid);
  for(int i=0; i<3; i++) {
    LatticeCoordinate(coor,i);	// fills coor with value of coord in i dir.
    rho_phases[i] = where((mod(coor,2)==(Integer) 1), minusOne, one);
  }

///////////////////////////////////////////////////////////////
//		A1  Staggered Phases
///////////////////////////////////////////////////////////////

  LatticeComplex a1_phases[3] = {&Grid, &Grid, &Grid}; 
  LatticeInteger coor_x(&Grid);
  LatticeCoordinate(coor_x,0);	

  LatticeInteger coor_y(&Grid);
  LatticeCoordinate(coor_y,1);	

  LatticeInteger coor_z(&Grid);
  LatticeCoordinate(coor_z,2);	

  LatticeInteger coor_t(&Grid);
  LatticeCoordinate(coor_t,3);	

  coor = coor_y+coor_z+coor_t ;
  a1_phases[0] = where((mod(coor,2)==(Integer) 1), minusOne, one);

  coor = coor_x+coor_z+coor_t ;
  a1_phases[1] = where((mod(coor,2)==(Integer) 1), minusOne, one);

  coor = coor_x+coor_y+coor_t ;
  a1_phases[2] = where((mod(coor,2)==(Integer) 1), minusOne, one);


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


  /**
     Rho meson
   **/

  LatticeComplex  c_rho(&Grid) ;
  cout << "Vector meson \n" ;
  for(int j=0; j<3; j++) 
    {
      c_rho = c * rho_phases[j] ;	
      sliceSum(c_rho  , rho_corr, Tp);
      for(int tt = 0 ; tt < nt ; ++tt)
	{
	  double ttt = real(rho_corr[tt]) ;
	  cout << "RHO[" << j <<  "] " << tt << " "  <<  ttt  << endl ;
	  rho_corr_av[tt]  += rho_corr[tt] ;
	}
    }
      for(int tt = 0 ; tt < nt ; ++tt)
	{
	  double ttt = real(rho_corr_av[tt]) ;
	  cout << "RHO[av] " << tt << " "  <<  ttt  << endl ;
	}

      /***
	  a1 meson
       ***/

  cout << "A1 meson \n" ;
  LatticeComplex  c_a1(&Grid) ;
  for(int j=0; j<3; j++) 
    {
      c_a1= c * a1_phases[j] ;	
      sliceSum(c_a1  , a1_corr, Tp);
      for(int tt = 0 ; tt < nt ; ++tt)
	{
	  double ttt = real(a1_corr[tt]) ;
	  cout << "A1[" << j <<  "] " << tt << " "  <<  ttt  << endl ;
	  a1_corr_av[tt]  += a1_corr[tt] ;
	}
    }
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      double ttt = real(a1_corr_av[tt]) ;
      cout << "A1[av] " << tt << " "  <<  ttt  << endl ;
    }





}

