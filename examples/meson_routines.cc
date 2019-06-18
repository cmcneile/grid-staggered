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


void symm_shift_b(LatticeGaugeField &Umu, LatticeStaggeredFermion &q, int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);

  q =  0.5 * Cshift(q, dir, shift) +  0.5 *Cshift(q, dir, -shift) ;

}


void symm_shift_a(LatticeGaugeField &Umu, LatticeStaggeredFermion &q, int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);

  q =  0.5 * U * Cshift(q, dir, shift) + 0.5 * adj(Cshift(U, dir, -shift))*Cshift(q, dir, -shift) ;

}


LatticeStaggeredFermion symm_shift_n(LatticeGaugeField &Umu, LatticeStaggeredFermion q, int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);
  LatticeStaggeredFermion tmp(grid);
  tmp = 0.5 * (U * Cshift(q, dir, shift) + adj(Cshift(U, dir, -shift))*Cshift(q, dir, -shift)) ;
  return tmp;
}



void symm_shift(LatticeGaugeField &Umu, LatticeStaggeredFermion &q, int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);

  q =  0.5 * U * Cshift(q, dir, shift) + 0.5 *adj(Cshift(U, dir, -shift))*Cshift(q, dir, -shift) ;

}

void symm_shift_c(LatticeGaugeField &Umu, LatticeStaggeredFermion &q, int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);

  q =  U * Cshift(q, dir, shift) ;

}

/**
   1-+ hybrid staggered operator.

 **/

LatticeStaggeredFermion hybrid_op(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs,
                                  int dir)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bi(grid); LatticeColourMatrix Bj(grid);
  int i, j;
  const int shift = 1 ;


  if(dir==0) { i=1; j=2;} 
  if(dir==1) {i=2; j=0;} 
  if(dir==2) {i=0; j=1;}

  // will need another conditional to ensure dir < i, dir < j.
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, j, dir); 
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bj, Umu, i, dir); 

  LatticeStaggeredFermion tmp(grid);

  tmp = signs[i]*symm_shift_n(Umu, Bj*q , i, shift) + Bj*signs[i]*symm_shift_n(Umu, q, i, shift)
    - signs[j]*symm_shift_n(Umu, Bi*q, j, shift)  - Bi*signs[j]*symm_shift_n(Umu, q, j, shift); 
  // generalise this and/or make it more compact 
  return tmp;
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


/*****

  Flavour singlet one link rho operator.

 ***/

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

  Lattice<iScalar<vInteger> > tL(&Grid); LatticeCoordinate(tL,3);


  const int XUP = 0 ;
  const int YUP = 1 ;
  const int ZUP = 2 ;
  const int TUP = 3 ;


  //      r =  coor[ZUP] + coor[TUP];
  //     r =  coor[ZUP] ; 
          r =  coor[TUP] ;   // no
  //     r =  coor[XUP] + coor[YUP] + coor[3] ;  // no
  //    r =  coor[XUP] + coor[YUP] ;  // no

  //       r =  coor[XUP] + coor[3] ;   // no
  //  r =  coor[XUP] + coor[YUP] +coor[ZUP] + coor[3] ;   //

  //          r = coor[0] + coor[1] + coor[2] + coor[3];
  //           r = coor[XUP] + coor[YUP] + coor[ZUP] ;

  //  phases = (-1)**(x+y+z+t)
  phases = where((mod(r,2)== (Integer) 0), one, minusOne);



  //    cout << phases << "\n" ;

//////////////////////////////////////////////////////////////



  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop[2] = {&Grid, &Grid}  ;

  Qprop[1] = Qprop[0] = zero ;
  const int shift_dir = 2 ;  // x-direction  debug z-direction

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
        if(k) symm_shift(Umu, local_src, shift_dir , 1); // do the unshifted qprop first

        Ds.Mdag(local_src, out) ; // apply Mdagger
        local_src = out;

        // invert 
        out = zero ;  // intial guess
        CG(HermOp,local_src,out);

        // apply the shifted operator to the sink
	if(k) symm_shift(Umu, out, shift_dir, 1); 

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
    //    c[j] = trace(adj(Qprop[1]) * Qprop[1]) ; 
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

  string dir_name[4] = {"X", "Y", "Z", "T"}; 

  { int m = shift_dir ;
    cout << "\nCorrelator in spin component " << dir_name[m]  << endl; 
    for(int tt = 0 ; tt < nt ; ++tt) {
      
        double ttt = real(corr[m][tt]) ;
        double ttt_img = imag(corr[m][tt]) ;
        cout << tt << " "  <<  ttt << "  "   << ttt_img  << endl ;
    }
  }


}






/*****

  1-+ hybrid operator

 ***/

void compute_onemp_hybrid(LatticeGaugeField & Umu, GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
			  int nt, int Tp)
{

///////////////////////////////////////////////////////////////
//		Staggered Sign Functions
///////////////////////////////////////////////////////////////

  LatticeInteger coor[4] = {&Grid, &Grid, &Grid, &Grid}; 
  LatticeInteger n[5] = {&Grid, &Grid, &Grid, &Grid, &Grid};
  LatticeComplex signs[5] =  {&Grid, &Grid, &Grid, &Grid, &Grid}; 
  LatticeComplex One(&Grid), minusOne(&Grid); One = 1; minusOne = -1;

  for(int m=0; m<4; m++) {
    LatticeCoordinate(coor[m], m);  
  }

  n[0] = 1; 
  n[1] = coor[0]; 
  n[2] = n[1] + coor[1]; 
  n[3] = n[2] + coor[2];
  n[4] = n[3] + coor[3];
  
  for(int i=0; i<5; i++) {
    signs[i] = where((mod(n[i],2)== (Integer) 1), minusOne, One) ;
  } 

// Decided to split the phases up a bit, here are the standard staggered sign functions, corresponding
// to the gamma matrices. the phase from the inversion of M is signs[4].

  enum dirs {XUP=0, YUP=1, ZUP=2, TUP=3};
  dirs shift_dir = XUP ; // index of hybrid operator

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
        std::vector<int> site({0,0,0,0});
        ColourVector cv = zero;
        cv()()(ic)=1.0;  
        local_src = zero;
        pokeSite(cv,local_src,site);
      
        // shift the source
        if(k) local_src = hybrid_op(Umu, local_src, signs, shift_dir); // do the unshifted qprop first

        Ds.Mdag(local_src, out) ; // apply Mdagger
        local_src = out;

        // invert 
        out = zero ;  // intial guess
        CG(HermOp,local_src,out);

        // apply the shifted operator to the sink
	if(k) out = hybrid_op(Umu, out, signs, shift_dir);

        // add solution to propagator structure
        FermToProp_s(Qprop[k], out , ic  ) ; 
      }
  }
  //  -----  Use the quark propagator to compute the correlator

  // correlator
   std::vector<TComplex> corr(nt) ;

  // contract the quark propagators
  LatticeComplex  c(&Grid);

  for(int j=0; j<3; j++) {
    c = trace(adj(Qprop[0]) * Qprop[1]) ; 
    c = c * signs[4];	// phase from inversion of M (aka epsilon)
  }

  //  this correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  //  cout << "\nTp = " << Tp  << "\n"; 
  sliceSum(c, corr, Tp);

  string dir_name[4] = {"X", "Y", "Z", "T"}; 

  // output the correlators
  cout << "\n\nSHIFTED Hybrid 1-+ dir= " << dir_name[shift_dir]  << "\n";

    for(int tt = 0 ; tt < nt ; ++tt) {
           double ttt = real(corr[tt]) ;
        double ttt_img = imag(corr[tt]) ;
        cout << tt << " "  <<  ttt << "  "   << ttt_img  << endl ;
    }
  


}




