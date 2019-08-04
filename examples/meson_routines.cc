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



LatticeStaggeredFermion symm_shift_n(LatticeGaugeField &Umu, LatticeStaggeredFermion q, int dir)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  const int shift = 1 ;


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

/**
   Remove the trace from a colour matrix.
 **/
void make_traceless(LatticeColourMatrix & B)
{
  GridBase *grid = B._grid;
 LatticeColourMatrix  con(grid);
 con = 1/3.0 ;

  LatticeComplex  tt(grid);
  tt = trace(B) ; 
 
  B -= tt * con ;

}

/**
   1-+ hybrid staggered operator using the (gamma_i cross 1) rho operator.

 **/

LatticeStaggeredFermion hybrid_op(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs,
                                  int dir)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bi(grid); LatticeColourMatrix Bj(grid);
  int i, j;
  const int shift = 1 ;


  if(dir==0) 
    { i=1; j=2;} 
  else if(dir==1) 
    {i=2; j=0;} 
  else if(dir==2) 
    {i=0; j=1;}
  else
    {
      cout << "hybrid_op:: Error dir = " << dir << " out of range\n" ;
      exit(0) ;
    }

  // compute the field strength tensor
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, j, dir); 
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bj, Umu, i, dir); 

   make_traceless(Bi) ;
   make_traceless(Bj) ;

  LatticeStaggeredFermion tmp(grid);

  double milc_nrm = 8.0 ; 
  tmp = signs[i]*symm_shift_n(Umu, Bj*q ,i) + signs[i]*Bj*symm_shift_n(Umu, q, i)
    -   signs[j]*symm_shift_n(Umu, Bi*q, j) - Bi*signs[j]*symm_shift_n(Umu, q, j); 

  // Field strength in Grid devides by a factor of 8, but the MILC code does not
  tmp *= milc_nrm ; 

  return tmp;
}



/**
   1-+ hybrid staggered operator using the (gamma_i cross gamma_i ) rho operator.

 **/

LatticeStaggeredFermion hybrid_localrho_op(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs,
                                  int dir)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bi(grid); LatticeColourMatrix Bj(grid);
  int i, j;
  const int shift = 1 ;


  if(dir==0) 
    { i=1; j=2;} 
  else if(dir==1) 
    {i=2; j=0;} 
  else if(dir==2) 
    {i=0; j=1;}
  else
    {
      cout << "hybrid_op:: Error dir = " << dir << " out of range\n" ;
      exit(0) ;
    }

  // compute the field strength tensor
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, j, dir); 
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bj, Umu, i, dir); 

   make_traceless(Bi) ;
   make_traceless(Bj) ;

  LatticeStaggeredFermion tmp(grid);

  double milc_nrm = 8.0 ; 
  tmp = signs[i]*Bj*q  - signs[j]*Bi*q ;

  // Field strength in Grid devides by a factor of 8, but the MILC code does not
  tmp *= milc_nrm ; 

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
  for(int i=0; i<3; i++) {
    LatticeCoordinate(coor,i);	// fills coor with value of coord in i dir.
    a1_phases[i] = where((mod(coor,2)==(Integer) 1), minusOne, one);
  }


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
  std::vector<TComplex> a1_corr(nt)  ;

  // contract the quark propagators
  LatticeComplex  c(&Grid)  ;

   c = trace(Qprop * adj(Qprop)) ; 

#if 0
  // contract the quark propagators
  LatticeComplex  c_rho[3] = {&Grid, &Grid, &Grid};

  for(int j=0; j<3; j++) {
        c_rho[j] = c * rho_phases[j] ;	
  }
#endif

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
	}
    }

  cout << "A1 meson \n" ;
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      double ttt = real(a1_corr[tt]) ;
      cout << "A1 " << tt << " "  <<  ttt  << endl ;
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

  1-+ hybrid operator using the one link rho operator

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

  //  n[0] = 1; 
  n[0] = 1; 
  n[1] = coor[0]; 
  n[2] = n[1] + coor[1]; 
  n[3] = n[2] + coor[2];
  n[4] = n[3] + coor[3];
  
  for(int i=0; i<5; i++) {
    signs[i] = where((mod(n[i],2)== (Integer) 1), minusOne, One) ;
  } 

// The standard staggered sign functions, corresponding
// to the gamma matrices. the phase from the inversion of M is signs[4].

  enum dirs {XUP=0, YUP=1, ZUP=2, TUP=3};
  dirs shift_dir = ZUP ; // index of hybrid operator

//////////////////////////////////////////////////////////////


  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop[2] = {&Grid, &Grid}  ;

  Qprop[1] = Qprop[0] = zero ;

  // Compute the staggered quark propagators
  for(int k=0; k<2; k++) 
    {

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

  c = trace(adj(Qprop[0]) * Qprop[1]) ; 
  c = c * signs[4];	// phase from inversion of M (aka epsilon)
    // This should be in, but it gets almost MILC correlators.


  //  The correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  //  cout << "\nTp = " << Tp  << "\n"; 
  sliceSum(c, corr, Tp);

  string dir_name[4] = {"X", "Y", "Z", "T"}; 

  // output the correlators
  cout << "\n\nSHIFTED Hybrid 1-+ dir= " << dir_name[shift_dir]  << "\n";

  const double sss = -1 ;
  for(int tt = 0 ; tt < nt ; ++tt) 
    {
      double ttt      = sss*real(corr[tt]) ;
      double ttt_img  = sss*imag(corr[tt]) ;

      cout << "ONEMP " << tt << " "  <<  ttt << "  "   << ttt_img  << endl ;
    }
  


}







/*****

  1-+ hybrid operator using the one link rho operator

 ***/

void compute_onemp_hybrid_BLOCK(LatticeGaugeField & Umu, GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
				BlockConjugateGradient<FermionField> & BCG,
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

  //  n[0] = 1; 
  n[0] = 1; 
  n[1] = coor[0]; 
  n[2] = n[1] + coor[1]; 
  n[3] = n[2] + coor[2];
  n[4] = n[3] + coor[3];
  
  for(int i=0; i<5; i++) {
    signs[i] = where((mod(n[i],2)== (Integer) 1), minusOne, One) ;
  } 

// The standard staggered sign functions, corresponding
// to the gamma matrices. the phase from the inversion of M is signs[4].

  enum dirs {XUP=0, YUP=1, ZUP=2, TUP=3};
  dirs shift_dir = ZUP ; // index of hybrid operator

//////////////////////////////////////////////////////////////


  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop[2] = {&Grid, &Grid}  ;

  Qprop[1] = Qprop[0] = zero ;

  /** use two sources  ***/
  int no_vec = 2 ;
  std::vector<LatticeStaggeredFermion> v_local_src(no_vec, &Grid) ;
  std::vector<LatticeStaggeredFermion> v_out(no_vec, &Grid) ;

  // Compute the staggered quark propagators
    {

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
	  
	  // local source
	  Ds.Mdag(local_src, v_local_src[0]) ; // apply Mdagger

	  // apply the hybrid operator onto point source
	  local_src = hybrid_op(Umu, local_src, signs, shift_dir); // do the unshifted qprop first
	  Ds.Mdag(local_src,v_local_src[1] ) ; // apply Mdagger

	  // invert 
	  v_out[0] = zero ;  // intial guess
	  v_out[1] = zero ;  // intial guess

	   BCG(HermOp,v_local_src,v_out);

	  // apply the hybrid operator to the sink
	  v_out[1] = hybrid_op(Umu, v_out[1], signs, shift_dir);

	  // add solution to propagator structure
	  FermToProp_s(Qprop[0], v_out[0] , ic  ) ; 
	  FermToProp_s(Qprop[1], v_out[1] , ic  ) ; 
	}
    }

  //  -----  Use the quark propagator to compute the correlator

  // correlator
   std::vector<TComplex> corr(nt) ;

  // contract the quark propagators
  LatticeComplex  c(&Grid);

  c = trace(adj(Qprop[0]) * Qprop[1]) ; 
  c = c * signs[4];	// phase from inversion of M (aka epsilon)
    // This should be in, but it gets almost MILC correlators.


  //  The correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  //  cout << "\nTp = " << Tp  << "\n"; 
  sliceSum(c, corr, Tp);

  string dir_name[4] = {"X", "Y", "Z", "T"}; 

  // output the correlators
  cout << "\n\nSHIFTED Hybrid 1-+ dir= " << dir_name[shift_dir]  << "\n";

  const double sss = -1 ;
  for(int tt = 0 ; tt < nt ; ++tt) 
    {
      double ttt      = sss*real(corr[tt]) ;
      double ttt_img  = sss*imag(corr[tt]) ;

      cout << "ONEMP-BLOCK " << tt << " "  <<  ttt << "  "   << ttt_img  << endl ;
    }
  


}





/*****

  1-+ hybrid operator using the local rho operator

 ***/

void compute_onemp_localrho_hybrid(LatticeGaugeField & Umu, GridCartesian & Grid,   
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

  //  n[0] = 1; 
  n[0] = 1; 
  n[1] = coor[0]; 
  n[2] = n[1] + coor[1]; 
  n[3] = n[2] + coor[2];
  n[4] = n[3] + coor[3];
  
  for(int i=0; i<5; i++) {
    signs[i] = where((mod(n[i],2)== (Integer) 1), minusOne, One) ;
  } 

// The standard staggered sign functions, corresponding
// to the gamma matrices. the phase from the inversion of M is signs[4].

  enum dirs {XUP=0, YUP=1, ZUP=2, TUP=3};
  dirs shift_dir = ZUP ; // index of hybrid operator

//////////////////////////////////////////////////////////////


  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop[2] = {&Grid, &Grid}  ;

  Qprop[1] = Qprop[0] = zero ;

  // Compute the staggered quark propagators
  for(int k=0; k<2; k++) 
    {

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
	  if(k) local_src = hybrid_localrho_op(Umu, local_src, signs, shift_dir); // do the unshifted qprop first
	  
	  Ds.Mdag(local_src, out) ; // apply Mdagger
	  local_src = out;

	  // invert 
	  out = zero ;  // intial guess
	  CG(HermOp,local_src,out);

	  // apply the shifted operator to the sink
	  if(k) out = hybrid_localrho_op(Umu, out, signs, shift_dir);

	  // add solution to propagator structure
	  FermToProp_s(Qprop[k], out , ic  ) ; 
	}
    }

  //  -----  Use the quark propagator to compute the correlator

  // correlator
   std::vector<TComplex> corr(nt) ;

  // contract the quark propagators
  LatticeComplex  c(&Grid);

  c = trace(adj(Qprop[0]) * Qprop[1]) ; 
  c = c * signs[4];	// phase from inversion of M (aka epsilon)
    // This should be in, but it gets almost MILC correlators.


  //  The correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  //  cout << "\nTp = " << Tp  << "\n"; 
  sliceSum(c, corr, Tp);

  string dir_name[4] = {"X", "Y", "Z", "T"}; 

  // output the correlators
  cout << "\n\nSHIFTED Hybrid 1-+ dir= " << dir_name[shift_dir]  << "\n";

  const double sss = -1 ;
  for(int tt = 0 ; tt < nt ; ++tt) 
    {
      double ttt      = sss*real(corr[tt]) ;
      double ttt_img  = sss*imag(corr[tt]) ;

      cout << "ONEMP-LOCALRHO-Z " << tt << " "  <<  ttt << "  "   << ttt_img  << endl ;
    }
  


}




