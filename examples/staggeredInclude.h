
  typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermionR::ComplexField ComplexField; 
//  typename ImprovedStaggeredFermionR::ImplParams params; 


//
//  FUNCTION protoptypes
//





void compute_local_mesons(GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
			  int nt, int Tp) ;



void compute_onelink_rho(LatticeGaugeField & Umu, GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
			 int nt, int Tp) ;


void compute_onemp_hybrid(LatticeGaugeField & Umu, GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
			  int nt, int Tp) ;



