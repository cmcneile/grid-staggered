
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



void compute_onemp_localrho_hybrid(LatticeGaugeField & Umu, GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
				   int nt, int Tp)  ;




void compute_onemp_hybrid_BLOCK(LatticeGaugeField & Umu, GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
				BlockConjugateGradient<FermionField> & BCG,
			  ImprovedStaggeredFermionR & Ds, 
				int nt, int Tp) ;



void QED_mesons(GridCartesian & Grid,   
			  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> & HermOp,
			  ConjugateGradient<FermionField> & CG , 
			  ImprovedStaggeredFermionR & Ds, 
		int nt, int Tp) ;
