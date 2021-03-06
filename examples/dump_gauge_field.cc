    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_nersc_io.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);

  std::cout <<GridLogMessage<< " main "<<std::endl;

  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> latt_size  ({8,8,8,8});

  int orthodir=3;
  int orthosz =latt_size[orthodir];
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);

  GridParallelRNG   pRNGa(&Fine);
  GridParallelRNG   pRNGb(&Fine);
  GridSerialRNG     sRNGa;
  //  GridSerialRNG     sRNGb;

  std::cout <<GridLogMessage<< " seeding... "<<std::endl;
  pRNGa.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  sRNGa.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  std::cout <<GridLogMessage<< " ...done "<<std::endl;

  LatticeGaugeField Umu(&Fine);
  LatticeGaugeField Umu_diff(&Fine);
  LatticeGaugeField Umu_saved(&Fine);

  std::vector<LatticeColourMatrix> U(4,&Fine);
  
  SU3::HotConfiguration(pRNGa,Umu);

  RealD WFlow_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  RealD WFlow_TC   = WilsonLoops<PeriodicGimplR>::TopologicalCharge(Umu);

  std::cout << GridLogMessage << "Plaquette          " << WFlow_plaq << std::endl;
  std::cout << GridLogMessage << "TopologicalCharge  " << WFlow_TC   << std::endl;

  FieldMetaData header;

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Writing out  ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::string file("./ckpoint_ildg.4000");
  IldgWriter _IldgWriter(Fine.IsBoss());
  _IldgWriter.open(file);
  _IldgWriter.writeConfiguration(Umu,4000,std::string("dummy_ildg_LFN"),std::string("dummy_config"));
  _IldgWriter.close();


  Grid_finalize();
#endif
}
