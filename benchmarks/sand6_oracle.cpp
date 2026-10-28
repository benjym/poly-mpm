#include <iostream>
#include <iomanip>
#include "simu/FormBuilder.hh"
#include "geo/MeshImpl.hh"
#include <bogus/Core/Block.impl.hpp>
int main() {
 using namespace d6;
 FormMat<6,3>::Type B; B.setRows(8); B.setCols(8);
 FormMat<3,3>::Type A; A.setRows(8); A.setCols(8);
 for(int i=0;i<8;++i) for(int j=0;j<8;++j){B.insertBack(i,j).setZero(); A.insertBack(i,j).setZero();}
 B.finalize(); A.finalize();
 int count; std::cin>>count;
 for(int p=0;p<count;++p){
  double x,y,z,v;std::cin>>x>>y>>z>>v;
  MeshType::Interpolation itp;MeshType::Derivatives grad;
  for(int i=0;i<8;++i){
   double xx=(i&1)?x:1-x, yy=(i&2)?y:1-y, zz=(i&4)?z:1-z;
   itp.nodes[i]=i;itp.coeffs[i]=xx*yy*zz;
   grad(i,0)=((i&1)?1:-1)*yy*zz;
   grad(i,1)=((i&2)?1:-1)*xx*zz;
   grad(i,2)=((i&4)?1:-1)*xx*yy;
  }
  std::vector<Index> ids={0,1,2,3,4,5,6,7};
  FormBuilder::addTauDu(B,v,itp,grad,ids,ids);
  FormBuilder::addDuDv(A,v,itp,grad,ids,ids);
 }
 std::cout<<std::setprecision(17);
 for(int i=0;i<8;++i)for(int a=0;a<6;++a){
  for(int j=0;j<8;++j)for(int b=0;b<3;++b)std::cout<<B.block(i,j)(a,b)<<' ';
  std::cout<<'\n';
 }
 for(int i=0;i<8;++i)for(int a=0;a<3;++a){
  for(int j=0;j<8;++j)for(int b=0;b<3;++b)std::cout<<A.block(i,j)(a,b)<<' ';
  std::cout<<'\n';
 }
}
