#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <vector>
#include "bandmatrix.cpp"

using namespace std;

class Spline
{
public:
  Spline();
  void set_spline(vector<double>, vector<double>, vector<double>&, vector<double>&, vector<double>&);
  double seval(vector<double>, vector<double>, vector<double>, vector<double>, vector<double>,double);
  double seval_derivative(vector<double>, vector<double>, vector<double>, vector<double>, vector<double>,double);


};


Spline::Spline()
{

}

void Spline::set_spline(vector<double> x, vector<double> y, vector<double>& a, vector<double>& b, vector<double>& c)
{
  assert(x.size()==y.size());
  int n=x.size();
  band_matrix A(n,1,1);
  vector<double>  rhs(n);
      for(int i=1; i<n-1; i++) {
         A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
         A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
         A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
         rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      }
  // boundary conditions, zero curvature b[0]=b[n-1]=0
  A(0,0)=2.0;
  A(0,1)=0.0;
  rhs[0]=0.0;
  A(n-1,n-1)=2.0;
  A(n-1,n-2)=0.0;
  rhs[n-1]=0.0;

  // solve the equation system to obtain the parameters b[]
  b=A.lu_solve(rhs);


  // calculate parameters a[] and c[] based on b[]
  a.resize(n);
  c.resize(n);
  for(int i=0; i<n-1; i++) {
      a[i]=1.0/3.0*(b[i+1]-b[i])/(x[i+1]-x[i]);
      c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-1.0/3.0*(2.0*b[i]+b[i+1])*(x[i+1]-x[i]);
     }
  // for the right boundary we define
  // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
  double h=x[n-1]-x[n-2];
  // b[n-1] is determined by the boundary condition
  a[n-1]=0.0;
  c[n-1]=3.0*a[n-2]*h*h+2.0*b[n-2]*h+c[n-2];   // = f'_{n-2}(x_{n-1})
}

double Spline::seval(vector<double> x, vector<double> y, vector<double> a, vector<double> b, vector<double> c,double u)
{
  assert(x.size()==y.size());
  int n=x.size();
  std::vector<double>::const_iterator it;
  it=std::lower_bound(x.begin(),x.end(),u);
  int idx=std::max( int(it-x.begin())-1, 0);
  double h=u-x[idx];
  double interpol;
  if(u<x[0]) {
      // extrapolation to the left
      interpol=((b[0])*h + c[0])*h + y[0];
   } else if(u>x[n-1]) {
      // extrapolation to the right
      interpol=((b[n-1])*h + c[n-1])*h + y[n-1];
   } else {
      // interpolation
      interpol=((a[idx]*h + b[idx])*h + c[idx])*h + y[idx];
   }
   return interpol;
}

double Spline::seval_derivative(vector<double> x, vector<double> y, vector<double> a, vector<double> b, vector<double> c,double u)
{
   double h=0.01;
   double fplus,fminus;
   double uplus,uminus;
   double deriv;
   uplus=u+h;
   uminus=u-h;
   fplus=seval(x,y,a,b,c,uplus);
   fminus=seval(x,y,a,b,c,uminus);
   deriv=(fplus-fminus)/(2.0*h);
   return deriv;
}
