#include<iostream>
#include<fstream>
#include<armadillo>
#include<cstdlib>
#include<ctime>
#include<complex>
#include<random>

using namespace arma;
using namespace std;


int main()
{
  time_t t_now=time(NULL);
  int seed = (int)t_now;
  mt19937 gen(seed);
  normal_distribution<double>d(0.0, 1.0);
  ofstream file;
  int N=1000;
  int n=N*N;

  cx_mat A(N, N);
  cx_vec eig_val;

  for(int i=0; i<n; i++)
    {
      A[i]={d(gen), 0};
    }

  //code for symmetrising
  A=0.5*(A+A.t());
  eig_val=eig_gen(A);

  
  file.open("eig.dat");
  file<<real(eig_val);
  file.close();

  
  return 0;
} 
