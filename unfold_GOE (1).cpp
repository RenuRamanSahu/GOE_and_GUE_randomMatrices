#include<iostream>
#include<cstdlib>
#include<ctime>
#include<random>
#include<fstream>
#include<armadillo>
#include<complex>

using namespace arma;
using namespace std;

int seed()
{
  int s;
  time_t current_time=time(NULL);

  s=(int)current_time;
  //  cout<<"The seed for this simulation is "<< s<<"\n";
  
  return s;
}


void run(int s)//s is the seed
{
  
  mt19937 gen(s);
  normal_distribution<double>d(0.0, 1.0);
  int N=100;
  int n=N*N;
  ofstream file;
  cx_mat A(N, N);
  cx_vec eig_val;

  
  
  
  for(int i=0; i<n; i++)
    {
      A[i]={d(gen), d(gen)};
    }

  A=0.5*(A+A.t());

  eig_val = eig_gen(A);
  file.open("eig_GOE.dat", ios::app);
  file<<real(eig_val);
  file.close();
  

}


int main()
{
  int T=1000;
  ofstream file;
  file.open("eig_GOE.dat");
  file<<"";
  file.close();

  int s= seed();
  for(int i=0; i<T; i++)
      {
        run(s+i);
        cout<<i<<"\nh";
      }

  return 0;
}
