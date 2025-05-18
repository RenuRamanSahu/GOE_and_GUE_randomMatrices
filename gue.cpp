#include<iostream>
#include<cstdlib>
#include<random>
#include<random>
#include<ctime>
#include<armadillo>
#include<complex>
#include<fstream>


using namespace arma;
using namespace std;

void run_simulation()
{
  cout<<"Hi\n";
  //generating seed for the PRNG(Pseudo random number generator)
  time_t current_time=time(NULL);
  int seed=(int)current_time;//The current time is the seed 
  cout<<"The time since epoch is "<<seed<<"\n";

  //The Mersenne twister Pseudo random generator is called
  mt19937 gen_real(seed), gen_imag(seed*2);
  double mu=0, sigma_sq=1.0;
  normal_distribution<double>distribution(mu, sigma_sq);

  int N;//Dimension of the matrix
  int T;//Number of matrices
  int m;
  ofstream ofile;

  
  N=200;
  m=N*N;
  T=500;
  

  
  cx_mat A(N, N);
  cx_vec eig_val;


  ofile.open("eig_vals_gue.dat");
  ofile<<"";
  ofile.close();

  ofile.open("eig_vals_gue.dat", ios::app);
  for(int i=0; i<T; i++)
    {
      for(int j=0; j<m; j++)
	{
          A[j]={distribution(gen_real), distribution(gen_imag)};
	}
      //      cout<<"Getting a matrix \n"<<A<<"\n";
      //cout<<"After hermitising..."<<i<<"\n";
      A=(A+A.t())/2.0;//A is now hermitian
      //      cout<<A<<endl;

      eig_val=eig_gen(A);
      ofile<<real(eig_val);

      //      cout<<"The eigen values of the hermitised matrix is\n";
      cout<<"Diagonalised Matrix ............"<<i<<"\n";
      
      
      
    }
  ofile.close();

  

}

int main()
{
  run_simulation();
  return 0;
}
