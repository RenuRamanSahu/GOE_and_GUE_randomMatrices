#include<iostream>
#include<fstream>
#include<ctime>
#include<random>
#include<armadillo>
#include<cstdlib>

using namespace std;
using namespace arma;

void run(int seed)
{
  ofstream file;
  
   //Calling the mersenne twister engine
   mt19937 gen_real(seed), gen_imag(seed+19937);
   double mu=1.0, sigma_sq=1.0;
   normal_distribution<double> distribution(mu, sigma_sq);

   

   
   int N=100;
   int n=N*N;
   cx_mat  A(2*N, 2*N);
   cx_vec eig_A(2*N);

   
   mat M01(2, 2, fill::zeros), M10(2, 2, fill::zeros);
   M01(0, 1)=1;
   M10(1, 0)=1;
   //   cout<<"The matrix M01 is \n"<<M01<<"\n and M10 is \n"<<M10;
   
   cx_mat W(N, N);
   for(int i=0; i<n; i++)
     {
       W[i]={distribution(gen_real), distribution(gen_imag)};
     }
   //  cout<<"the matrix W is \n"<<W<<"\n";


   A=kron(M01, W)+kron(M10, W.t());
   eig_A=eig_gen(A);

   //   cout<<"The matrix A is \n"<<A<<"\n";
   //   cout<<"The eigen values of this martix are \n"<<eig_A;



   vec ch_eig(N);
   vec eig = sort(real(eig_A));



   //eig.dat contains all the eigen values of the matrix A
   file.open("eig.dat", ios::app);
   file<<eig;
   file.close();


   //ch_eig.dat contains only the positive eigen values of the matrix 
   //   cout<< " The eigen values of A are \n " <<eig;
   for(int i=0; i<N; i++)
     {
       ch_eig[i]=eig[N+i];
     }

   file.open("ch_eig.dat", ios::app);
   file<<ch_eig;
   file.close();
   
   /*
   cx_mat matrix=W*W.t();
   cx_vec eig_WW;
   cout<<" The product W*W^{dagger} is \n "<< matrix<<"\n\n";

   eig_WW = eig_gen(matrix);
   cout<<"The eigen values of WW* is \n"<<eig_WW;

   vec eig_W_real = sqrt(real(eig_WW));

      cout<<"The square root of eigen values of W*W are\n "<<eig_W_real;

   */
}




int main()
{

  ofstream file;
  file.open("ch_eig.dat");
  file<<"";
  file.close();



   time_t current_time=time(NULL);
   int seed = (int) current_time;
   cout<<"The seed for the simulation is "<<seed<<"\n";


   int T=1000;

  for(int i=0; i<T; i++)
    {
      run(seed+17*i);
      cout<<i<<"\n";
    }

  return 0;
}
