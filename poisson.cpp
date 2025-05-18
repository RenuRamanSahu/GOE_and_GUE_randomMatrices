#include<iostream>
#include<fstream>
#include<armadillo>



using namespace std;
using namespace arma;


void run()
{
  int N=1000;
  mat D(N, N);


  for(int i=0; i<N; i++)
    {
      D(i, i)=i+1;
      for(int j=i+1; j<N; j++)
	{
          D(i, j)=i+1;
	  D(j, i)=i+1;
	}
    }


  /********Eigen values and Eigen vectors of D**********/

  vec eig_C(N);
  vec eig_D(N);
  mat eig_vec_D(N, N);

  eig_sym(eig_D, eig_vec_D, D);


  for(int i=0; i<N; i++)
    {
      eig_C[i]=1/(eig_D[i]);
    }

  
  //  cout<<"The matrix D is \n"<<D<<"\n";
  // cout<<"The eigen values of D are \n"<<eig_D<<"\n";
  // cout<<"the eigen vectors of D are \n"<<eig_vec_D<<"\n";
  // cout<<"The eigen values of C are \n"<<eig_C<<"\n";
  

  ofstream file;
  file.open("eig_C.dat");
  file<<eig_C;  
  file.close();

  file.open("eig_D.dat");
  file<<eig_D;  
  file.close();

  file.open("eig_vec_D.dat");
  file<<eig_vec_D;  
  file.close();




}


int main()
{
  run();
  return 0;
}
