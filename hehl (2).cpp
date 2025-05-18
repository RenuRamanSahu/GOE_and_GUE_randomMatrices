#include<iostream>
#include<ctime>
#include<random>
#include<cstdlib>
#include<armadillo>
#include<fstream>

using namespace std;
using namespace arma;



void run_simulation(int seed)
{

  double kappa=0.5;
  int seedA=seed; 
  int seedB=seedA*13;
 
  //  cout<<"The seed that will be used is "<<seedA<<"\n";
  
  //Calling the mersenne twister engine
  mt19937 genA(seedA), genB(seedB);
  double mu_A=0.0, mu_B=0.0, sigmasq_A=0.08, sigmasq_B=0.32;
  normal_distribution<double> distributionA(mu_A, sigmasq_A), distributionB(mu_B, sigmasq_B);


  //The following code will generate 8by8 random matrix with complex entries
  cx_mat B(8, 8);
  for(int i=0; i<64; i++)
    {
      B[i]={distributionB(genB), 0*distributionB(genB)};
    }
  B=(B+B.t())*0.5;
  //THE DAGGER OF B WILL BE GIVEN BY B.t()

  cx_mat A_tmp(2, 2);
  for(int i=0; i<4; i++)
    {
      A_tmp[i]={distributionA(genA), 0*distributionA(genA)};
    }
  //the dagger of A_tmp will be given by A_tmp.t()


 
  
  cx_mat I(4, 4);
  I.eye();
  cx_mat A(8,8);
  A=kron(A_tmp, I);//Direct product of A_tmp and I.


  /*for(int i=0; i<64; i++)
    {
      A[i]={distributionB(genB), 0*distributionB(genB)};
      }*/
  //THE DAGGER OF B WILL BE GIVEN BY B.t()


  //symmetrising
  A=(A+A.t())*0.5;


  //Making of the 32 by 32 matrix
  cx_mat M00(4, 4, fill::zeros);
   M00(0, 0)=1.0;
   M00(1,1)=-1.0;
   M00(2,2)=1.0;
   M00(3,3)=-1.0;
  
  

  cx_mat M02(4, 4, fill::zeros);
  M02(0, 2)=1.0;
  M02(1, 3)=1.0;
  
  

  cx_mat M03(4, 4, fill::zeros);
  M03(0, 3)=1.0;
  M03(1, 2)=1.0;
  


  cx_mat Wilson_mat(32,32), Identity(8, 8);
  Identity.eye();

  Wilson_mat=kron(M00,0.5*(1./(2*kappa))*Identity)+kron(M02, (Identity-A)) + kron(M03, B);
  //cout<<Wilson_mat(0, 30)<<"\n";
  

  Wilson_mat=Wilson_mat+Wilson_mat.t();
  //cout<<Wilson_mat(0, 0)<<"\n";

  //cout<<Wilson_mat(0, 30)<<"\n";


  //Following code is for the consistensy check of the  Wilson_matrix

  cx_mat Wm_transpose(32, 32);
  Wm_transpose=Wilson_mat.t();

  int status=0;
  for(int i=0; i<8; i++)
    { 
      for(int j=0; j<8; j++)
	{
          if(Wilson_mat(i, j)!=Wm_transpose(i, j))
	    {cout<<"Something is wrong. Wilson mat is not Hermitian\n";
	      break;
	      status=1;}
	}
      
    }

  /* if(status!=1)
    {
      // cout<<"Hermiticity successfully checked!\n";
      }*/


  //generating the eigenvalues
  cx_vec eig_wilson;
  eig_wilson=eig_gen(Wilson_mat);

  
  
  

 
  ofstream file;
  file.open("wilson_eig.dat", ios::app);
  file<<real(eig_wilson);
  file.close();


  //cout<<"The eigen values are \n";
  //cout<<eig_wilson<<"\n"; 
}


int main()
{
  time_t current_time= time(NULL);
  int seed=(int)current_time;



  ofstream file;
  file.open("wilson_eig.dat");
  file<<"";
  file.close();

  for(int i=0; i<1000; i++)
    {run_simulation(seed+17*i);
      cout<<i<<"\n";}

  return 0;
}
