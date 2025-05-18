#include<iostream>
#include<fstream>
#include<ctime>
#include<random>
#include<armadillo>
#include<complex>


using namespace std;
using namespace arma;

//function for unfolding
/********

This function comes by using the follwoing procedure
1) Generate a large number of instances of N*N GOE (N=100) in this case
2) Diagonalise and store the eigen values
3) Plot the histogram
4) get the envelope of the histogram
5) Fit the function with wigner semi circle law to obtain the density function
6) Multiply the density function by N to obtain R_1(x)
7) Integrate it to get the staircase function
8) Fit with appropriate parameters to obtain A, and R
9) use A, and R in the following function.

*********/
double unfold(double x)
{
  const double pi=3.14159265359;
  double A=0.153886;
  double R=20.5276;
  double term1 = (A/2.0)*x*sqrt(R*R-x*x);
  double term2 = (A/2.0)*R*R*(asin(x/R)+pi/2);
  double ret = term1+term2;
  return ret;
}


void generate_data(double seed, int N, int T)
{
  ofstream file; 
  //calling the pseudo random number generator
  mt19937 gen(seed);
  //calling the gaussian distribution function
  double mu=0.0, sigma_sq=1.0;
  normal_distribution<double> d(mu, sigma_sq);

  
  mat data(T, N);
  mat unfolded_data(T, N);
  mat D(N, N);
  mat u_D(N, N);//u_D stands for unfolded_D

  cx_mat A(N, N);
  cx_vec eig_A(N);
  vec eig_vals;
  
  int n=N*N;

  cout<<"Generating eigen data ... \n";
  for(int j=0; j<T; j++)
    {
      cout<<j<<"\n";
        //Matrix A has been created in which all the elements are from N(0, 1)
        for(int i=0; i<n; i++)
	      {
		A[i]={d(gen), 0.0};
	      }
        //The following operation symmetrises A
	A=(A+A.t())*0.5;

  
	//This matrix has to be diagonalised to get the eigen values
  
	//The command eig_sym() is a command in armadillo which
	//returns the eigen values in ascending order
  
	//	cout<<real(A);

	vec eig_vals;
	eig_A=eig_gen(A);
	eig_vals=sort(real(eig_A));

	for(int i=0; i<N; i++)
	  {
	    data(j, i)=eig_vals[i];
	  }
	
  
    }
  
  cout<<"Eigen data has been generated\n";
  /************************************ We have obtained the data *******************************************/

    
   file.open("data.dat");
  file<<data;
  file.close();
    cout<<"Eigen data has been stored \n";

  
  

  
  /********  This code is for getting the ensemble averages and computing D_ij in non-unfolded scale  *********/

  vec avg_x(N, fill::zeros);
  mat xixj(N, N);//<x_i><x_j>
  mat xx_ij(N, N);//<x_ix_j>

  xixj.zeros();
  xx_ij.zeros();
  
  for(int i=0; i<T; i ++ )
    {
      for(int j=0; j<N; j++)
	{
	  avg_x[j]=avg_x[j]+data(i, j);
	  xx_ij(j, j)= xx_ij(j, j) +  data(i, j)*data(i, j);
          for(int l=j+1; l<N; l++)
	    {
	      xx_ij(j, l)=xx_ij(j, l)+data(i, j)*data(i, l);
	      xx_ij(l, j)=  xx_ij(j, l);
	    }
	}
    }
  xx_ij=xx_ij/T;
//  cout<<"The matrix xx_ij is \n"<<xx_ij<<"\n";

  
  
  avg_x = avg_x/T;
  //cout<<" The average of the eigen values is\n "<<avg_x;

  
  
  for(int i=0; i<N; i++)
    {
      xixj(i, i)  = avg_x[i]*avg_x[i];
      for(int j=i+1; j<N; j++)
	{
	  xixj(i, j) = avg_x[i]*avg_x[j];
	  xixj(j, i) = xixj(i, j);
	}
    }


  //cout<<"The matrix x_ix_j is \n "<< xixj  <<"\n";


  D= xx_ij - xixj;
  cout<<" The matrix D is computed \n";

  
  /*******     The code ends here ... D obtained         **********/
  /**** The following code will compute the eigen values of D * **/

  vec eig_D(N), eig_C(N);
  mat eig_vec_D;
  eig_sym(eig_D, eig_vec_D, D);
  cout<<"Eigen values of D are computed\n";

  cout<<"Eigen vectors of D are computed\n";

  //*********Storing the eigen values of D and eigen vectors of D*********//

  file.open("eig_D.dat");
  file<<eig_D;
  file.close();

  file.open("eig_vec_D.dat");
  file<<eig_vec_D;
  file.close();

  

  
  for(int i=0; i<N; i++)
    {eig_C[i]=1.0/eig_D[i];}

  cout<<"The eigen values of C are computed ";
  
  file.open("eig_C.dat");
  file<<eig_C;
  file.close();

  
  

  /**********   The following code is for unfolding   ***********/

  for(int i=0; i<T; i++)
    {
      for(int j=0; j<N; j++)
	{
	  unfolded_data(i, j)=unfold(data(i, j));          
	}
    }


  cout<<"\n\nThe eigen data has been unfolded  \n ";


  file.open("unfolded_data.dat");
  file<<unfolded_data;
  file.close();

  
  //  cout<<"After unfolding the data is \n"<<unfolded_data;
  /***********      Unfolding completed     **************/



  /**************************************************/
  /**************** Analysis of Unfolded Data  *****************/
  /**************************************************/
  vec u_avg_x(N, fill::zeros);
  mat u_xixj(N, N);//<x_i><x_j>
  mat u_xx_ij(N, N);//<x_ix_j>

  u_xixj.zeros();
  u_xx_ij.zeros();
  
  for(int i=0; i<T; i ++ )
    {
      for(int j=0; j<N; j++)
	{
	  u_avg_x[j]=u_avg_x[j]+unfolded_data(i, j);
	  u_xx_ij(j, j)= u_xx_ij(j, j) +  unfolded_data(i, j)*unfolded_data(i, j);
          for(int l=j+1; l<N; l++)
	    {
	      u_xx_ij(j, l)=u_xx_ij(j, l)+unfolded_data(i, j)*unfolded_data(i, l);
	      u_xx_ij(l, j)=  u_xx_ij(j, l);
	    }
	}
    }
  u_xx_ij=u_xx_ij/T;
  //cout<<"The matrix u_xx_ij is \n"<<u_xx_ij<<"\n";

  
  
  u_avg_x = u_avg_x/T;
  //cout<<" The average of the unfolded eigen values is\n "<<u_avg_x;

  
  
  for(int i=0; i<N; i++)
    {
      u_xixj(i, i)  = u_avg_x[i]*u_avg_x[i];
      for(int j=i+1; j<N; j++)
	{
	  u_xixj(i, j) = u_avg_x[i]*u_avg_x[j];
	  u_xixj(j, i) = u_xixj(i, j);
	}
    }


  //cout<<"The matrix u_x_ix_j is \n "<< u_xixj  <<"\n";


  u_D= u_xx_ij - u_xixj;
  cout<<" The matrix D for the unfolded data is computed \n";

  
  /*******     The code ends here ... u_D obtained         **********/
  /**** The following code will compute the eigen values of D * **/

  vec u_eig_D(N), u_eig_C(N);
  mat u_eig_vec_D;
  eig_sym(u_eig_D, u_eig_vec_D, u_D);
  cout<<"Eigen values of unfolded D are computed\n";

  cout<<"Eigen vectors of unfolded D are \n";

  //*********Storing the eigen values of u_D and eigen vectors of u_D*********//

  file.open("u_eig_D.dat");
  file<<u_eig_D;
  file.close();

  file.open("u_eig_vec_D.dat");
  file<<u_eig_vec_D;
  file.close();

  

  
  for(int i=0; i<N; i++)
    {u_eig_C[i]=1.0/u_eig_D[i];}

  cout<<"The eigen values of unfolded C are computed\n ";
  
  file.open("u_eig_C.dat");
  file<<u_eig_C;
  file.close();


  cout<<"Run Completed!!";

  
}


int main()
{
  //creates an empty file named data.dat
  //This file will store the array of eigen values
  //the rows will be the N eigen values
  //there will be T rows, where T is the number of rows
 
  
  time_t t=time(NULL);
  int seed = (int )t;

 
  
    
  int N=100;//N is the dimension of N*N matrix.
  int T=1000;
  // stores the data in the two files created above row by row;

  generate_data(seed, N, T);


  

  

  
  return 0;
}
