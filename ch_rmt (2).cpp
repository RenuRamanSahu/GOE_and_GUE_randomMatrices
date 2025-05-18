#include<iostream>
#include<fstream>
#include<ctime>
#include<random>
#include<armadillo>
#include<cstdlib>

using namespace std;
using namespace arma;


double unfold(double x)
{
  int N=100;
   double A, R, term1, term2, ret;
   A=0.165856;
   R=27.6536;
   
   
   if(x<R)
     {
       term1=(A/2)*x*sqrt(R*R-x*x);
       term2=(A/2)*R*R*asin(x/R);
       ret=term1+term2;
       return ret;
     }
   if(x>R)
     {
       ret = N;
     }

}




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



   vec ch_eig(N), unfolded_ch_eig(N);
   vec eig = sort(real(eig_A));

   //ch_eig.dat contains only the positive eigen values of the matrix 
   //   cout<< " The eigen values of A are \n " <<eig;
   for(int i=0; i<N; i++)
     {
       ch_eig[i]=eig[N+i];
       unfolded_ch_eig[i]=unfold(ch_eig[i]);
     }

   file.open("ch_eig.dat", ios::app);
   file<<ch_eig;
   file.close();


   file.open("unfolded_ch_eig.dat", ios::app);
   file<<unfolded_ch_eig;
   file.close();   
   
   /******************
   cx_mat matrix=W*W.t();
   cx_vec eig_WW;
   cout<<" The product W*W^{dagger} is \n "<< matrix<<"\n\n";

   eig_WW = eig_gen(matrix);
   cout<<"The eigen values of WW* is \n"<<eig_WW;

   vec eig_W_real = sqrt(real(eig_WW));

      cout<<"The square root of eigen values of W*W are\n "<<eig_W_real;

   ********************************/
}

void analysis()
{
  ofstream file;
  int T=1000;
  int N=100;
  /***** non-unfolded data analysis *****/
  mat data(T, N);
  mat unfolded_data(T, N);
  mat D(N, N);
  mat u_D(N, N);//u_D stands for unfolded_D




  ifstream infile;
  infile.open("ch_eig.dat");
  for(int i=0; i<T; i++)
    {
      for(int j=0; j<N; j++)
	{
          infile>>data(i, j);
	}
    }
  infile.close();


 
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

  
  


  
  infile.open("unfolded_ch_eig.dat");
  for(int i=0; i<T; i++)
    {
      for(int j=0; j<N; j++)
	{
          infile>>unfolded_data(i, j);
	}
    }
  infile.close();

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

  cout<<"Eigen vectors of unfolded D are computed\n";

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

  ofstream file;
  file.open("ch_eig.dat");
  file<<"";
  file.close();

  file.open("unfolded_ch_eig.dat");
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


  analysis();

  
  return 0;
}
