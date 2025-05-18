#include<iostream>
#include<fstream>
#include<random>
#include<ctime>
#include<cstdlib>
#include<armadillo>

using namespace arma; 
using namespace std;

void run_simulation(int seed)
{
  ofstream file;
  //Calling the mersenne twister pseudo random number generator
  mt19937 gen(seed);
  normal_distribution<double> d(0.0, 1.0);


  int N=5000;//Number of IID RV to generate
  vec RV(N);
  vec spacing(N-1);  

  for(int i=0; i<N; i++)
    {
      RV[i]=d(gen);
    }
 
 RV=sort(RV);
 
 file.open("RV_spacing.dat");
 file<<"";
 file.close();
 for(int i=0; i<N-1; i++)
   {
     spacing[i]=RV[i+1]-RV[i];
   }
 file.open("RV_spacing.dat");
 file<<spacing;
 file.close();
 
}

int main()
{
  time_t current_time = time(NULL);
  int seed=(int) current_time;
  
  run_simulation(seed) ;
  return 0;
}
