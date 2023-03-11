#include <iostream>
#include <vector>
#include <fstream>
#include<chrono>
#include<random>
using namespace std;
double ar[1000][9];


int main()
{
/* Create random engine with the help of seed */
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
    default_random_engine e (seed);

    /* declaring normal distribution object 'distN' and initializing its mean and standard deviation fields. */
  /* Mean and standard deviation are distribution parameters of Normal distribution.
  Here, we have used mean=5, and standard deviation=2. You can take mean and standard deviation as per your choice */
    normal_distribution<double> distN0(0.29,0.000166); // for X
    normal_distribution<double> distN1(7.5e21,8.33e20); // for NA
    normal_distribution<double> distN2(90,3.333333); // for SRV
    normal_distribution<double> distN3(10.0e-6,6.66666667e-7); // for Diffusion
    normal_distribution<double> distN4(10.0e-6,3.33333333e-7); // for thickness
    int n;
    cout<<"\nEnter no. of samples: ";
    cin>>n;
    cout<<"\nNormal distribution for "<<n<<" samples  =>\n";
// VALUE OF n is taken 1000
    for (int i=0; i<n; i++)
    {
      ar[i][0]=distN0(e); // for X
      ar[i][1]=distN1(e);//for NA
      ar[i][2]=distN2(e);//for SRV
      ar[i][3]=distN3(e);//for diffusion
      ar[i][4]=distN4(e);//for thickness
    }
    cout<<"RANDOM VALUES ARE"<<endl;
    for(int i=0;i<n;i++)
    {
        cout<<i<<"\t"<<ar[i][0]<<"\t"<<ar[i][1]<<"\t"<<ar[i][2]<<"\t"<<ar[i][3]<<"\t"<<ar[i][4]<<endl;
    }

    // open a file in write mode.
   ofstream outfile;
   outfile.open("check1.csv");
   outfile<<"X"<<","<<"NA"<<","<<"SRV"<<","<<"Diffusion"<<","<<"Thickness"<<","<<"RA"<<",";
   outfile<<endl;

   cout << "Writing to the file" << endl;
   for(int i=0;i<n;i++)
   {
       for(int j=0;j<5;j++)
       {
           outfile<<ar[i][j]<<",";
       }
       outfile<<endl;
   }


   // write inputted data into the file.
   /*outfile << ar[0]<<",";
   outfile << ar[1] << endl;
   outfile<<ar[2]<<endl; */


   // close the opened file.
   outfile.close();


return 0;
}
