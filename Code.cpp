#include <iostream>
#include <vector>
#include<chrono>
#include<random>
#include <cmath>
#include <float.h>
#include <fstream>
#include <sstream>
using namespace std;
double Pai	= 3.1415926;
double Q	= 1.6e-19;
double K 	= 1.38e-23;
double MO 	= 9.1E-31;
double H  	= 6.625E-34;
double Hcr 	= 1.0544E-34;
double EO 	= 8.854E-12;
double c	= 3e8;
double MCT_lattice = 6.46e-10;
double ar[1000][9];
//All parameters are in SI units. Except Lamda (Wavelength) which is in microns.
void netd(double F_num, double Aj, double Dstar, double lam1, double lam2, double tint, double Tscene,int p);
double eg(double X, double T)
{
	double eg1 = - 0.302 + 1.93 * X + 0.000535 * T * (1.0 - 2.0 * X);
	double eg2 = - 0.81 * X * X + 0.832 * X * X * X;
	//printf("%f\n",(eg1+eg2));
	return (eg1 + eg2);
}

double ni(double X, double T)
{
	double EG  = eg(X, T);
	double NI1 = 5.585 - 3.82 * X + 0.001753 * T - 0.001364 * X * T;
	double NI2  = 1E+14 * pow(EG,0.75) * pow (T, 1.5);
	double NI3   = exp(-EG / (2.0 * K * T / Q));
	double NII = NI1 * NI2 * NI3;
	//printf("%e\t%e\n", EG, NII);
	return (NI1 * NI2 * NI3 * 1000000);
}

double di(double X)
{
	return 20.5 - 15.5 * X + 5.7 * X * X;
}
double W_dep(double NA, double ND, double X, double T)
{
	double	NI	= ni(X, T),
		DI	= di(X),
		VBI	= K * T / Q * log(NA * ND / (NI * NI));
//		VR	= VBI - V;
	return sqrt((2 * DI * VBI * EO * (NA + ND)) / (Q * NA * ND));
}
double mun(double X, double T)
{
	double MUN1 = pow(0.2 / X, 7.5);
	double MUN2 = 2 * pow(0.2 / X, 0.6);
	double MUN3 = pow(T, MUN2);
	return (9.0E4 * MUN1/MUN3); //MUN = 100000 * .0001;
}
double FLUX_mono (double Lmda, double Tscene)
{
	double C1 = 1.8836e23, C2    = 14388.0;
	//, TH    = 34.0*(PI/180);
	double den   = pow(Lmda, 4.0) * (exp(C2 / (Lmda * Tscene)) - 1);
	return C1 / den;
}
double FLUX_band (double Lmda1, double Lmda2, double TinC)
{
	double T, FL, FU, S=0.0, i, N = 500, Hight;
	T = 273.0 + TinC;
	Hight = (Lmda2 - Lmda1) / N;
	FL = FLUX_mono(Lmda1, T);
	FU = FLUX_mono(Lmda2, T);
	for(i = (Lmda1 + Hight); i<= (Lmda2 - Hight); i+=Hight)  S+= FLUX_mono(i, T);
	double fb = Hight / 2 * (FL + 2 * S + FU) * 1e4;

	return fb;
}
double alpha(double x, double T, double Lamda)
{
	double	EG	= eg(x, T);

	double	SIGMA	= 3.267e4 * (1 + x),
		EGO	= -0.3424 + 1.838 * x + 0.148 * pow(x, 4.0),
		ALFA0	= exp(-18.88 + 53.61 * x),
		W	= 81.9 + T,
		ALFAT	= 100 + 5000 * x;
	double	ENRG	= 1.24 / Lamda,
		ALFA1	= SIGMA * (ENRG - EGO) / W,
		ALFA	= ALFA0 * exp(ALFA1);
	if(ALFA<=ALFAT) return ALFA;
	double	ET	= W * log(ALFAT/ALFA0) / SIGMA + EGO;
	double	BETA	= ALFAT / sqrt(ET - EG);

	return BETA * sqrt(ENRG - EG);

}
double qe1(double Dep_W, double L, double Wavelength, double Thickness, double S, double X, double T) //Main quantum Eff
{
	double	MUN	= mun(X, T);
	double	D	= MUN * K * T / Q;
	double	S1	= (S * L / D),
		ALFA	= alpha(X, T, Wavelength)*100,
		ALFA1	= L * ALFA,
			GP = S1 + ALFA1,
			GN = S1 - ALFA1,
			BP = (1 + S1) * exp( Thickness / L),
			BN = (1 - S1) * exp(-Thickness / L),
				FAC1 = ALFA1 / (ALFA1 * ALFA1 - 1),
				FAC2 = 2 * GP - (BP - BN) * exp(-ALFA * Thickness),
				FAC3 = exp(-ALFA * Thickness) / ALFA1,
					QE1 = FAC1 * (FAC2 / (BP + BN) - FAC3) - exp(-ALFA * Thickness),
					QE2 = exp(-ALFA * Thickness) * (1 - exp(-ALFA * Dep_W)),
						QE = QE1 + QE2;
	return QE;
}


void outputs(double X,double NA,double S,double L,double Th,double RA,int i)
{
    int p=i;
	double 	Lam, Lmst=3.2, Lmend=4.8, Lmstp=0.1;
	double		AREA_junct = (15e-6)*(15e-6);
	double		F	= 3.25;
	//double		RA	= 5000.0;//CHANGE LATER ON
	double		ND	= 1e18*1e6;
	double		EPD = 1e6*1e4;
	double		T	= 80;
	double		Trans_opt= 0.9409;
		double 	Wd = W_dep(NA, ND, X, T);
	double loops = 1+(Lmend-Lmst)/Lmstp;
	double RV=0.0, Resp=0.0, Noise=0.0, Dstar=0.0;
	double Q_e=0.0,Q_eff=0.0;



double 	FLUX=FLUX_band(Lmst, Lmend, 27);
	double Jsh = 1e-12;
    //Lam=4e6;

			for(Lam=Lmst; Lam<=Lmend; Lam+=Lmstp)
			{
				Q_e		= qe1(Wd, L, Lam, Th, S, X, T);
				Q_eff  += Q_e;
				RV     += Q_e*Lam;
				Noise  += sqrt(4*K*T/RA + 2*Q_e*Q*Q*FLUX*1/(4*F*F+1) + 1.25e-23*1e4 + 1e-4*Jsh*Jsh*AREA_junct);
			}
			ar[i][7]=Q_eff;
			Resp = (RV/loops)*Q/(H*c)*1e-6;
			Dstar= Resp/(Noise/loops);
			ar[i][6]=Dstar;
			//variable declaration for netd;
			double lamda1=3,lamda2=5,tint=4e-3,tscene=300.0;
			netd(F,AREA_junct,ar[i][6],lamda1,lamda2,tint,tscene,p);

}
double POWER(double Lam, double Tscene) // Calibrated with Derniak (pg: 285 for 300K)
{
	/*----------------* IMPERICAL RELATION *----------------------*\
	double 	B	= 8.7317e-9,
			X1	= 14388/(LamCut*Tscene);
	return B*pow(Tscene, 4)*(X1*X1*X1 + 3*X1*X1 + 6*X1 + 6)*exp(-X1);
	***************************************************************/
	double	C1 = 37400,
			C2 = 14388,
			var= C2/(Lam*Tscene),
			den= pow(Lam, 5)*(exp(var)-1);
	return	C1/den;
}
double dPdT_mono(double Lam /*microns*/, double Tscene /*Kelvin*/)	//W(cm^-2)(K^-1)
{
	double	C1 = 37400,
			C2 = 14388,
			var= C2/(Lam*Tscene),
			num= C2*exp(var)*POWER(Lam, Tscene),
			den= Lam*Tscene*Tscene*(exp(var)-1);
	return	num/den;
}

double dPdT_band(double Lam1, double Lam2, double Tscene)	//W(cm^-2)(K^-1)
{
	double stp = (Lam2-Lam1)/100;
	double dPdT= 0.0;
	for(double Lam=Lam1; Lam<=Lam2; Lam+=stp)	dPdT += dPdT_mono(Lam, Tscene)*stp;
	return dPdT;
}
void netd(double F_num, double Aj, double Dstar, double lam1, double lam2, double tint, double Tscene,int p)
{
	double dpdT = dPdT_band(lam1, lam2, Tscene)/Pai;	//in Wcm^-2/St
	ar[p][8]=4*F_num*F_num/(Pai*sqrt(Aj)*Dstar*dpdT*sqrt(2*tint));
}
int main()
{
string fname="check1.csv";


vector<vector<string>> content;
vector<string> row;
string line, word;

fstream file (fname, ios::in);
if(file.is_open())
{
while(getline(file, line))
{

row.clear();

stringstream str(line);

while(getline(str, word, ','))
row.push_back(word);
content.push_back(row);
}
}
else
cout<<"Could not open the file\n";

for(int i=0;i<content.size();i++)
{
for(int j=0;j<content[i].size();j++)
{
//cout<<content[i][j]<<" ";
ar[i][j]=strtod(content[i][j].c_str(), nullptr);//string to double
}
cout<<"\n";
}

for(int i=1;i<1002;i++)
{
    outputs(ar[i][0],ar[i][1],ar[i][2],ar[i][3],ar[i][4],ar[i][5],i);
}
//The below is just for printing the array values
for(int i=1;i<1001;i++)
{

    for(int j=0;j<9;j++)
    {

    cout<<ar[i][j]<<",";
    }
    cout<<endl;
}
ofstream outfile;
   outfile.open("output.csv");
   outfile<<"Dstar"<<","<<"Qeff"<<","<<"Netd"<<",";
   outfile<<endl;

   cout << "Writing to the file" << endl;
   for(int i=1;i<1001;i++)
   {
       for(int j=6;j<9;j++)
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


