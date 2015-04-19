// NGDE- Aerosol Dynamics solver.
// Authors: Anand Prakash, Ameya P. Bapat and Michael R. Zachariah, Mar 27, 2003.

// This C code is used to solve for nucleation, coagulation and surface growth problems. This example problem solves for characteristics of an Aluminum aerosol. However, it can be used for any other material whose properties are known. There are four different sections in this code: 
//1) Coagulation, 
//2) Nucleation (classical) + coagulation, 
//3) Surface growth 
//4) Unified GDE with all the three phenomena combined. 

// The above phenomena constitute four different sections of the code. Each of the four sections are independent of each other. They have especially been written in such manner so that, one can identify the contribution of each phenomena to the GDE. Also, if one requires to use nucleation, coagulation or surface growth alone, it can be easily done.   

// The solution algorithm involves a new approach to solve the Aerosol GDE, where the particle volume space is divided into nodes of zero width. Using nodes allows us to cover the entire size range of aerosol (1 nm to 10 micrometer) using just 40 nodes.A more detailed description of the theory and algorithm can be found in the reference:
//A. Prakash, A. P. Bapat and M.R. Zachariah,  Aerosol Science and Technology

//The main theoretical constraints of this implementation are:
//1. Use of the Self Consistent Classical nucleation Model
//2. Free molecule collision kernel
//3. Free molecule surface growth

// The main variables in the code are v[MAX] and N [MAX]. These arrays are respectively the size (volume), and the number concentration of particles corresponding to that size. v1 is the monomer volume, which we assume, is equal to the volume of a molecule. We have tried to comment the code sufficiently to enable the user to make changes as needed. The main variables used are listed in the table below.

//**************************Note on compiling and running the code**************************
//To compile the code, use command "gcc gde.c -lm". To run the code the user should edit the input file as required and then the command "a.out > outputfile" should be used to direct the output to a file.
//******************************************************************************************

// *************************************NOTE ON UNITS***************************************
//                            SI Units are used throught the code.
//******************************************************************************************
/*                                   ------------                                            */
/* 				     NOMENCLATURE                                            */
/* 				     ------------                                            */

/* __________________________________________________________________________________________*/
/* Variables            	            Meaning (Units)                                  */
/* __________________________________________________________________________________________*/
/* v[i]         	Volume of particles at ith node (m3/m3 of aerosol)                   */
/* N[i]         	Number of particles at ith node (#/m3 of aerosol)                    */
/* X[i][j][k]   	Splitting operator for coagulation (dimensionless)                   */
/* K[i][j]      	Collision frequency function(m3/s)                                   */
/* N1s[I]       	Saturation monomer concentration over particle of size I(#/m3)       */
/* sum1,sum2    	Addition and subtraction terms in Nk due to collision(#/(m3s))       */
/* Vav          	Average volume of particles (m3/m3 of aerosol)                       */
/* Vtotal       	Total volume of particles and monomer (m3/m3 of aerosol)             */
/* T                    Temperature	(K)                                                  */
/* coolrate     	Cooling rate	(K/s)                                                */
/* t            	Time	(s)                                                          */
/* step	                Timestep for integration  (s)                                        */
/* Ninf	                Total number of particles as got by solving GDE	(#/m3 of aerosol)    */
/* Ninf_eq       	Total number of particles by eq 7.77, Friedlander (#/m3 of aerosol)  */
/* fi           	Total volume of particles (m3/m3 of aerosol)                         */
/* Vtot	                Total volume of particles (excluding monomers)	(m3/m3 of aerosol)   */
/* addterm,subterm	Addition and subtraction terms in Nk due to surface growth(#/(m3s))  */
/* Ntot         	Total number of particles (#/m3 of aerosol)                          */
/* dpav         	Mean Diameter	(m)                                                  */
/* v1           	Monomer volume	(m3/m3 of aerosol)                                   */
/* sigma        	Surface Tension	(N/m)                                                */
/* Ps           	Saturation Pressure (Pa)                                             */
/* ns           	Monomer concentration at saturation(#/m3 of aerosol)                 */
/* kstar        	Critical node size (dimensionless)                                   */
/* dpstar       	Critical particle size (m)                                           */
/* vstar        	Volume of particles in the critical bin	(m3/m3 of aerosol)           */
/* s1           	Surface area of monomer	(m2)                                         */
/* S            	Saturation ratio (dimensionless)                                     */
/* m1           	Monomer mass(Kg)                                                     */
/* zeta[i]      	Splitting operator for nucleation at the ith node (dimensionless)    */
/* t_SPD                Time requred to reach a self-preserving distribution (s)             */

//*****************************EXAMPLE PROBLEM*********************************
// Aluminum vapor at a high temperature (1773 K) is being cooled at a rate ~1000 K/s. As the system cools down, supersaturation of vapor causes nucleation of particles, which then coagulate to form larger particles. Also, due to high concentration of monomer units (molecules), supersaturation causes surface growth of existing particles. The monomers either evaporate from the surface of existing particles, or condense on them.
//******************************************************************************

// Authors: Anand Prakash, Ameya P. Bapat and Michael R. Zachariah, Nov 27, 2002.
// Beginning of the code
#include<stdio.h>
#include<math.h>

#define	pi	3.1415926	//pi
#define kb  	1.38e-23	//Boltzmann constant
#define R       8.314           //Universal gas constant Joule/(K.mol)
#define Na      6.0225e23       //Avogadro's Number
//#define MAX1    45
#define MAX1    200
#define CONSR  8.31447 
#define Q 1.6021892E-19 //one electron charge
#define	STEPT 0.0001 //delta x to calculate temperature gradient
#define ALFA 1.142
#define	BETA 0.558
#define	GAMA 0.999
#define	FI 0.491
#define	Cm 1.14
#define	Cs 1.17
#define	Ct 2.18
#define M 2
#define MAX2 180 //the maximum array size for this program 180

double v[MAX1 + 1], VV[MAX1], SA[MAX1], Sa[MAX1], N[MAX1], X[MAX1][MAX1][MAX1],
		K[MAX1][MAX1], N1s[MAX1], SURF_MONO[MAX1], sum1, sum2, sum22, sum3,
		sum4, temp, temp1, temp2, Vav, Vtotal, fs_tio2A[MAX1];
double Y[3], DD[3];
double EF, ROU_P, DIA_P, EPS;
//NSTEP is the number of grids in X(i) position-----for 20,30torr case
//  int NSTEP=148;//for 40torr case
int NSTEP = 171;

//hong------------------functions----------------------------------------------
//       one dimension explode fuction
double enlgr(double x[], double y[], int n, double rrt)

//       x---x value being given
//       y---y value being given
//       n---total mumber of the known function
//       t---x value of the unknown point 
//       z---y value of the unknown point
{
	int i, j, k, m;
	double s, z;

	z = 0.0;

	if (n <= 0)
		return (z);
	if (n == 1) {
		z = y[1];
		return (z);
	}
	if (n == 2) {
		z = (y[1] * (rrt - x[2]) - y[2] * (rrt - x[1])) / (x[1] - x[2]);
		return (z);
	}

	i = 1;

	loop: if (x[i] < rrt) {
		i++;
		if (i <= n)
			goto loop;
	}
	k = i - 4;
	if (k < 1)
		k = 1;
	m = i + 3;
	if (m > n)
		m = n;
	for (i = k; i <= m; i++) {
		s = 1.0;
		for (j = k; j <= m; j++) {
			if (j != i)
				s = s * (rrt - x[j]) / (x[i] - x[j]);
		}
		z = z + s * y[i];
	}
	return (z);
}
//----------------------------------------------------------------------------------------------------------------------

//Testing GIT commits

void F(double z_P[], double z_DE[], double z_V[], double z_T[], double z_DV[],
		double z_TC[], double z_WM[], int nn) {
	double LAMDA, factK, Cdrag, MIU, ROU_F, VELM, POSI_P, TEMP_P, VELO_P,
			POSI_PT, TEMP_PT, G_TEMP, CONDF, CONDP, WMO;
	double forceSD, forceES, forceTP, C1, factor1, factorK, factor2;

//-----MIU, VELM, CONDF,CONDP ARE UNKNOWN YET------------------------
//	MIU=2E-5
// 	CONDF=0.06
//-----thermal conductivity of TiO2
	CONDP = 3.36;
//-----thermal conductivity of Al2O3
//      CONDP=6.28;

	C1 = pi * pow(DIA_P, 3.0) * (ROU_P) / 6.0;

//-----%OBTAIN CORRECTION FACTOR C FOR DRAG FORCE%---------------
	POSI_P = Y[2];
	MIU = enlgr(z_P, z_DV, NSTEP + 1, POSI_P);
	CONDF = enlgr(z_P, z_TC, NSTEP + 1, POSI_P);
	WMO = enlgr(z_P, z_WM, NSTEP + 1, POSI_P);
	VELO_P = enlgr(z_P, z_V, NSTEP + 1, POSI_P);
	ROU_F = enlgr(z_P, z_DE, NSTEP + 1, POSI_P);
	TEMP_P = enlgr(z_P, z_T, NSTEP + 1, POSI_P);

	VELM = sqrt(8.0 * CONSR * TEMP_P / pi / (WMO * 0.001));

	LAMDA = MIU / (FI * ROU_F * VELM);
//    printf("\n%e",LAMDA);
	factK = 2.0 * LAMDA / DIA_P;

	factor1 = factK;
//	write(*,*) VELM, MIU, LAMDA, factK
	Cdrag = 1.0 + factK * (ALFA + BETA * exp(-GAMA / factK));

	factor2 = Cdrag;
	forceSD = -3.0 * MIU * DIA_P * (Y[1] - VELO_P) / Cdrag;

//-----%OBTAIN THERMOPHEROTIC FORCE%-----------------------------

	POSI_PT = POSI_P + STEPT;
	TEMP_PT = enlgr(z_P, z_T, NSTEP + 1, POSI_PT);
	G_TEMP = (TEMP_PT - TEMP_P) / STEPT;
	forceTP = -6.0 * pi * MIU * MIU / ROU_F * DIA_P * Cs
			* (CONDF / CONDP + Ct * factK) * G_TEMP / TEMP_P
			/ (1.0 + 3.0 * Cm * factK)
			/ (1.0 + 2.0 * CONDF / CONDP + 2.0 * Ct * factK);

//-----% OBTAIN ELECTROSTATIC FORCE%-----------------------------
	forceES = -Q * EF;

	DD[1] = (forceSD + forceTP + forceES) / C1;

//	WRITE(*,*) forceTP, forceES, forceSD
	DD[2] = Y[1];

}

//*********************************************************************************************
void RKT3(double T, double H, double y_P[], double y_DE[], double y_V[],
		double y_T[], double y_DV[], double y_TC[], double y_WM[], int mm) {
	double A[4 + 1], B[M + 1], C[M + 1], G[M + 1], P, X, HH, DT, TT, QQ;
	int N, i, j, k;

	HH = H;
	N = 1;
	P = 1.0 + EPS;
	X = T;
	for (i = 1; i <= M; i++)
		C[i] = Y[i];
	while (P >= EPS) {
		A[1] = HH / 2.0;
		A[2] = A[1];
		A[3] = HH;
		A[4] = HH;
		for (i = 1; i <= M; i++) {
			G[i] = Y[i];
			Y[i] = C[i];
		}

		DT = H / N;
		T = X;

		for (j = 1; j <= N; j++) {
			F(y_P, y_DE, y_V, y_T, y_DV, y_TC, y_WM, NSTEP + 1);
//	  printf("\n%e\t%e\t%e\t%e",Y[1],Y[2],DD[1],DD[2]);
			for (i = 1; i <= M; i++)
				B[i] = Y[i];
			for (k = 1; k <= 3; k++) {
				for (i = 1; i <= M; i++) {
					Y[i] = Y[i] + A[k] * DD[i];
					B[i] = B[i] + A[k + 1] * DD[i] / 3.0;
				}

				TT = T + A[k];
				F(y_P, y_DE, y_V, y_T, y_DV, y_TC, y_WM, NSTEP + 1);
			}
			for (i = 1; i <= M; i++)
				Y[i] = B[i] + HH * DD[i] / 6.0;
			T = T + DT;
		}
		P = 0.0;
		for (i = 1; i <= M; i++) {
			QQ = fabs(Y[i] - G[i]);
			if (QQ > P)
				P = QQ;
		}

		HH = HH / 2.0;
		N = N + N;
	}
	T = X;
}

//----------------main program--------------------------------------------------------------------------------------  
int main() {
	// Declaration of variables
	int i, j, jj, k, choice, counter, counter2, opt, printopt, flag, MAX,
			beta_option, FLAG_RK;
	int nn1, nn0, nn_count, mstep;
//  double v[MAX1+1],VV[MAX1],SA[MAX1],Sa[MAX1],N[MAX1],X[MAX1][MAX1][MAX1],K[MAX1][MAX1],N1s[MAX1],SURF_MONO[MAX1],sum1,sum2,sum22,sum3,sum4,temp,temp1,temp2,Vav,Vtotal;
	double t, step, temp3, Ninf, n, Vtot, Stot, addterm, subterm, Ntot, dpav,
			dpap, v1, fi, t_SPD, step_coag, q, P;
	double sigma, Ps, ns, theta, a, b, Jk, kstar, dpstar, vstar, s1, S, m1, c,
			zeta[MAX1], dp[MAX1], t1, t2, t3, t4, Kmin;
	double A, B, C, D, E, MW, rho, Ninf_eq, T, coolrate, d;
	double lambda, kn1, kn2, mu, D1, D2, m[MAX1], c1, c2, g1, g2, l1, l2, A_mu,
			B_mu;
	FILE *fptr, *fptr1, *fptr2, *fcoag, *fcoagnu, *fconusurf, *fsurf, *fdata,
			*fparticle, *fextra1, *fextra2;
//hong-----------------------------------------------------------------------
	//double SN1,SN2,SN3,SN4,SN5,SN6,SN7,SN8,SN9,SN10,SN11;
	double SN1, SN2, SN3, SN4, SN5, SN6, SN7, SN8, SN9, SN10;
	double POSIF[MAX2], DENSITYF[MAX2], VELOF[MAX2], TEMPF[MAX2], DV[MAX2],
			TC[MAX2], WM[MAX2];
	double tRKT, Nmax, Vmax, t0, SIGMAgeo, Vgeo;
	double POSI_PM, MIU_M, CONDF_M, WMO_M, VELO_PM, ROU_FM, TEMP_PM, VELM_M;
	double T_FLAME, V_FLAME, RHO_FLAME, RK, RKg, RKs, te_mono, te_mono1,
			time_mono, A_vk, CNT_MONO1, CNT_MONO2, YY_TEMP, CNT_MONO0, CNT_SUM;
	double fs_tio2, ppd, temNp, coeff, coeff1, coalstep, Df;
//  double fs_tio2A[MAX1];

	Df = 3.0;
	if ((fptr = fopen("main.inp.txt", "r+")) == NULL) //Opening the input file that contains property data for the aerosol material.
	{
		printf("Warning: Could not open file");
		return (-1);
	}
	//Scanning the property data from the input file "main.inp"
	fscanf(fptr,
			"Number of nodes: %d Starting Temperature(Kelvin): %lf Cooling Rate(Kelvin/s): %lf Pressure(Atm): %lf Molecular Weight(Kg/mol): %lf Density(Kg/m3): %lf Surface Tension in the form A-BT(dyne/cm): A=%lf B=%lf Saturation Vapor Pressure of material in form Ps=exp(C-D/T)P: C=%lf D=%lf Enter your choice - 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:%d Choice of collision frequency function:1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime:%d Sutherland's constants for carrier gas in the form A*T^1.5/(B+T): A=%lf B=%lf",
			&MAX, &T, &coolrate, &P, &MW, &rho, &A, &B, &C, &D, &choice,
			&beta_option, &A_mu, &B_mu);
	fclose(fptr);  //Closing the main property data file.
//hong---read the database from spin output

	fdata = fopen("DATATEST2.dat", "r");
	fparticle = fopen("PARTICL.dat", "w");
	fextra1 = fopen("primary.dat", "w");
	fprintf(fparticle,
			"\nPOSITION(m)\t\tDENSITY(Kg/m3)\t\tVELOCITY(m/s)\t\tTEMP(K)\t\tMWEIGHT\t\tD_VISCOSITY(Kg/m/S)\t\tTHERMO_CONDUC(W/m/K)");

	for (j = 1; j <= NSTEP; j++) {
		fscanf(fdata, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",
				&POSIF[j], &DENSITYF[j], &VELOF[j], &TEMPF[j], &SN1, &SN2, &SN3,
				&SN4, &SN5, &SN6, &SN7, &SN8, &SN9, &SN10, &DV[j], &TC[j]);

		
		//fscanf(fdata,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",&POSIF[j],&DENSITYF[j],&VELOF[j],&TEMPF[j],&SN1,&SN2,&SN3,&SN4,&SN5,&SN6,&SN7,&SN8,&SN9,&DV[j],&TC[j],&SN10,&SN11);
		POSIF[j] = POSIF[j] / 100.0;
		VELOF[j] = VELOF[j] / 100.0;
		DENSITYF[j] = DENSITYF[j] * 1000.0;
		WM[j] = SN1 * 2.0 + SN2 * 32.0 + SN3 * 18.0 + SN4 * 1.0 + SN5 * 16.0
				+ SN6 * 17.0 + SN7 * 33.0 + SN8 * 34.0 + SN9 * 28.0
				+ SN10 * 2.0;
		DV[j] = 0.1 * DV[j];
		TC[j] = 1E-5 * TC[j];
		fprintf(fparticle, "\n%e\t%e\t%e\t%e\t%e\t%e\t%e", POSIF[j],
				DENSITYF[j], VELOF[j], TEMPF[j], WM[j], DV[j], TC[j]);
	}
	fclose(fdata);
//end reading------------------------------------------------------------------------------     
	q = pow(10.0, (8.0 / (MAX - 2))); //Calculation of geometric spacing factor which depends on the number of nodes.
	v1 = MW / (rho * Na); // Volume of a monomer unit (a molecule); Na is the Avogadro's Number.
//hong---------initialization--------------------------------------------------------------
	EF = 0.0 / 0.04;  //electric field for 4cm gap
//      EF=0.0/0.02; 
	ROU_P = 3.84e3; //density of anatase TiO2
//	ROU_P=3.7E3; //density of gamma AL2O3
	tRKT = 0.0;
	Y[1] = -2.5;  //initial velocity
//	Y[1]=-1.85;
	Y[2] = 3.850000e-02;  //initial  position
	CNT_MONO0 = 9.0e21;
//	CNT_MONO0=1.1927e2;
	CNT_MONO1 = CNT_MONO0;
	CNT_SUM = 0;
	YY_TEMP = Y[2];
	nn0 = 0;
	nn1 = 0;
	nn_count = 0;
	t0 = 0.0;
	FLAG_RK = 0;
//	H=0.000001;
	EPS = 1.0e-06;  //------tolerence--
	Nmax = 0.0;
	Vmax = 0.0;

	//**********************************************
	//           PURE COAGULATION
	//**********************************************
	// This part of the code solves for coagulation problems. The code verifies that we obtain a self preserving distribution. It verifies the problem stated in Friedlander, 2000 - page 218, line 4 - and shows that SPD is obtained. And note that in this case node 1 is not the monomer. Its particle of size 1 nm. We do not have to worry about monomers, since we are concerned with coagulation only.
	if (choice == 1) {
		fcoag = fopen("coag.dat", "w");
		if ((fptr1 = fopen("coag.inp", "r+")) == NULL) //opening the input file for coagulation probelm.
		{
			printf("Warning: Could not open file");
			return (-1);
		}
		for (i = 1; i <= MAX - 1; i++) {
			N[i] = 0.0;
			SA[i] = 0.0;
			VV[i] = 0.0;
		}

		//Scanning the initial monomer size and number concentration of the monodisperse distribution for coagulation.
		fscanf(fptr1,
				"\n Initial diameter of monodisperse particles(nm): %lf  Initial number concentration(#/m3): %lf Temperature(Kelvin): %lf",
				&d, &N[1], &T);
		fprintf(fparticle,
				"\nt(s)\t\tt0(s)\t\tY[1](m/s)\t\tY[2](m)\t\tT(K)\t\tDIA_P(m)\t\tstep(s)\t\tSIGMAgeo");
		v[1] = 1e-27 * pi * pow(d, 3) / 6; // Setting diameter of particle that the user enters as the smallest node( coagulation would never lead to decrease in particle size).
		dp[1] = pow((6 * v[1] / pi), 0.3333333);
		m[1] = v[1] * rho;
		Sa[1] = 1e-18 * pi * d * d;
		SA[1] = N[1] * Sa[1];
		VV[1] = N[1] * v[1];
		for (i = 2; i <= MAX; i++) // Calculating volumes of larger nodes using a geometric spacing factor of q.
				{
			v[i] = v[1] * pow(q, (i - 1));
			dp[i] = pow((6 * v[i] / pi), 0.3333333);
			m[i] = v[i] * rho;
			Sa[i] = Sa[1] * pow(q, (i - 1));
		}
		//Calculation of size splitting operators.
		for (k = 1; k <= MAX - 1; k++) {
			for (i = 1; i <= MAX - 1; i++) {
				for (j = 1; j <= MAX - 1; j++) {
					//Conditions in parentheses check if the combined volume of colliding particles is between k and k+1.
					if (v[k] <= (v[i] + v[j]) && (v[i] + v[j]) < v[k + 1]) {
						X[i][j][k] = (v[k + 1] - v[i] - v[j])
								/ (v[k + 1] - v[k]);
					} else {
						if (v[k - 1] <= (v[i] + v[j]) && (v[i] + v[j]) < v[k]) {
							X[i][j][k] = (v[i] + v[j] - v[k - 1])
									/ (v[k] - v[k - 1]);
						} else
							X[i][j][k] = 0;
					}
				}
			}
		}
		t = 0.0; // Initializing time.
		T = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
		step = 1e-12; //Initializing timestep for integration.
		t_SPD = 5
				/ (pow(3.0 / (4.0 * pi), 0.1666666)
						* pow((6.0 * kb * T / rho), 0.5) * pow(v[1], 0.1666666)
						* N[1]); //Estimation of time required to reach SPD.
//      printf("\n\nEstimated time required to reach SPD=%es\n",t_SPD);
		fprintf(fcoag, "\n\nEstimated time required to reach SPD=%es\n", t_SPD);
		Ninf_eq = 1e24; // A temporary variable that calculates N_total according to Friedlander, 2000, equation 7.77.
//      printf("\nt\t\tNinf\t\tN[1]\t\tN[5]\t\tN[10]\t\tN[15]");
		fprintf(fcoag,
				"\nt\t\tNinf\t\tN[1]\t\tN[5]\t\tN[10]\t\tN[15]\tN[20]\t\tN[25]\t\tN[30]\t\tN[35]\t\tN[40]");
//      while(t<t_SPD)//Running the coagulation mode of the code until SPD is reached. The dimensionless number distribution remains same after SPD has reached.
//      while(t<1e-4)
		while (Y[2] >= 0.00005) {
			Nmax = 0.0;
			Vmax = 0.0;
			//Calculation of collision frequency function ( same as in " nucleation + coagulation" section of the code).
			Kmin = 1e-9; //setting an arbitrary value for Kmin to start with.
			T_FLAME = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
			V_FLAME = enlgr(POSIF, VELOF, NSTEP + 1, Y[2]);
			if (beta_option == 1) {
				for (i = 1; i <= MAX - 1; i++) {
					for (j = 1; j <= MAX - 1; j++) {

						temp1 = 1 / v[i] + 1 / v[j];
						temp2 = pow(v[i], 0.333333) + pow(v[j], 0.333333);
						K[i][j] = pow(3.0 / (4.0 * pi), 0.1666666)
								* pow((6.0 * kb * T / rho), 0.5)
								* pow(temp1, 0.5) * pow(temp2, 2.0);
						if (K[i][j] < Kmin)
							Kmin = K[i][j]; //Calculating the smallest collision frequency function to decide the characteristic coagulation time.

					}
				}
			}
			if (beta_option == 2) {
				mu = A_mu * pow(T, 1.5) / (B_mu + T);
				lambda = (mu / (P * 101325.0))
						* sqrt(pi * R * T / (2.0 * 0.04));

				for (i = 1; i <= MAX - 1; i++) {
					for (j = 1; j <= MAX - 1; j++) {

						kn1 = (2.0 * lambda) / dp[i];
						kn2 = (2.0 * lambda) / dp[j];
						D1 = (kb * T) / (3.0 * pi * mu * dp[i])
								* ((5.0 + 4.0 * kn1 + 6.0 * kn1 * kn1
										+ 18.0 * kn1 * kn1 * kn1)
										/ (5.0 - kn1 + (8.0 + pi) * kn1 * kn1));
						D2 = (kb * T) / (3.0 * pi * mu * dp[j])
								* ((5.0 + 4.0 * kn2 + 6.0 * kn2 * kn2
										+ 18.0 * kn2 * kn2 * kn2)
										/ (5.0 - kn2 + (8.0 + pi) * kn2 * kn2));
						c1 = sqrt((8.0 * kb * T) / (pi * m[i]));
						c2 = sqrt((8.0 * kb * T) / (pi * m[j]));
						l1 = (8.0 * D1) / (pi * c1);
						l2 = (8.0 * D2) / (pi * c2);
						g1 = (pow((dp[i] + l1), 3)
								- pow((dp[i] * dp[i] + l1 * l1), 1.5))
								/ (3.0 * dp[i] * l1) - dp[i];
						g2 = (pow((dp[j] + l2), 3)
								- pow((dp[j] * dp[j] + l2 * l2), 1.5))
								/ (3.0 * dp[j] * l2) - dp[j];
						K[i][j] = 2.0 * pi * (D1 + D2) * (dp[i] + dp[j])
								/ ((dp[i] + dp[j])
										/ (dp[i] + dp[j]
												+ 2.0 * sqrt(g1 * g1 + g2 * g2))
										+ (8.0 * (D1 + D2))
												/ (sqrt(c1 * c1 + c2 * c2)
														* (dp[i] + dp[j])));
						//  printf("kn %e\n",kn2);
						if (K[i][j] < Kmin)
							Kmin = K[i][j]; //Calculating the smallest collision frequency function to decide the characteristic coagulation time.

					}
				}
			}

			//Calculating the gain and loss terms due to coagulation.
			nn_count = nn_count + 1;

			for (k = 1; k <= MAX - 1; k++) {
				sum1 = 0; //Addition term when i and j collide to form a k sized particle.
				sum2 = 0; //Subtraction term when k collides with any other particle.
				sum3 = 0; //addition term to surface area.
				sum22 = 0; //subtraction term to surface area.
				sum4 = 0; //subtraction term to surface area due to coalescence.
				for (i = 1; i <= MAX - 1; i++) {
					sum2 = sum2 + K[k][i] * N[i];
					sum22 = sum22 + K[k][i] * N[i];
					for (j = 1; j <= k; j++) {
						sum1 = sum1 + X[i][j][k] * K[i][j] * N[i] * N[j];
						sum3 = sum3
								+ X[i][j][k] * K[i][j] * N[i] * N[j] * Sa[j]
										* pow(q, (k - j)) * (Sa[i] + Sa[j])
										/ Sa[k]
										/ (pow(q, (i - k)) + pow(q, (j - k)));
					}
				}

				printf("\n\n%i\t%e\t%e\t%e", k, Sa[k] * N[k], SA[k],
						6.0 * VV[k] / SA[k]);
				ppd = 6.0 * v[k] / Sa[k];
				fs_tio2 = 7.4375e16 * pow(ppd, 4) * T_FLAME
						* exp(31030.0 / T_FLAME);
				sum4 = 1.0
						* (SA[k]
								- N[k] * pi
										* pow((6.0 * v[k]) / pi, 0.666666667))
						/ fs_tio2;
				if (k == 1)
					sum4 = 0;
				if (sum4 < 0)
					sum4 = 0;

				SA[k] = SA[k] + step * (0.5 * sum3 - Sa[k] * N[k] * sum22);
				SA[k] = SA[k] - step * sum4;
				VV[k] = VV[k] + step * (0.5 * sum1 - N[k] * sum2) * v[k];

				N[k] = N[k] + step * (0.5 * sum1 - N[k] * sum2);
				if (k > 1) {
					if (SA[k] < N[k] * pi * dp[k] * dp[k]) {
						SA[k] = N[k] * pi * dp[k] * dp[k];
					}
					if (SA[k] > N[k] * Sa[1] * pow(q, (k - 1))) {
						SA[k] = N[k] * Sa[1] * pow(q, (k - 1));
					}
					if (SA[k] < N[k] * Sa[k - 1] * pow(q, 0.6666666667)) {
						SA[k] = N[k] * Sa[k - 1] * pow(q, 0.6666666667);
					}
					if (SA[k] > N[k] * Sa[k - 1] * q) {
						SA[k] = N[k] * Sa[k - 1] * q;
					}
				}
				printf("\n%i\t%e\t%e\t%e\t%e", k, ppd, SA[k], N[k], Y[2]);

			}

			for (i = 2; i <= MAX - 1; i++) {
				if (N[i] != 0)
					Sa[i] = SA[i] / N[i];
			}
			if (nn_count == 1000) {
				for (i = 1; i <= MAX - 1; i++)
					fprintf(fcoag, "\n%i\t%e\t%e\t%e\t%e\t%e", i,
							6.0 * v[i] / Sa[i], Sa[i], Sa[i + 1] / Sa[i], N[i],
							Y[2]);
				nn_count = 0;
			}
//hong---addition of monomers from precursor pyrolyze

			if (CNT_SUM <= 0.99) {
				RK = 3.96e5 * exp(-8479.7 / T_FLAME) * 10.0;
				time_mono = -(YY_TEMP - Y[2]) / V_FLAME;
				N[1] = N[1] + RK * CNT_MONO1 * time_mono;
				SA[1] = SA[1] + Sa[1] * RK * CNT_MONO1 * time_mono;
				VV[1] = VV[1] + v[1] * RK * CNT_MONO1 * time_mono;
				CNT_MONO2 = CNT_MONO1 * (1.0 - RK * time_mono);

				CNT_SUM = CNT_SUM + RK * CNT_MONO1 * time_mono / CNT_MONO0;
				CNT_MONO1 = CNT_MONO2;
				YY_TEMP = Y[2];
//	  fprintf(fcoag,"\n%e\t%e\t%e\t%e",t,CNT_SUM,CNT_MONO1,Y[2]);

			}
			Ninf = 0; // Total number of particles, which is zero initially.
			for (i = 1; i <= MAX - 1; i++)
				Ninf = Ninf + N[i]; //N_infinity as we get by solving the coagulation part of GDE.
			step = 1e-3 / (Kmin * Ninf); //Adaptive timestep for integration, based on characteristic coagulation time.
			if (step >= 1e-6)
				step = 1e-6;
//      step=1e-8;
			fi = 0;
			for (i = 1; i <= MAX - 1; i++) {
				fi = fi + N[i] * v[i];
				if (Nmax < N[i]) {
					Nmax = N[i];
					Vmax = v[i];
				}
			}
			//N_infinity according to equation 7.77 in Friedlander.
			Ninf_eq = Ninf_eq
					- step
							* (3.33 * pow((3 / (4 * pi)), 0.1666666)
									* pow(6 * kb * T / rho, 0.5)
									* pow(fi, 0.1666666)
									* pow(Ninf_eq, (11.0 / 6.0)));
			Vtot = 0;
			printf("\n%e\t%e\t%e\t%e\t%e\t%e", t, Ninf, N[1], N[5], N[10],
					N[15]);
//	  fprintf(fcoag,"\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",t,Ninf,N[1],N[3],N[5],N[8],N[10],N[13],N[15],N[18],N[20],N[23],N[25],N[28],N[30],N[33],N[35],N[38],N[40]);
			t = t + step;
//hong call the subroutine of RKT to calculate the particle temperature and velocity

//	  printf("%e\t%e",Vmax,Nmax);
//	  DIA_P=pow((6*Vmax/pi),0.3333333);	//most probable diameter
			Ntot = 0.0;
			Vtot = 0.0;
			Vgeo = 0.0;
			SIGMAgeo = 0.0; // Initializing total number of particles to zero, as there are no particles prior to nucleation.
			for (i = 2; i <= MAX - 1; i++) {
				Ntot = Ntot + N[i]; //Total number of particles (does not include monomers).
				Vtot = Vtot + N[i] * v[i]; //Total volume of particles. Note that loop runs from i=2 because i=1 corresponds to monomers, which we do not count as particles.
//        DIA_P=DIA_P+pow(pow((6*v[i]/pi),0.3333333),2)*N[i];//RMS diameter
				Vgeo = Vgeo + v[i] * N[i] * log(v[i]);
			}
			Vav = Vtot / Ntot; // Average volume of particles ( number average)
			DIA_P = pow((6.0 * Vav / pi), 0.3333333); //Volume based mean diameter of particles
//      DIA_P=sqrt(DIA_P/Ntot);//RMS diameter
			Vgeo = exp(Vgeo / Vtot);
			for (i = 1; i <= MAX - 1; i++) {
				SIGMAgeo = SIGMAgeo
						+ v[i] * N[i] * pow((log(v[i] / Vgeo)), 2.0);
			}
			SIGMAgeo = exp(sqrt(1.0 / 9.0 * SIGMAgeo / Vtot));
			tRKT = tRKT + step;

			if ((tRKT >= 1e-12) && (DIA_P >= 1e-9)) {
				nn0 = nn0 + 1;
				if (nn0 == 1)
					t0 = t;
				else
					t0 = t0 + step;
//		  printf("\n%e\t%e\t%e",DIA_P,T,Y[2]);
				T = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
				nn1 = nn1 + 1;
				if (nn1 == 1000) {
					fprintf(fparticle, "\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", t,
							t0, Y[1], Y[2], T, DIA_P, step, SIGMAgeo);
					nn1 = 0;
				}
//		  printf("\n%e\t%e\t%e\t%e",t,Y[1],Y[2],T);
//		  tRKT=1e-6;
//		  RKT3(t0,tRKT,POSIF,DENSITYF,VELOF,TEMPF,DV,TC,WM,NSTEP+1);

				RKT3(t0, step, POSIF, DENSITYF, VELOF, TEMPF, DV, TC, WM,
						NSTEP + 1);
//		  tRKT=0.0;
			} else {
				Y[1] = enlgr(POSIF, VELOF, NSTEP + 1, Y[2]);
				T = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
				Y[2] = Y[2] + Y[1] * step;
				nn1 = nn1 + 1;
				if (nn1 == 10) {
					fprintf(fparticle, "\n%e\t%e\t%e\t%e\t%e", t, Y[1], Y[2], T,
							DIA_P);
					nn1 = 0;
				}
			}

		}
		// Printing the final number distribution after SPD is reached.
//      printf("\n\n***Size Distribution after reaching SPD****\n");
//      printf("\nvolume\t\tnumber");
		fprintf(fcoag, "\n\n***Size Distribution after reaching SPD****\n");
		fprintf(fcoag, "\nvolume\t\tnumber");
		for (i = 1; i <= MAX - 1; i++) {
			printf("\n%e\t%e\t%e", v[i], N[i], Sa[i]);
			fprintf(fcoag, "\n%e\t%e\t%e", v[i], N[i], Sa[i]);
		}
		fclose(fptr1);
		fclose(fcoag);
//	  fclose(2)
	}
	//************************************************************
	//             Nucleation + Coagulation
	//************************************************************

	else if (choice == 2) //If user chooses 2, code runs for nucleation + coagulation.
			{
		fcoagnu = fopen("coagnu.dat", "w");
		v[1] = v1; //Assigning first node to monomer. Note that these are not considered as particles.
		for (i = 1; i <= MAX - 1; i++)
			N[i] = 0.0; //Setting initial number of particles = 0 in all the nodes.
		for (i = 2; i <= MAX; i++)
		//Defining volume spacing of bins. Using a geometric factor of 2 to cover 1nm to 10 micrometer size particles.
				{
			v[i] = v[1] * pow(q, (i - 1));
			m[i] = rho * v[i];
			dp[i] = pow(6.0 * v[i] / pi, 0.3333333);
			fprintf(fcoagnu, "\n%e", dp[i]);
		}
		dp[1] = pow(6.0 * v[1] / pi, 0.3333333);
		m[1] = rho * v[1];
		//Calculation of splitting operator for coagulation, X[i][j][k].
		for (k = 1; k <= MAX - 1; k++) {
			for (i = 1; i <= MAX - 1; i++) {
				for (j = 1; j <= MAX - 1; j++) {
					//Conditions in parentheses check if the combined volume of colliding particles is larger than k.
					if (v[k] <= (v[i] + v[j]) && (v[i] + v[j]) < v[k + 1]) {
						X[i][j][k] = (v[k + 1] - v[i] - v[j])
								/ (v[k + 1] - v[k]);
					} else {
						if (v[k - 1] <= (v[i] + v[j]) && (v[i] + v[j]) < v[k]) {
							X[i][j][k] = (v[i] + v[j] - v[k - 1])
									/ (v[k] - v[k - 1]);
						} else
							X[i][j][k] = 0;
					}
				}
			}
		}
		t = 0.0; //Initializing time.
//      step=5e-4;//Timestep1 for integration.
//      step=1e-4; //Timestep2 for integration.
		step = 1e-7; //Timestep3 for integration.
//      step=1e-7; //Timestep4 for integration.
		S = 1.001; // Setting initial saturation ratio, a little larger than 1. Hereafter the saturation ratio is determined by cooling or heating rate of the aerosol. A heating rate would require a negative value.

		// Calculating surface area of the monomer, given the volume of monomer v1.
		temp1 = 6 * v1 / pi;
		s1 = pi * pow(temp1, 0.66666666667);
		//Calculating mass of the monomer unit.
		m1 = rho * v1;
		//Calculating saturation vapor pressure for the Aerosol.
//      Ps=exp(13.07-36373.0/T)*101325.0*1.0;//hong:note that data for Al was used here.
		Ps = exp(18.8836 - 69076.78448 / T) * 101325.0 * 1.0;
		//Calculating number concentration of monomers at saturation.
		ns = Ps / (kb * T);
		// Calculating monomer concentration ( multiplying saturation ratio by saturation concentration of monomers).
		N[1] = S * ns;
		counter = 0; // Temporary variable for printing after a large number of iterations.
		printf("\ntime\t\tMonomer\t\tJk\t\tS\t\tDp_mean\t\tk*\t\tNtot");
		fprintf(fcoagnu,
				"\ntime\t\tMonomer\t\tJk\t\tS\t\tDp_mean\t\tk*\t\tNtot");
//      while(T>300) // Calculating aerosol properties until the system cools down to 27 C.
		while (Y[2] >= 0.00005) {
			sigma = (A - B * T) * 1e-3; // Surface tension.
			Ps = exp(C - D / T) * 101325.0 * 1.0; // Saturation pressure.
			ns = Ps / (kb * T); // Saturation concentration of monomers using ideal gas law.
			S = N[1] / ns; // Saturation ratio.
			theta = (s1 * sigma) / (kb * T); // Dimensionless surface tension theta.
			a = (2 * sigma) / (pi * m1); //Temporary variable.
			b = theta - (4 * pow(theta, 3)) / (27 * pow(log(S), 2)); //Temporary variable.
			Jk = pow(ns, 2) * S * v1 * pow(a, 0.5) * exp(b); //Nucleation rate using the classical SCC model.
			c = 0.6666667 * theta / log(S); // Temporary variable.
			kstar = pow(c, 3.0); //Calculating critical cluster size that will determine the node at which nucleation occurs.
			dpstar = 4 * sigma * v1 / (kb * T * log(S)); //Size of particles corresponding to the critical node.
			vstar = pi * pow(dpstar, 3) / 6; // Volume of particle corresponding to the critical node.
			Ntot = 0.0; // Initializing total number of particles to zero, as there are no particles prior to nucleation.
			for (i = 2; i <= MAX - 1; i++)
				Ntot = Ntot + N[i]; //Total number of particles (does not include monomers).
			Vtot = 0.0;
			for (i = 2; i <= MAX - 1; i++)
				Vtot = Vtot + N[i] * v[i]; //Total volume of particles. Note that loop runs from i=2 because i=1 corresponds to monomers, which we do not count as particles.
			Vtotal = Vtot + N[1] * v[1]; //Total volume for mass conservation check (includes monomers)
			Vav = Vtot / Ntot; // Average volume of particles ( number average)
			dpav = pow((6 * Vav / pi), 0.3333333); //Volume based mean diameter of particles printing after every 50 times the loop runs.
			counter = counter + 1;
			if (counter == 5000) {
				printf("\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", t, N[1], Jk, S,
						dpav, kstar, Ntot, Y[1], Y[2]);
				fprintf(fcoagnu, "\n%e\t%e\t%e\t%e\t%e\t%e\t%e", t, N[1], Jk, S,
						dpav, kstar, Ntot);
				counter = 0;
			}
			//Calculation of collision frequency function beta(i,j)=K(i,j).

			if (beta_option == 1) {
				for (i = 1; i <= MAX - 1; i++) {
					for (j = 1; j <= MAX - 1; j++) {

						temp1 = 1 / v[i] + 1 / v[j];
						temp2 = pow(v[i], 0.333333) + pow(v[j], 0.333333);
						K[i][j] = pow(3.0 / (4.0 * pi), 0.1666666)
								* pow((6.0 * kb * T / rho), 0.5)
								* pow(temp1, 0.5) * pow(temp2, 2.0);
						if (K[i][j] < Kmin)
							Kmin = K[i][j]; //Calculating the smallest collision frequency function to decide the characteristic coagulation time.

					}
				}
			}
			if (beta_option == 2) {
				mu = A_mu * pow(T, 1.5) / (B_mu + T);
				lambda = (mu / (P * 101325.0))
						* sqrt(pi * R * T / (2.0 * 0.04));
				for (i = 1; i <= MAX - 1; i++) {
					for (j = 1; j <= MAX - 1; j++) {
						kn1 = (2.0 * lambda) / dp[i];
						kn2 = (2.0 * lambda) / dp[j];
						D1 = (kb * T) / (3.0 * pi * mu * dp[i])
								* ((5.0 + 4.0 * kn1 + 6.0 * kn1 * kn1
										+ 18.0 * kn1 * kn1 * kn1)
										/ (5.0 - kn1 + (8.0 + pi) * kn1 * kn1));
						D2 = (kb * T) / (3.0 * pi * mu * dp[j])
								* ((5.0 + 4.0 * kn2 + 6.0 * kn2 * kn2
										+ 18.0 * kn2 * kn2 * kn2)
										/ (5.0 - kn2 + (8.0 + pi) * kn2 * kn2));
						c1 = sqrt((8.0 * kb * T) / (pi * m[i]));
						c2 = sqrt((8.0 * kb * T) / (pi * m[j]));
						l1 = (8.0 * D1) / (pi * c1);
						l2 = (8.0 * D2) / (pi * c2);
						g1 = (pow((dp[i] + l1), 3)
								- pow((dp[i] * dp[i] + l1 * l1), 1.5))
								/ (3.0 * dp[i] * l1) - dp[i];
						g2 = (pow((dp[j] + l2), 3)
								- pow((dp[j] * dp[j] + l2 * l2), 1.5))
								/ (3.0 * dp[j] * l2) - dp[j];
						K[i][j] = 2.0 * pi * (D1 + D2) * (dp[i] + dp[j])
								/ ((dp[i] + dp[j])
										/ (dp[i] + dp[j]
												+ 2.0 * sqrt(g1 * g1 + g2 * g2))
										+ (8.0 * (D1 + D2))
												/ (sqrt(c1 * c1 + c2 * c2)
														* (dp[i] + dp[j])));
						if (K[i][j] < Kmin)
							Kmin = K[i][j]; //Calculating the smallest collision frequency function to decide the characteristic coagulation time.
					}
				}
			}

			//Operator to put nucleated particles in the bin just higher than k*.
			for (k = 2; k <= MAX - 1; k++) {
				if (vstar < v1)
					zeta[2] = vstar / v[2]; //Putting particles formed smaller than monomers in the smallest particle node (node 2). This situation arises when k* falls below monomer size.
				else if (v[k - 1] <= vstar && vstar < v[k])
					zeta[k] = vstar / v[k]; // Putting particles in node just larger than k*.
				else
					zeta[k] = 0.0;
			}
			//Calculation of gain and loss terms for Nk due to coagulation.
			for (k = 2; k <= MAX - 1; k++) {
				//Initializing gain and loss terms for coagulation.
				sum1 = 0; // Production term due to collision of two smaller particles to form a k size particle.
				sum2 = 0; // Loss term due to collision of k size particle with any other particle.
				for (i = 2; i <= MAX - 1; i++) {
					sum2 = sum2 + K[k][i] * N[i]; // k collides with any other particle to get out of node k, thus loss term.
					for (j = 2; j <= k; j++)
						sum1 = sum1 + X[i][j][k] * K[i][j] * N[i] * N[j]; // i and j collide to form particle of size k which is then multiplied by the size splitting operator to put the particles in the adjacent nodes after adjusting the volume.
				}
				// Change in PSD due to nucleation + coagulation.
				N[k] = N[k] + step * (0.5 * sum1 - N[k] * sum2 + Jk * zeta[k]);
			}
			//Monomer balance. accounting for the monomers loss due to nucleation.
			N[1] = N[1] - Jk * kstar * step;
			DIA_P = dpav;
			tRKT = tRKT + step;
			if ((tRKT >= 1e-12) && (DIA_P >= 1e-9)) {
				nn0 = nn0 + 1;
				if (nn0 == 1)
					t0 = t;
				else
					t0 = t0 + step;
//		  printf("\n%e\t%e\t%e",DIA_P,T,Y[2]);
				T = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
				fprintf(fcoagnu, "\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", t, t0,
						Y[1], Y[2], T, DIA_P, step);
				RKT3(t0, step, POSIF, DENSITYF, VELOF, TEMPF, DV, TC, WM,
						NSTEP + 1);
			} else {
				Y[1] = enlgr(POSIF, VELOF, NSTEP + 1, Y[2]);
				T = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
				Y[2] = Y[2] + Y[1] * step;
				nn1 = nn1 + 1;
				if (nn1 == 10) {
					fprintf(fcoagnu, "\n%e\t%e\t%e\t%e\t%e", t, Y[1], Y[2], T,
							DIA_P);
					nn1 = 0;
				}
			}

//	  T=T-step*coolrate;// Temperature would decrease as time progresses, due to cooling.
			t = t + step; //Time increment.
		}
		printf("\n\nFinal Particle Size Distribution\n\n");
		printf("\nvolume\t\tnumber");
		fprintf(fcoagnu, "\n\nFinal Particle Size Distribution\n\n");
		fprintf(fcoagnu, "\nvolume\t\tnumber");
		for (i = 2; i <= MAX - 1; i++) {
			printf("\n%e\t%e", v[i], N[i]); // Printing the final number distribution obtained due to nucleation and coagulation. (number concentration or particles against nodes).
			fprintf(fcoagnu, "\n%e\t%e", v[i], N[i]); // Printing the final number distribution obtained due to nucleation and coagulation. (number concentration or particles against nodes).
		}
		fclose(fcoagnu);
	}

	//*****************************************
	//SURFACE GROWTH + NUCLEATION + COAGULATION
	//*****************************************
	else if (choice == 3) // THE COMPLETE GDE
			{
		fconusurf = fopen("conusurf.dat", "w");
		fextra2 = fopen("extra2.dat", "w");
		T_FLAME = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);

		for (i = 1; i <= MAX - 1; i++) {
			N[i] = 0.0; // Setting number of particles = 0 at all the nodes initially.
			N1s[i] = 0.0; // Setting number concentration of monomers over different size particles = 0 initially.
			SURF_MONO[i] = 0.0;
			SA[i] = 0.0;
			VV[i] = 0.0;

		}
		v[1] = v1; // Setting volume of the first node equal to the monomer volume, i.e. volume of a molecule, which is v1.
		dp[1] = pow((6.0 * v[1] / pi), 0.333333);
		m[1] = rho * v[1];
		Sa[1] = pi * dp[1] * dp[1];
		SA[1] = N[1] * Sa[1];
		VV[1] = N[1] * v[1];
		for (i = 2; i <= MAX; i++) {
			v[i] = v[1] * pow(q, (i - 1)); // Calculating volume of all larger nodes using a geometric factor of 2. Note that volume of node 1 is the volume of a molecule.
			dp[i] = pow((6.0 * v[i] / pi), 0.333333);
			fprintf(fconusurf, "\n%i\t%e", i, dp[i]);
			m[i] = rho * v[i];
			Sa[i] = Sa[1] * pow(q, (i - 1));

		}

		//Size splitting operator for coagulation ( Same as in "coagulation + nucleation" section of the code).
		for (k = 1; k <= MAX - 1; k++) {
			for (i = 1; i <= MAX - 1; i++) {
				for (j = 1; j <= MAX - 1; j++) {
					if (v[k] <= (v[i] + v[j]) && (v[i] + v[j]) < v[k + 1]) {
						X[i][j][k] = (v[k + 1] - v[i] - v[j])
								/ (v[k + 1] - v[k]);
					} else {
						if (v[k - 1] <= (v[i] + v[j]) && (v[i] + v[j]) < v[k]) {
							X[i][j][k] = (v[i] + v[j] - v[k - 1])
									/ (v[k] - v[k - 1]);
						} else
							X[i][j][k] = 0;
					}
				}
			}
		}
		t = 0.0; // Initializing time.
//      step=1e-4; // Initial time step for integration.
		step = 5e-8;
		S = 1.001; // Initial saturation ratio. Assuming that the vapor is just saturated initially.
		temp1 = 6 * v1 / pi; // Temporary variable.
		s1 = pi * pow(temp1, 0.66666666667); // Surface area of monomer.
		m1 = rho * v1; // Mass of a monomer.
		Ps = exp(C - D / T) * 101325.0 * 1.0; // Saturation vapor pressure at temperature T.
		ns = Ps / (kb * T); // Saturation concentration of monomers using ideal gas law.
		N[1] = S * ns; // Monomer concentration.
//      printf("\nTime\t\tN[1]\t\tJk\t\tS\t\tdiameter\tParticle Volume\tTotal Number");// Printing result heading.
		fprintf(fconusurf,
				"\nTime\t\tN[1]\t\tJk\t\tS\t\tdiameter\tParticle Volume\tTotal Number"); // Printing result heading.

		// Counters below are temporary variables that determine number of time the loop runs before every print step. Generally time steps are small and printing results is done after every 4000 steps.
		counter = 1;
		counter2 = 1;
//      while(T>300) // Observing changes until temperature drops down to room temperature ~ 27 C.
		while (Y[2] >= 0.00005) {
			// Terms below have been explained in "nucleation + coagulation" section of the code above.
			if (Y[2] < 0.028) {
				EF = -700.0 / 0.04;
			}
			T_FLAME = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
			sigma = (A - B * T) * 1e-3;
			Ps = exp(C - D / T) * 101325.0 * 1.0;
			ns = Ps / (kb * T);
			S = N[1] / ns;
			theta = (s1 * sigma) / (kb * T);
			a = (2 * sigma) / (pi * m1);
			b = theta - (4.0 * pow(theta, 3)) / (27.0 * pow(log(S), 2));
			Jk = pow(ns, 2) * S * v1 * pow(a, 0.5) * exp(b);
//	  if(t>=9.4000e-005) printf("\n%e\t%e",Jk,b);
			c = 0.6666667 * theta / log(S);
			kstar = pow(c, 3.0);
			dpstar = 4 * sigma * v1 / (kb * T * log(S));
			vstar = pi * pow(dpstar, 3) / 6;
//	  RHO_FLAME=enlgr(POSIF,DENSITYF,NSTEP+1,Y[2]);
			FLAG_RK = 0;
			nn_count = nn_count + 1;
			Ntot = 0.0;
			te_mono = 0.0;
			te_mono1 = 0.0;
			fs_tio2A[1] = 7.4375e16 * T_FLAME * pow(dp[1], 4)
					* exp(3.103e4 / T_FLAME);
			for (i = 2; i <= MAX - 1; i++) {
				Ntot = Ntot + N[i]; //Total number of particles.
				te_mono = te_mono + N[i] * pi * pow(dp[i], 2);
				te_mono1 = te_mono1 + SA[i];
				fs_tio2A[i] = 7.4375e16 * T_FLAME * pow(dp[i], 4)
						* exp(3.103e4 / T_FLAME);

//	  fprintf(fconusurf,"\n%e\t%e\t%e\t%e",te_mono,N[i],RHO_FLAME,SURF_MONO[i]);
			}

			Vtot = 0.0;
			Vtotal = 0.0;
			Stot = 0.0;
			Vgeo = 0.0;
			SIGMAgeo = 0.0;
			for (i = 2; i <= MAX - 1; i++) {
				Vtot = Vtot + N[i] * v[i]; //Total volume of particles . This volume is used to evaluate the mean volume and thus mean diameter of particles.
				Stot = Stot + N[i] * Sa[i]; //total surface area os particles.
				Vgeo = Vgeo + N[i] * log(v[i]);
			}
			Vtotal = Vtot + N[1] * v[1]; //Total mass volume (including monomers) This volume is used to check for mass conservation.
			Vav = Vtot / Ntot; // Mean volume of particles.
			dpav = pow((6 * Vav / pi), 0.3333333); //Mean diameter of particles.
			dpap = 6.0 * Vtot / Stot; //mean primary diameter
			Vgeo = exp(Vgeo / Ntot);
			for (i = 2; i <= MAX - 1; i++) {
				SIGMAgeo = SIGMAgeo + N[i] * pow((log(v[i] / Vgeo)), 2.0);
			}
			SIGMAgeo = exp(sqrt(1.0 / 9.0 * SIGMAgeo / Ntot)); //geometric standard deviation based on number weighted rather than volume weighted

//	  if(Ntot>100.0)//Reducing the timestep for integration when there are sufficient number of particles such that surface growth takes over and the code can no longer run for large time-steps.
//	    step=1e-5;

			counter = counter + 1;
			counter2 = counter2 + 1;
			if (counter == 100) // Printing mean diameter etc. after every 100 steps.
					{
				printf("\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", t, N[1], Jk, S,
						dpav, Vtot, Ntot, Y[1], Y[2]);
				fprintf(fextra2, "\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",
						t, N[1], Jk, S, dpav, dpap, SIGMAgeo, Vtot, Ntot, Y[1],
						Y[2]);
				counter = 1;
			}
			if (counter2 == 40000) // Printing PSD after every 40,000 steps.
					{
//	      for (i=1;i<=MAX-1;i++)
//		  {
//		  printf("\nN%d\t%e",i,N[i]);
//		  fprintf(fconusurf,"\nN%d\t%e",i,N[i]);}
				counter2 = 1;
			}
			//Collision frequency function for coagulation.
			if (beta_option == 1) {
				for (i = 1; i <= MAX - 1; i++) {
					for (j = 1; j <= MAX - 1; j++) {

						temp1 = 1 / v[i] + 1 / v[j];
						temp2 = pow(v[i], 1.0 / Df) + pow(v[j], 1.0 / Df);
						K[i][j] = pow(3.0 / (4.0 * pi), (2.0 / Df - 0.5))
								* pow((6.0 * kb * T / rho), 0.5)
								* pow(temp1, 0.5) * pow(temp2, 2.0);
						if (K[i][j] < Kmin)
							Kmin = K[i][j]; //Calculating the smallest collision frequency function to decide the characteristic coagulation time.
// 	     printf("\n%e",Kmin);
					}
				}
			}
			if (beta_option == 2) {
				mu = A_mu * pow(T, 1.5) / (B_mu + T);
				lambda = (mu / (P * 101325.0))
						* sqrt(pi * R * T / (2.0 * 0.04));
//		  printf("\n%e",lambda);
				for (i = 1; i <= MAX - 1; i++) {
					for (j = 1; j <= MAX - 1; j++) {
						kn1 = (2.0 * lambda) / dp[i];
						kn2 = (2.0 * lambda) / dp[j];
						D1 = (kb * T) / (3.0 * pi * mu * dp[i])
								* ((5.0 + 4.0 * kn1 + 6.0 * kn1 * kn1
										+ 18.0 * kn1 * kn1 * kn1)
										/ (5.0 - kn1 + (8.0 + pi) * kn1 * kn1));
						D2 = (kb * T) / (3.0 * pi * mu * dp[j])
								* ((5.0 + 4.0 * kn2 + 6.0 * kn2 * kn2
										+ 18.0 * kn2 * kn2 * kn2)
										/ (5.0 - kn2 + (8.0 + pi) * kn2 * kn2));
						c1 = sqrt((8.0 * kb * T) / (pi * m[i]));
						c2 = sqrt((8.0 * kb * T) / (pi * m[j]));
						l1 = (8.0 * D1) / (pi * c1);
						l2 = (8.0 * D2) / (pi * c2);
						g1 = (pow((dp[i] + l1), 3)
								- pow((dp[i] * dp[i] + l1 * l1), 1.5))
								/ (3.0 * dp[i] * l1) - dp[i];
						g2 = (pow((dp[j] + l2), 3)
								- pow((dp[j] * dp[j] + l2 * l2), 1.5))
								/ (3.0 * dp[j] * l2) - dp[j];
						K[i][j] = 2.0 * pi * (D1 + D2) * (dp[i] + dp[j])
								/ ((dp[i] + dp[j])
										/ (dp[i] + dp[j]
												+ 2.0 * sqrt(g1 * g1 + g2 * g2))
										+ (8.0 * (D1 + D2))
												/ (sqrt(c1 * c1 + c2 * c2)
														* (dp[i] + dp[j])));
						if (K[i][j] < Kmin)
							Kmin = K[i][j]; //Calculating the smallest collision frequency function to decide the characteristic coagulation time.
// 		     printf("\n%e",Kmin);
					}
				}
			}
			//Putting nucleated particles in bins just larger than k* ( explained in section " coagulation + nucleation" above).
			for (k = 2; k <= MAX - 1; k++) {
				if (vstar < v1)
					zeta[2] = vstar / v[2];
				else if (v[k - 1] <= vstar && vstar < v[k])
					zeta[k] = vstar / v[k];
				else
					zeta[k] = 0.0;
			}
			//Calculating the saturation monomer concentration over an i sized particle that has a diameter dp[i].
			for (i = 1; i <= MAX - 1; i++) {
				dp[i] = pow((6 * v[i] / pi), 0.333333); // Diameter of particle at the ith node.
				N1s[i] = ns * exp(4 * sigma * MW / (R * T * rho * dp[i])); // Calculating saturation monomer concentration over a i size particle including the kelvin effect.
			}
			// t1 through t4 are terms used for monomer balance. These are loss and gain terms due to condensation and evaporation respectively. Each of these terms has been explained below as I use them. Initially setting these terms = 0.
			t1 = 0;
			t2 = 0;
			t3 = 0;
			t4 = 0;
			// Note that in the following part "addterm" and "subterm" are addition and subtraction term in particle size distribution due to condensation or evaporation (surface growth), where as terms t1 through t4 are addition and subtraction term for monomer balance.
			for (k = 2; k <= MAX - 1; k++) {
				sum1 = 0; //Addition term due to coagulation of size i+j to form k.
				sum2 = 0; //Loss term due to coagulation of k with any other particle.
				sum3 = 0.0;
				sum22 = 0.0;
				sum4 = 0.0;
				if (step != 1e-4) //Activating surface growth when sufficient number of particles have nucleated. At smaller times there are very small number of particles and surface growth hardly changes PSD. Doing so makes the code run faster.
						{
					// The following four "if" statements check if the vapor pressure of monomer in the system is larger or smaller than saturation vapor pressure over the k size particle. The difference in the two vapor pressures is the driving force for evaporation or condensation.
					if (N[1] > N1s[k - 1]) // i.e. if actual monomer concentration in the aerosol is greater than saturation concentration over a k-1 size particle.
							{
						if (k == 2)
							addterm = 0.0; // There are no particles smaller than 2nd node. so particles in 2nd node can not be added due to growth of smaller particles.
						else {
							addterm = (v[1] / (v[k] - v[k - 1])) * K[1][k - 1]
									* (N[1] - N1s[k - 1]) * N[k - 1]; //Growth of k due to condensation of monomers on k-1.
							t1 = t1
									+ K[1][k - 1] * (N[1] - N1s[k - 1])
											* N[k - 1]; //Loss of monomers that have condensed.
						}
					}
					// Similarly, we can explain the following terms below.
					if (N[1] < N1s[k + 1]) {
						addterm = -(v[1] / (v[k + 1] - v[k])) * K[1][k + 1]
								* (N[1] - N1s[k + 1]) * N[k + 1]; //Growth of k due to evaporation of monomers from k+1.
						t2 = t2 + K[1][k + 1] * (-N[1] + N1s[k + 1]) * N[k + 1]; //Gain of monomers that have evaporated.
					}
					if (N[1] < N1s[k]) {
						subterm = -(v[1] / (v[k] - v[k - 1])) * K[1][k]
								* (N[1] - N1s[k]) * N[k]; //Loss of k due to evaporation of monomers from k.
						t3 = t3 + K[1][k] * (-N[1] + N1s[k]) * N[k]; //Gain of monomers that have evaporated.
					}
					if (N[1] > N1s[k]) {
						subterm = (v[1] / (v[k + 1] - v[k])) * K[1][k]
								* (N[1] - N1s[k]) * N[k]; //Loss of k due to condensation of monomers on k.
						t4 = t4 + K[1][k] * (N[1] - N1s[k]) * N[k]; //Loss of monomers that have condensed.
					}
				}
				// Terms due to coagulation.
				for (i = 2; i <= MAX - 1; i++) {
					sum2 = sum2 + K[k][i] * N[i];
					sum22 = sum22 + K[k][i] * N[i];
					for (j = 2; j <= k; j++) {
						sum1 = sum1 + X[i][j][k] * K[i][j] * N[i] * N[j];
						sum3 = sum3
								+ X[i][j][k] * K[i][j] * N[i] * N[j] * Sa[j]
										* pow(q, (k - j)) * (Sa[i] + Sa[j])
										/ Sa[k]
										/ (pow(q, (i - k)) + pow(q, (j - k)));
					}
				}
//		  printf("\n\n%i\t%e\t%e\t%e",k,Sa[k]*N[k],SA[k],6.0*VV[k]/SA[k]);
//		  ppd=6.0*v[k]/Sa[k];
//		  fs_tio2=7.4375e16*pow(ppd,4)*T_FLAME*exp(31030.0/T_FLAME);
				temNp = pow(Sa[k], 3) / 36.0 / pi / v[k] / v[k];
				if (temNp <= 1.0)
					sum4 = 0.0;
				else
					sum4 = 1.0
							* (SA[k]
									- N[k] * pi
											* pow((6.0 * v[k]) / pi,
													0.666666667)) / fs_tio2A[k];

				if (sum4 < 0)
					sum4 = 0;
//        sum4=0;
				SA[k] = SA[k] + step * (0.5 * sum3 - Sa[k] * N[k] * sum22);
				SA[k] = SA[k] - step * sum4;
//        VV[k]=VV[k]+step*(0.5*sum1-N[k]*sum2)*v[k];

//		  N[k]= N[k] + step*(0.5*sum1 - N[k]*sum2);

// Change in particle size due to nucleation, coagulation and surface growth.
				N[k] = N[k]
						+ step
								* (0.5 * sum1 - N[k] * sum2 + Jk * zeta[k]
										+ addterm - subterm); //PSD when all the three phenomena are activated.
				SA[k] = SA[k]
						+ step * Sa[k] * (Jk * zeta[k] + addterm - subterm);
// change in particle number concentration due to precursor pyrolysis and monomer attachment to particles
				if (CNT_SUM <= 0.99) {
					T_FLAME = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
					V_FLAME = enlgr(POSIF, VELOF, NSTEP + 1, Y[2]);
					A_vk = N[k] * pi * pow(6.0 * v[k] / pi, 0.66666667);
//	  RHO_FLAME=enlgr(POSIF,DENSITYF,NSTEP+1,Y[2]);
					RK = 3.96e5 * exp(-8479.7 / T_FLAME) * 10.0; // reaction rate of precursor decomposition
					RKs = 1.0e9 * exp(-15155.16 / T_FLAME);
					RKg = RK - te_mono * RKs;

					if (RKg <= 0.0) {
						RKg = 0.0;
						FLAG_RK = 1;
						RKs = RK / te_mono;
					}
//	  fprintf(fconusurf,"\n%e\t%e\t%e",RK,RKg,RKs);

					time_mono = -(YY_TEMP - Y[2]) / V_FLAME;
					SURF_MONO[k] = RKs * CNT_MONO1 * time_mono * A_vk;
					CNT_MONO2 = CNT_MONO1 - SURF_MONO[k];

					CNT_SUM = CNT_SUM + SURF_MONO[k] / CNT_MONO0;
					CNT_MONO1 = CNT_MONO2;
					N[k] = N[k] + v[1] * SURF_MONO[k - 1] / (v[k] - v[k - 1])
							- v[1] * SURF_MONO[k] / (v[k + 1] - v[k]);
					SA[k] = SA[k]
							+ Sa[k]
									* (v[1] * SURF_MONO[k - 1]
											/ (v[k] - v[k - 1])
											- v[1] * SURF_MONO[k]
													/ (v[k + 1] - v[k]));
//      printf("\n%i\t%e\t%e\t%e",k,N[k],A_vk,SURF_MONO[k]);
				}
				VV[k] = N[k] * v[k];

				if (k > 1) {
					if (SA[k] < N[k] * pi * pow(6.0 * v[k] / pi, 0.66666667)) {
						SA[k] = N[k] * pi * dp[k] * dp[k];
					}
					if (SA[k] > N[k] * Sa[1] * pow(q, (k - 1))) {
						SA[k] = N[k] * Sa[1] * pow(q, (k - 1));
					}
//		if (SA[k]<N[k]*Sa[k-1]*pow(q,0.6666666667))
//		{ SA[k]=N[k]*Sa[k-1]*pow(q,0.6666666667);}
//		if (SA[k]>N[k]*Sa[k-1]*q)
//		{ SA[k]=N[k]*Sa[k-1]*q;}
				}
//         printf("\n%i\t%e\t%e\t%e",k,6.0*VV[k]/SA[k],N[k],Y[2]);
			}

			for (i = 2; i <= MAX - 1; i++) {
				if (N[i] != 0)
					Sa[i] = SA[i] / N[i];
				coeff1 = 1.0;
				coeff = 1.0;
				mstep = 0;
				if (N[i] > 0.0) {
					temNp = pow(Sa[i], 3) / 36.0 / pi / v[i] / v[i];
					ppd = 6.0 * v[i] / Sa[i];
					if (temNp > 3.0) {
						coalstep = log(temNp) / log(3.0);
						mstep = (int) (log(temNp) / log(3.0)); // convert to int type.

						for (jj = 1; jj <= mstep - 1; jj++)
							coeff1 = coeff1 + pow(3.0, 4.0 / 3.0 * jj);

						if (coalstep > mstep)
							coeff = coeff1
									+ pow(3.0, 4.0 / 3.0 * mstep)
											* (coalstep - mstep); // convert to double type
						else
							coeff = coeff1
									+ pow(3.0, 4.0 / 3.0 * (mstep - 1))
											* (coalstep - mstep); // convert to double type

					}
					fs_tio2A[i] = coeff * 7.4375e16 * T_FLAME * pow(ppd, 4)
							* exp(3.103e4 / T_FLAME);
				}
//       printf("\n%i\t%i\t%e\t%e\t%e",i,mstep,coeff1,coeff,fs_tio2A[i]);
			}
			if (nn_count == 1000) {
				for (i = 1; i <= MAX - 1; i++)
					fprintf(fextra1, "\n%i\t%e\t%e\t%e\t%e\t%e\t%e", i, dp[i],
							N[i], SA[i], 6.0 * v[i] / Sa[i], Sa[i + 1] / Sa[i],
							Y[2]);
				nn_count = 0;
			}
//hong--monomer N[1] concentration change because of precursor pyrolyze----
			if (CNT_SUM <= 0.99) {
//		  T_FLAME=enlgr(POSIF,TEMPF,NSTEP+1,Y[2]);
//      V_FLAME=enlgr(POSIF,VELOF,NSTEP+1,Y[2]);
//	  RK=3.96e5*exp(-8479.7/T_FLAME)*10.0;
//----RKg was calculated above---------------
				time_mono = -(YY_TEMP - Y[2]) / V_FLAME;
				N[1] = N[1] + RKg * CNT_MONO1 * time_mono;
				CNT_MONO2 = CNT_MONO1 * (1.0 - RKg * time_mono);

				CNT_SUM = CNT_SUM + RKg * CNT_MONO1 * time_mono / CNT_MONO0;
				CNT_MONO1 = CNT_MONO2;
				YY_TEMP = Y[2];
//	  fprintf(fconusurf,"\n%e\t%e\t%e\t%e\t%i",t,CNT_SUM,CNT_MONO1,Y[2],FLAG_RK);
//	  printf("\n%e\t%e\t%e\t%e\t%e\t%i",t,CNT_SUM,CNT_MONO1,Y[2],te_mono*RKs,FLAG_RK);

			}
			// Change in the concentration of monomers as they are used up in nucleation and surface growth.
			N[1] = N[1] - Jk * kstar * step - step * (t1 + t4 - t3 - t2); //Monomer balance equation.
			SA[1] = Sa[1] * N[1];
			VV[1] = v[1] * N[1];
//	  T=T-step*coolrate;//Reducing the temperature at the given cooling rate.
			DIA_P = dpav;
			tRKT = tRKT + step;
			if ((tRKT >= 1e-12) && (DIA_P >= 1e-10)) {
				nn0 = nn0 + 1;
				if (nn0 == 1)
					t0 = t;
				else
					t0 = t0 + step;
//		  printf("\n%e\t%e\t%e",DIA_P,T,Y[2]);
				T = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
//		  printf("\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",t,t0,Y[1],Y[2],T,DIA_P,step);
				RKT3(t0, step, POSIF, DENSITYF, VELOF, TEMPF, DV, TC, WM,
						NSTEP + 1);
			} else {
				Y[1] = enlgr(POSIF, VELOF, NSTEP + 1, Y[2]);
				T = enlgr(POSIF, TEMPF, NSTEP + 1, Y[2]);
				Y[2] = Y[2] + Y[1] * step;
//		  nn1=nn1+1;
//		  if(nn1==10) {
//			  printf("\n%e\t%e\t%e\t%e\t%e",t,Y[1],Y[2],T,DIA_P); nn1=0;}
			}

			t = t + step; // Time increment.
		}
		for (i = 1; i <= MAX - 1; i++) {
			fprintf(fconusurf, "\nN%d\t%e", i, N[i]);
		}
		fprintf(fextra2, "\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", t,
				N[1], Jk, S, dpav, dpap, SIGMAgeo, Vtot, Ntot, Y[1], Y[2]);
		for (i = 1; i <= MAX - 1; i++) {
//		  printf("\n%e\t%e\t%e",v[i],N[i],Sa[i]);
			fprintf(fextra1, "\n%i\t%e\t%e\t%e\t%e", i, dp[i], N[i], SA[i],
					6.0 * v[i] / Sa[i]);
		}

		fclose(fconusurf);
		fclose(fextra1);
		fclose(fextra2);
	}

	//*****************************************
	//         PURE SURFACE GROWTH
	//*****************************************

	else if (choice == 4) {
		fsurf = fopen("surf.dat", "w");
		for (i = 1; i <= MAX - 1; i++) {
			N[i] = 0.0;
			N1s[i] = 0.0;
		}
		v[1] = v1;
		for (i = 2; i <= MAX; i++)
			v[i] = v[1] * pow(q, (i - 1));
		if ((fptr2 = fopen("grow.inp", "r+")) == NULL) //Opening the input file for surface growth problems
		{
			printf("Warning: Could not open file");
			return (-1);
		}
		//Scanning the initial size and number concentration of particles whose surface growth is to be studied.
		fscanf(fptr2,
				"Enter the node in which you want to put in particles:%d Enter the particle number concentration in #/m3 initially:%lf",
				&i, &temp);
		fclose(fptr2);
		N[i] = temp;
		t = 0.0;
		step = 1e-8;
		S = 1.001;
		temp1 = 6 * v1 / pi;
		s1 = pi * pow(temp1, 0.66666666667);
		m1 = rho * v1;
		Ps = exp(C - D / T) * 101325.0 * 1.0;
		ns = Ps / (kb * T);
		N[1] = S * ns;
		counter = 0;
		while (T > 300) {
			sigma = (A - B * T) * 1e-3;
			Ps = exp(C - D / T) * 101325.0 * 1.0;
			ns = Ps / (kb * T);
			S = N[1] / ns;
			theta = (s1 * sigma) / (kb * T);
			a = (2 * sigma) / (pi * m1);
			b = theta - (4 * pow(theta, 3)) / (27 * pow(log(S), 2));
			Jk = pow(ns, 2) * S * v1 * pow(a, 0.5) * exp(b);
			c = 0.6666667 * theta / log(S);
			kstar = pow(c, 3.0);
			if (counter == 1000)  //printing after every 1000 steps
					{
				printf("\n\ntime=%es\n", t);
				fprintf(fsurf, "\n\ntime=%es\n", t);
				for (i = 1; i <= MAX - 2; i++) {
					printf("\nN[%d]\t%e", i, N[i]);
					fprintf(fsurf, "\nN[%d]\t%e", i, N[i]);
				}
				counter = 0;
			}
			counter = counter + 1;
			//collision frequency functions for the monomers
			for (j = 1; j <= MAX - 1; j++) {
				temp1 = 1 / v[1] + 1 / v[j];
				temp2 = pow(v[1], 0.333333) + pow(v[j], 0.333333);
				K[1][j] = pow(3.0 / (4.0 * pi), 0.1666666)
						* pow((6.0 * kb * T / rho), 0.5) * pow(temp1, 0.5)
						* pow(temp2, 2.0);
			}
			//calculation of saturation concentration of monomers over an i sized particle
			for (i = 1; i <= MAX - 1; i++) {
				dp[i] = pow((6 * v[i] / pi), 0.333333);
				N1s[i] = ns * exp(4 * sigma * MW / (R * T * rho * dp[i]));
			}
			for (k = 2; k <= MAX - 2; k++) {
				if (step != 1e-4) {
					if (N[1] > N1s[k - 1]) {
						if (k == 2)
							addterm = 0.0;
						else
							addterm = (v[1] / (v[k] - v[k - 1])) * K[1][k - 1]
									* (N[1] - N1s[k - 1]) * N[k - 1]; //growth of k due to condensation of monomers on k-1
					}
					if (N[1] < N1s[k + 1])
						addterm = -(v[1] / (v[k + 1] - v[k])) * K[1][k + 1]
								* (N[1] - N1s[k + 1]) * N[k + 1]; //growth of k due to evaporation of monomers from k+1
					if (N[1] < N1s[k])
						subterm = -(v[1] / (v[k] - v[k - 1])) * K[1][k]
								* (N[1] - N1s[k]) * N[k]; //loss of k due to evaporation of monomers from k
					if (N[1] > N1s[k])
						subterm = (v[1] / (v[k + 1] - v[k])) * K[1][k]
								* (N[1] - N1s[k]) * N[k]; //loss of k due to condensation of monomers on k
				}
				N[k] = N[k] + step * (addterm - subterm); //changing PSD due to surface growth/evaporation only
			}
			T = T - step * coolrate;
			t = t + step;
		}
		fclose(fsurf);
	}

	fclose(fparticle);
}

