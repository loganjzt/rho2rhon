#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <malloc.h>

#define PI 3.14159265358979323846 /* pi */

    using namespace std;

int main(int argc, char **argv){
	
	/*variables*/
	FILE *data;	
	char f[40];
	char fTmp[100];

	int nBin;
	double *rho1;
	double *rho2;
	double *z;
	double zTmp,rhoTmp;

	double xBox = 5.0;
	double yBox = 5.0;
	double zBox = 30.0;

	int natoms = 4500;

	double rho_max, rho_min;
	double a_1, a_2, b;		// parameter for tanh fn fitting one slab
	double a1_1, a1_2, a2_1, a2_2;		// parameter for tanh fn two slabs
	double Atmp,Ctmp;
	
	/* cal weight */
	int wnMax = 20;
	double weight1[wnMax],weight2[wnMax];
	double weight1_im[wnMax],weight2_im[wnMax];
	double weight1_max,weight2_max;
	double cosTmp1,cosTmp2;
	double sinTmp1,sinTmp2;
	double x1;

	a_1 = 8.09041;		
	a_2 = zBox - a_1;

	b   = 0.863882;

	rho_max = 11.592851;
	rho_min = 1.218462;

	/* read rho(z) data from file */
	sprintf(f,"rho_z.dat");
	data = fopen(f,"r");
	fgets(fTmp,100,data);
	nBin = 0;
	while(!feof(data)){
		fgets(fTmp,100,data);
		nBin ++;
	}	

	nBin = nBin - 5;
	fclose(data);

	printf("# number of bins in Z = %d\n",nBin);

	rho1 =(double *)malloc(sizeof(double)*nBin);
	rho2 =(double *)malloc(sizeof(double)*nBin);
	z =(double *)malloc(sizeof(double)*nBin);
	for(int i=0;i<nBin;i++){ 
		rho1[i] = 0.0;
		rho2[i] = 0.0;
		z[i] = 0.0;
	}
	
	data = fopen(f,"r");
	fgets(fTmp,100,data);
	for(int i =0;i<nBin;i++){
		fscanf(data,"%lf %lf %lf \n",&zTmp,&rhoTmp,&x1);
		z[i] = zTmp;
		rho1[i] = rhoTmp;
	}	
	fclose(data);
	printf("# Finish reading rho data\n");

	/* get a_1 a_2 */
	Atmp = sqrt(sqrt(  cosh( (zBox/2.0-a_1)/b ) * cosh( (zBox/2.0-a_2)/b ) / cosh(-a_1/b) / cosh( (zBox-a_2) / b ) ) );
	Ctmp = b/2.0 * log( (Atmp*exp(zBox/4.0/b)-1.0)/(1.0-Atmp*exp(-zBox/4.0/b)) );

	printf("# Atmp = %.4f\n",Atmp);
	printf("# Ctmp = %.4f\n",Ctmp);

	a1_1 = zBox/4.0 - Ctmp;
	a1_2 = zBox/4.0 + Ctmp;
	a2_1 = a1_1 + zBox/2.0;
	a2_2 = a1_2 + zBox/2.0;

	/* using f(x) instead of using rho from simulation */

	for(int i=0;i<nBin;i++){
		if(z[i] < zBox / 2.0)				rho1[i] = (rho_max - rho_min)/2.0*( tanh ( (z[i]-a_1)/b ) ) + (rho_max + rho_min)/2.0;
		else if(z[i] <= zBox )				rho1[i] = (rho_max - rho_min)/2.0*( tanh ( (-z[i]+a_2)/b ) ) + (rho_max + rho_min)/2.0;
	}

	for(int i=0;i<nBin;i++){
		if(z[i] < zBox / 4.0)				rho2[i] = (rho_max - rho_min)/2.0*( tanh ( (z[i]-a1_1)/b ) ) + (rho_max + rho_min)/2.0;
		else if(z[i] < zBox / 2.0)			rho2[i] = (rho_max - rho_min)/2.0*( tanh ( (-z[i]+a1_2)/b ) ) + (rho_max + rho_min)/2.0;
		else if(z[i] < zBox / 4.0 * 3.0)	rho2[i] = (rho_max - rho_min)/2.0*( tanh ( (z[i]-a2_1)/b ) ) + (rho_max + rho_min)/2.0;
		else if(z[i] <= zBox)				rho2[i] = (rho_max - rho_min)/2.0*( tanh ( (-z[i]+a2_2)/b ) ) + (rho_max + rho_min)/2.0;
	}

	/* get weight */
	
	for(int wn = 1; wn <= wnMax; wn ++){
		weight1[wn] = 0.0;
		weight1_im[wn] = 0.0;
		weight2[wn] = 0.0;
		weight2_im[wn] = 0.0;
		
		for(int i=0;i<nBin-1;i++){
			cosTmp1 = cos(2.0*PI*wn/zBox*z[i]);
			cosTmp2 = cos(2.0*PI*wn/zBox*z[i+1]);
			
			weight1[wn] += (rho1[i]*cosTmp1 + rho1[i+1]*cosTmp2 ) * (z[i+1]-z[i]) *xBox*yBox/2.0 ;
			weight2[wn] += (rho2[i]*cosTmp1 + rho2[i+1]*cosTmp2 ) * (z[i+1]-z[i])  *xBox*yBox/2.0;

			sinTmp1 = sin(2.0*PI*wn/zBox*z[i]);
			sinTmp2 = sin(2.0*PI*wn/zBox*z[i+1]);

			weight1_im[wn] += (rho1[i]*sinTmp1 + rho1[i+1]*sinTmp2 ) * (z[i+1]-z[i]) *xBox*yBox/2.0;
			weight2_im[wn] += (rho2[i]*sinTmp1 + rho2[i+1]*sinTmp2 ) * (z[i+1]-z[i]) *xBox*yBox/2.0;

		}
		if(wn == 1)	weight1_max = abs(weight1[wn]);
		if(wn == 2)	weight2_max = abs(weight2[wn]);
	}

	sprintf(f,"weight.dat");
	data = fopen(f,"w");

	fprintf(data,"# wn\tweight1Re\tweight1Re/ReabsMax\tweight1im\tweight1im/ReabsMax\n");
	for(int wn=1;wn<=wnMax;wn++) fprintf(data,"%d\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\n",wn,weight1[wn],weight1[wn]/weight1_max,weight1_im[wn],weight1_im[wn]/weight1_max,weight2[wn],weight2[wn]/weight2_max,weight2_im[wn],weight2_im[wn]/weight2_max);

	
	/* get the range of phi */

	double weight1Total1=0.0;
	double weight1Total2=0.0;
	double weight2Total1=0.0;
	double weight2Total2=0.0;

	for(int wn=1;wn<=wnMax;wn++){

		weight1Total1 += weight1[wn]/weight1_max;
		weight1Total2 += weight1[wn]/weight1_max * pow(-1,wn);
		weight2Total1 += weight2[wn]/weight2_max;
		weight2Total2 += weight2[wn]/weight2_max * cos(0.50*3.1415926*wn);

	}
	printf("# One Slab phi range(kJ/mol 300K)   %.3f : %.3f (%.3f) \n",1.0/weight1Total1*300*0.008314,1.0/weight1Total2*300*0.008314,-2.0*weight1_max/30.0*300*0.008314);
	printf("# Two Slabs phi range(kJ/mol 300K)  %.3f : %.3f (%.3f) \n",1.0/weight2Total1*300*0.008314,1.0/weight2Total2*300*0.008314,-2.0*weight2_max/30.0*300*0.008314);


	/* rho_z FFT rho_z for 1 or 2 slabs */	
	

	sprintf(f,"rho_FFT.dat");

	if(argc == 2) {
		wnMax = atoi(argv[1]);
		sprintf(f,"rho_FFT_%d.dat",wnMax);
	}
	data = fopen(f,"w");
	fprintf(data,"# %d wave number\n",wnMax);
	fprintf(data,"# z\trho1\trho1_fft\trho2\trho2_fft\tintRho1\tintRno2\n");
	double frho1,frho2;
	double sum_rho1,sum_rho2;

	double frho1_im,frho2_im;	// density corresponding to sin func (im weight)

	sum_rho1 = 0.0;
	sum_rho2 = 0.0;
	for(int i=0;i<nBin;i++){
		frho1 = natoms/zBox/xBox/yBox;
		frho1_im = 0.0;
		frho1 = natoms/zBox/xBox/yBox;
		frho2_im = 0.0;
	
		for(int wn=1;wn<=wnMax;wn++){
			frho1 += weight1[wn]*cos(2.0*PI*wn/zBox*z[i]) *2.0/zBox/xBox/yBox;
			frho1_im += weight1_im[wn]*sin(2.0*PI*wn/zBox*z[i])*2.0/zBox/xBox/yBox;
		}

		for(int wn=1;wn<=wnMax;wn++){
			frho2 += weight2[wn]*cos(2.0*PI*wn/zBox*z[i]) * 2.0/zBox/xBox/yBox;
			frho2_im += weight2_im[wn]*sin(2.0*PI*wn/zBox*z[i])*2.0/zBox/xBox/yBox;
		}
		if(i>0){
			 sum_rho1 += (rho1[i-1]+rho1[i])/2.0*(z[i]-z[i-1]);
			 sum_rho2 += (rho2[i-1]+rho2[i])/2.0*(z[i]-z[i-1]);
		}

		fprintf(data,"%.4f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",z[i],rho1[i],frho1,rho2[i],frho2,sum_rho1,sum_rho2,frho1+frho1_im,frho2+frho2_im);
	}


	fclose(data);

	return 0;
}
