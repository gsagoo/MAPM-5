
/*****  deleted comment
******  Coded by: G Sagoo
******  deleted comment
******  Seriesamp2.c 4/7/00  header comments above deleted in November 2020
******  Series Solution For A Stream Function
******  LU Factorisation
*****/
/*This program declares Matrix Space Statically (Not dynamically)
Only the PSI is declared as a MAPM -hopefully this will speed things up*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MAPM/m_apm.h"
/*#define TURNINGDATA*/

#define I 400
#define J 15
#define N 3  /*Number of stream function terms*/
#define DT 0.01
#define REYNOLDS 2.0 /*Keep this number less than 10 to stop the programm from running too long*/
#define MAX_DIGITS 32

void forward_step(MAPM (*psi)[I+3][3], double (*A)[5], int n, int j);
void setup_rhs(MAPM (*psi)[I+3][3], double dn, double dt, double R,int n, int j);
void comp_turn_pt(MAPM (*psi)[I+3][3],MAPM (*tp)[5],double time,double dt,int i, int j,int n,int low_hi);
void integrate_psi2(MAPM (*psi)[I+3][3],MAPM (*integrl_data)[2],int count,double dn, int j);
int converge_test(MAPM (*integrl_data)[2]);
void printmatrix(double (*A)[5],FILE *fp);
void printpsi(MAPM (*psi)[I+3][3],int n, int j,double dt,int min, int max, FILE *fp);
MAPM factorial(int);
void LUdecomposition(double (*A)[5]);
void printMAPM(const char *caption,MAPM m);
void outMAPM(MAPM m);

int main(void)
{
	 m_apm_cpp_precision(34);
	 int i, j, n,count,min_max[N];
	 double time, t_track,dn=2.0/(double)I,dt=DT,R=REYNOLDS,c[5],A[I-1][5];
	 MAPM tp[N][5],integrl_data[N][2],psi[N][I+3][3];
	 char file_name[20]="Results.txt";
    FILE *fp;
	 printf("\n#Steps in Eta [-1,+1]: %d",I);
	 printf("\nStepsize in Eta  Dn = \t%lf",dn);
	 printf("\nTime duration of run:  %d",J);
	 printf("\nStepsize in Time Dt = \t%lf",dt);
	 printf("\n#Functions in series solution:\t%d",N);
	 printf("\nReynolds number:\t\t\t%lf",R);
	 printf("\nOutput filename:\t\t%s",file_name);
	
	/*Initialising Coefficient Matrix,
	NOTE: In the EQUATION I multiplied through by R*dn^4*/
    c[0]=c[4]=-0.5;
    c[1]=c[3]=(2.0+(R*dn*dn/dt));
	 c[2]=(-3.0-(2.0*R*dn*dn/dt));
	 for(i=0; i<I-1; ++i)
		for(j=0; j<5; ++j)
			A[i][j]=c[j];
	 A[0][0]= A[0][1]=0.0 ; A[0][2]= (c[0]+c[2]); /*Overwrites the incorrect entries*/
	 A[1][0]= 0.0; A[I-3][4]=0.0;
	 A[I-2][2]=(c[2]+c[4]); A[I-2][3]=A[I-2][4]=0.0;

	for(n=0;n<N;++n){
		min_max[n]=0;
		for(i=0; i<5;++i)
			tp[n][i]=0.0;}

	 /*Initialise psi with BOUNDARY and INITIAL conditions*/
	 for(n=0; n<N; ++n){
		 for(i=0; i<I+3; ++i)
			psi[n][i][0]=0.0;/*Initial condition psi@(t=0) = 0*/
		 for(j=0; j<3;++j)
		 { /*Boundary condition psi@(Eta=-1,+1) =0*/
			psi[n][1][j]=0.0;/*The other boundary is enforced in setup_rhs*/
			psi[n][I+1][j]=0.0;}}
	for(n=0;n<N;++n)
		for(i=0;i<I+3;++i)
			for(j=0;j<3;++j)
				psi[n][i][j]=0.0;
	LUdecomposition(A);	/*Reduces A[I-1][5] to LU form*/
	i=50;/*Position in fluid where turning pt is monitored- (-1.0+i*dn)*/
	printf("\nIntegral of PSI[i]^2 at Eta = %lf ",(-1.0+(i*dn)));
	/*printf("\nPsi\tMinimum\t\tMaximum\t\tTime_min\tTime_max\tPeriod\n");*/

	for(time=0.0,t_track=0.25,j=0,count=0;time<(double)J; j=(j+1)%3)
	{
		time+=DT;
		for(n=0; n<N; ++n)
			{	
				setup_rhs(psi,dn, dt, R, n, j);/*psi[n][2->I][j+1%3] will store RHS of equation*/
				forward_step(psi,A,n,(j+1)%3); /*Final PSI*/

#ifdef TURNINGDATA/*THE && OPERATOR ISN'T DEFINED FOR MAPM*/
switch(min_max[n]){
	case 0:{/*(j-1)%3 is wrong and (j-1+3)%3 is correct*/
		if((psi[n][i+2][(j+2)%3]>psi[n][i+2][j])&&(psi[n][i+2][(j+1)%3]>psi[n][i+2][j]))
			{
				comp_turn_pt(psi,tp,dt,time-dt,j,n,min_max[n]);
				min_max[n]=1;
			}
		break;}
	case 1:{
		if((psi[n][i+2][j]>psi[n][i+2][(j+2)%3])&&(psi[n][i+2][j]>psi[n][i+2][(j+1)%3]))
			{
				comp_turn_pt(psi,tp,dt,time-DT,j,n,min_max[n]);
				min_max[n]=0;/*THE PRINT STATEMENT BELOW IS WRONG SHOULD BE printMAPM*/
				printf("\n%d\t%lf\t%lf\t%lf\t%lf\t%lf",n,tp[n][3],tp[n][4],tp[n][1],tp[n][2],(tp[n][1]-tp[n][0]));
				tp[n][0]=tp[n][1];
				break;
			}}}
#endif
			}
		if(fabs(time-t_track)<0.000001)
			{	
				integrate_psi2(psi,integrl_data,count,dn,(j+1)%3);
				for(i=0;i<N;++i)
					{
						printf("\npsi_%d\ttime=%lf\t",i,time);
						outMAPM(integrl_data[i][count]);
					}
				if(converge_test(integrl_data))
					printf("\nAll terms converged\n");
				printf("Significant figures %d",psi[1][10][j].significant_digits());
					/*break;*/
				count=(count+1)%2;
				t_track+=0.25;
			}
	 }
	if( (fp = fopen(file_name, "w") ) == NULL)
		{
			fprintf(stderr,"\nError opening file...exiting\n %s.", file_name);
			exit(1);
		}
	 /*printmatrix(A,fp);*/
	 printpsi(psi,N,j,dt,0, I+3,fp);
	 fclose(fp);
	 puts("\n\tProgram terminated sucessfully ");
	 puts("\t******************************");
	 printf("\a\a\a");
	 return(0);
}

void forward_step(MAPM (*psi)[I+3][3], double (*A)[5], int n, int j)
{
	 /*This function evaluate the psi[N][I-1][3] at the next time step.
	 This is done by back substitution using A[I-1][5] (in LU form).
	 This is an implace method, the RHS is overwritten*/
	 int i,u,p,q;
	 MAPM sum,keep;
	 p=q=2;
  /*The backsubstitution with L*/
	 for(i=0;i<I-1;++i)
		{
			sum=0.0;
			for(u=0;(u<p)&&(u<i);++u)
			{
				keep=A[i][(p-1)-u]*psi[n][i-u+1][j];
				sum+=keep.round(MAX_DIGITS);
			}
			psi[n][i+2][j]+=-1.0*sum;
		}

	/*The backsubstitution with U*/
	for(i=I-2;i>-1;--i)
		{
			sum=0.0;
			for(u=0;(u<q)&&(u<I-2-i);++u)
			{	
				keep=A[i][p+1+u]*psi[n][i+3+u][j];
				sum+=keep.round(MAX_DIGITS);
			}
			psi[n][i+2][j]+= -1.0*sum;
			keep=psi[n][i+2][j]/A[i][p];
			psi[n][i+2][j]=keep.round(MAX_DIGITS);
		}
	 /*This then finds PSI at the next time step*/
}

void setup_rhs(MAPM (*psi)[I+3][3], double dn, double dt,double R,int n, int j)
{   /*1.This function sets up the RHS as it would be before operations act on it*/
	 /*2.RHS at time step j is known*/
	 /*3. (-1)^(n+1) = (1.0-2.0*((n+1)%2))*/
	 int i, m;
	 MAPM keep,temp, a,b,c,d,e;
	 /*This is the derivative boundary condition on psi@ (x=-1,+1)*/
	 keep=psi[n][2][j]-2.0*dn*(1.0-2.0*((n+1)%2))/factorial(2*n+1);
	 psi[n][0][j]= keep.round(MAX_DIGITS);
	 keep=psi[n][I][j]+2.0*dn*(1.0-2.0*((n+1)%2))/factorial(2*n+1);
	 psi[n][I+2][j]=keep.round(MAX_DIGITS);
	 for(i=0; i<I-1; ++i)
	{
		 temp=0.0;
		 for(m=0; m<(n+1);++m)/*Sigma part of the RHS*/
			{	
				a=psi[n-m][i+4][j]-psi[n-m][i][j];
				b=psi[n-m][i+3][j]-psi[n-m][i+1][j];
				c=psi[m][i+3][j]+psi[m][i+1][j];
				keep=psi[m][i+2][j]*a;
				d=keep.round(MAX_DIGITS);
				keep=d-b*c;
				e=keep.round(MAX_DIGITS);
				keep=(2.0*m+1.0)*e*R*dn*0.5;
				temp+=keep.round(MAX_DIGITS);
			}
		 keep=(psi[n][i+3][j]-2.0*psi[n][i+2][j]+psi[n][i+1][j])*(R*dn*dn/dt);
		 a=keep.round(MAX_DIGITS);
		 psi[n][i+2][(j+1)%3]=temp+a
		+(0.5*psi[n][i+4][j]-2.0*psi[n][i+3][j]+3.0*psi[n][i+2][j]
		-2.0*psi[n][i+1][j]+0.5*psi[n][i][j]);
	}
	 /*However there is an additional contribution to the RHS via boundary data*/
	 /*See notes*/
	 keep=(-dn*(1.0-2.0*((n+1)%2))/factorial(2*n+1));
	 psi[n][2][(j+1)%3]+=keep.round(MAX_DIGITS);
	 keep=(dn*(1.0-2.0*((n+1)%2))/factorial(2*n+1));
	 psi[n][I][(j+1)%3]+=keep.round(MAX_DIGITS);
}

void comp_turn_pt(MAPM (*psi)[I+3][3],MAPM (*tp)[5],double time,double dt,int i, int j,int n,int low_hi)
{
	/*This function computes the turning point (in the time variable) of PSI
	using a polynomial fit (y(t)= a*t^2 + b*t + c on [0, 2*dt]) on the three values
	y(0) = psi[n][i][(j-1)%3],y(dt)=psi[n][i][j%3],y(2*dt)=psi[n][i][(j+1)%3].
	Data relating to a maximum is stored in tp[n][2,4] and tp[n][0] is used
	externally to calculate PERIOD*/
	MAPM a,b,c;
	c=psi[n][i][(j+2)%3];/*The notion of (j-1)%3 is WRONG! (j-1+3)%3 is right!*/
	b=(psi[n][i][(j+1)%3]-4.0*psi[n][i][j]+3.0*c)/(-2.0*dt);
	a=(psi[n][i][(j+1)%3]-c-2.0*dt*b)/(4.0*dt*dt);
	tp[n][1+low_hi]=time-(b*0.5/a);/*Time for the turning point to occur*/
	tp[n][3+low_hi]=c-(b*b*0.25/a);/*Value of PSI at turning point*/
}

void integrate_psi2(MAPM (*psi)[I+3][3],MAPM (*integrl_data)[2],int count,double dn, int j)
{
	/*This function uses Simpson's Rule to compute the Integral of PSI^2
	NOTE: Area = (h/3)*(Y(0)+4*Y(h)+2*Y(2*h)+4*Y(3*h)+.2*Y()+4*Y()+Y() )*/
	int k,i;
	MAPM sigma,keep;
	for(k=0;k<N;++k)
	{
		sigma=0.0;
		for(i=2;i<I-1;i+=2)
		{	
			keep=4.0*(psi[k][i][j]*psi[k][i][j]);
			sigma+=keep.round(MAX_DIGITS);
			keep=2.0*(psi[k][i+1][j]*psi[k][i+1][j]);
			sigma+=keep.round(MAX_DIGITS);
		}
		keep=dn*(sigma+(psi[k][1][j]*psi[k][1][j])+
					4.0*(psi[k][I][j]*psi[k][I][j])+(psi[k][I+1][j]*psi[k][I+1][j]))/3.0;
		integrl_data[k][count]=keep.round(MAX_DIGITS);
	}
}

int converge_test(MAPM (*integrl_data)[2])
{
	int n;
	MAPM a="1.0e-18",b,c=0.99999,d=1.00001;
	for(n=0;n<N;++n)
		{
		if(integrl_data[n][1]<a)
			continue;
		b=integrl_data[n][0]/integrl_data[n][1];
		if((b<c)||(b>d))
			break;
		}
	if(n==N)
		return(1);
	else
		{
		printf("\n#Psi's converged: ");
		if(n==0)
			printf("NONE!\n");
		for(n-=1;n>-1;--n)
			printf("%d ",n);
		return(0);
		}
}

void printmatrix(double (*A)[5],FILE *fp)
{
	int i,j;
	for(i=0; i<I-1; ++i)
		{
	    	for(j=0; j<5; ++j)
				fprintf(fp,"%lf\t",A[i][j]);
	    	fprintf(fp,"\n");
		}
	fprintf(fp,"\n\n");
}

void printpsi(MAPM (*psi)[I+3][3],int n, int j,double dt,int min, int max, FILE *fp)
{
	int i,k;
	char mBuf[1000];
	/*fprintf(fp,"\n%lf\t",dt*j);*/
	fprintf(fp,"\n");
	for(i=min; i<max; ++i)
	{
		for(k=0;k<n;++k)
		{
			psi[k][i][j].toString(mBuf,20);
			fprintf(fp,"%s\t",mBuf);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n\n");

}

MAPM factorial(int n)
{
	 unsigned long i;
	 MAPM ans;
	 ans=n;
	 if(n<2)
		{	
			ans=1.0;
			return(ans);
		}
	 for(i=(n-1); i>0; --i)
		ans*=(ans-1.0);
	 return(ans);
}

void LUdecomposition(double (*A)[5])
{
	 int i,j,k,p,q;
	 p=q=2;
	 for(i=0;i<(I-2);++i)
		{
			for(k=i+1;(k<i+p+1)&&(k<I-1);++k)
				{	
					A[k][p+(i-k)]/=A[i][p];
					for(j=i+1;(j<i+q+1)&&(j<I-1);++j)
						A[k][p+(j-k)]+= -1.0*A[k][p+(i-k)]*A[i][p+(j-i)];
				}
		}
}

void printMAPM(const char *caption,MAPM m)
{
	char mBuf[100];
	m.toString(mBuf,20);
	printf("\n%s%s",caption,mBuf);
}

void outMAPM(MAPM m)
{
	char mBuf[100];
	m.toString(mBuf,20);
	printf("Integral=%s",mBuf);
}