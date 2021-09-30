#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846
#define ALPHA_INCREMENT 0.1
#define ELL_INCREMENT 0.1
#define PRECISION_INTERVAL 0.000000000000001

static double sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double drand48();


//////////////////////////
/*Função de Levy dependente de beta*/

double rng_levy48(double alpha, double rr, double beta){
double ee, phi;
double mu=alpha;
double xmu=1/mu;
double xmu1=xmu-1;
double ksi;
phi=(drand48()-0.5)*PI;
ee=-log(drand48());
double zeta=-beta*tan(PI*mu/2.0);
if((1.0-PRECISION_INTERVAL)<mu && mu<(1.0+PRECISION_INTERVAL)){
	ksi=PI/2;
	return (1.0/ksi)*((PI/2.0+beta*phi)*tan(phi)-beta*log((PI/2.0)*ee*cos(phi)/(PI/2.0+beta*phi)));
}else{
	ksi=(1/mu)*atan(-zeta);
	return pow((1+SQR(zeta)),0.5*xmu)*sin(mu*(phi+ksi))/pow(cos(phi),xmu)*pow(cos(phi-mu*(phi+ksi))/ee, xmu1);
}
}



int main(){

	// o alpha assume o nome de index no nosso programa
	double alpha;
	double ell;
	int values;
	double rry;

/* Gerando todos os valores de ell e mu*/
	FILE *arq1;
	arq1 = fopen("dados_levy_teste.csv", "w+");

	for(alpha=0.1;alpha<2.1;alpha+=ALPHA_INCREMENT){
		for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
			if((2.0-PRECISION_INTERVAL<alpha && alpha<2.0+PRECISION_INTERVAL) && (10.0-PRECISION_INTERVAL<ell && ell<10.0+PRECISION_INTERVAL)){
				fprintf(arq1, "(alpha=%.2lf x R0=%.2lf)\n", alpha, ell);
			}
			else{
			fprintf(arq1, "(alpha=%.2lf x R0=%.2lf),", alpha, ell);
			}
		}

	}

	for(values=1;values<10000;values++){
	for(alpha=0.1; alpha<2.1; alpha+=ALPHA_INCREMENT){
	for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
		if((2.0-PRECISION_INTERVAL<alpha && alpha<2.0+PRECISION_INTERVAL) && (10.0-PRECISION_INTERVAL<ell && ell<10.0+PRECISION_INTERVAL)){
			rry = rng_levy48(alpha, ell, 1.00);

			fprintf(arq1, "%lf \n", rry);
		}
		else{

			rry = rng_levy48(alpha, ell, 1.00);
			fprintf(arq1, "%lf,", rry);
		}
	}
	}
	}

	fclose(arq1);

	return 0;
}
