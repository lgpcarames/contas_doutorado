#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846
#define INDEX_INCREMENT 0.5
#define INDEX_VALUES 4
#define ELL_VALUES 3
#define ELL_INCREMENT 0.1

double drand48();

double rng_levy48(double alpha, double rr){
  double ee, phi;
  double mu=alpha;
  double mu1=mu-1;
  double xmu=1/mu;
  double xmu1=xmu-1;
  phi=(drand48()-0.5)*PI;
  ee=-log(drand48());
  return rr*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}

int main(){

	// o alpha assume o nome de index no nosso programa
	double index;
	double ell;
	int values;


	FILE *arq1;
	arq1 = fopen("dados_levy.csv", "w+");

	for(index=0.5;index<2.5;index+=INDEX_INCREMENT){
		for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
			if(index == 2.0 && ell == 10){
				fprintf(arq1, "(alpha=%.2lf x ell=%.2lf)\n", index, ell);
			}
			else{
			fprintf(arq1, "(alpha=%.2lf x ell=%.2lf),", index, ell);
			}
		}

	}

	for(values=1;values<20000;values++){
	for(index=0.5; index<2.5; index+=INDEX_INCREMENT){
	for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
		if(index == 2.0 && ell == 10){
			fprintf(arq1, "%lf \n", rng_levy48(index, ell));
		}
		else{

			fprintf(arq1, "%lf,", rng_levy48(index, ell));
		}
	}
	}
	}


	/*
	FILE *arq1;
	arq1 = fopen("arquivotestea.csv", "w+");

	for(values=1;values<10000;values++){
	for(index=0.5; index<2.5; index+=INDEX_INCREMENT){
	for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
		if(index == 2.0 && ell == 10){
			fprintf(arq1, "%lf\n", rng_levy48(index, ell));
		}
		else{

			fprintf(arq1, "%lf,", rng_levy48(index, ell));
		}
	}
	}
	}
*/


	fclose(arq1);

	return 0;
}
