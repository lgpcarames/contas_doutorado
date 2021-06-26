#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846
#define INDEX_INCREMENT 1
#define INDEX_VALUES 4
#define ELL_VALUES 3
#define ELL_INCREMENT 0.1
#define FUNC_LIMIT 9999.0

double drand48();

void valor(void){
	printf("%lf\n", drand48());
}


/*Função de levy modificada*/
double rng_levy48(double alpha, double rr){
double ee, phi;
double mu=alpha-1;
double mu1=mu-1;
double xmu=1/mu;
double xmu1=xmu-1;
double f;
phi=(drand48()-0.5)*PI;
ee=-log(drand48());
f = rr*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
//int condicao = f>10000 || f<-10000;
if(f>10000 || f<-10000){
printf("Deu fora do alcance, f = %lf\n", f);
do{
phi=(drand48()-0.5)*PI;
ee=-log(drand48());
f = rr*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}while(f>10000 || f<-10000);
printf("Novo valor de f = %lf\n", f);
}
//printf("valor f=%lf\n", f);
return f;
}


/*  Função rng_levy48 dependendo de ell
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
*/

/* Função rng_levy48 com ell=1
double rng_levy48(double alpha){
  double ee, phi;
  double mu=alpha+1;
  double mu1=mu-1;
  double xmu=1/mu;
  double xmu1=xmu-1;
  phi=(drand48()-0.5)*PI;
  ee=-log(drand48());
  return sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}
*/
/* Função rng_levy48 com ell=0.1
double rng_levy48(double alpha){
  double ee, phi;
  double mu=alpha;
  double mu1=mu-1;
  double xmu=1/mu;
  double xmu1=xmu-1;
  phi=(drand48()-0.5)*PI;
  ee=-log(drand48());
  return 0.1*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}
*/

/* Função rng_levy48 com ell=10
double rng_levy48(double alpha){
  double ee, phi;
  double mu=alpha;
  double mu1=mu-1;
  double xmu=1/mu;
  double xmu1=xmu-1;
  phi=(drand48()-0.5)*PI;
  ee=-log(drand48());
  return 10*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}*/


int main(){

	// o alpha assume o nome de index no nosso programa
	double index;
	double ell;
	int values;

/* Gerando todos os valores de ell e mu*/
	FILE *arq1;
	arq1 = fopen("dados_levymod3.csv", "w+");

	for(index=1;index<31;index+=INDEX_INCREMENT){
		for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
			if(index == 30 && ell == 10){
				fprintf(arq1, "(alpha=%.2lf x ell=%.2lf)\n", index/10.0, ell);
			}
			else{
			fprintf(arq1, "(alpha=%.2lf x ell=%.2lf),", index/10.0, ell);
			}
		}

	}

	for(values=1;values<10000;values++){
	for(index=1; index<31; index+=INDEX_INCREMENT){
	for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
		if(index == 30 && ell == 10){
			fprintf(arq1, "%lf \n", rng_levy48(index/10.0, ell));
		}
		else{

			fprintf(arq1, "%lf,", rng_levy48(index/10.0, ell));
		}
	}
	}
	}




	/* Gerando valores de ell = 1

	FILE *arq1;
	arq1 = fopen("arquivoteste1.csv", "w+");


	for(index=0.5;index<2.5;index+=INDEX_INCREMENT){

			if(index == 2.0){
				fprintf(arq1, "(alpha=%.2lf)\n", index);
			}
			else{
			fprintf(arq1, "(alpha=%.2lf),", index);
			}


	}


	for(values=1;values<10000;values++){
	for(index=0.5; index<2.5; index+=INDEX_INCREMENT){

		if(index == 2.0){
			fprintf(arq1, "%lf\n", rng_levy48(index));
		}
		else{

			fprintf(arq1, "%lf,", rng_levy48(index));
		}

	}
	}
*/


/* Gerando valores de ell=0.1
	FILE *arq1;
	arq1 = fopen("arquivoteste2.csv", "w+");


	for(index=0.5;index<2.5;index+=INDEX_INCREMENT){

			if(index == 2.0){
				fprintf(arq1, "(alpha=%.2lf)\n", index);
			}
			else{
			fprintf(arq1, "(alpha=%.2lf),", index);
			}


	}


	for(values=1;values<10000;values++){
	for(index=0.5; index<2.5; index+=INDEX_INCREMENT){

		if(index == 2.0){
			fprintf(arq1, "%lf\n", rng_levy48(index));
		}
		else{

			fprintf(arq1, "%lf,", rng_levy48(index));
		}

	}
	}
	*/

	/*
	FILE *arq1;
	arq1 = fopen("arquivoteste3.csv", "w+");


	for(index=0.5;index<2.5;index+=INDEX_INCREMENT){

			if(index == 2.0){
				fprintf(arq1, "(alpha=%.2lf)\n", index);
			}
			else{
			fprintf(arq1, "(alpha=%.2lf),", index);
			}


	}


	for(values=1;values<10000;values++){
	for(index=0.5; index<2.5; index+=INDEX_INCREMENT){

		if(index == 2.0){
			fprintf(arq1, "%lf\n", rng_levy48(index));
		}
		else{

			fprintf(arq1, "%lf,", rng_levy48(index));
		}

	}
	}
*/

	fclose(arq1);

	return 0;
}
