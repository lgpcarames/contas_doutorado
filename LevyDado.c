#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846
#define INDEX_INCREMENT 0.1
#define INDEX_VALUES 4
#define ELL_VALUES 3
#define ELL_INCREMENT 0.1
#define FUNC_LIMIT 9999.0

static double sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double drand48();

#define LARGESTFLIGHT 999999999999999


/*Função de levy modificada
double rng_levy48(double alpha, double rr){
double ee, phi;
double mu=alpha-1;
double mu1=mu-1;
double xmu=1/mu;
double xmu1=xmu-1;
double f;

do{
phi=(drand48()-0.5)*PI;
ee=-log(drand48());
f = rr*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}while(f<0);

return f;
}
*/


/*  Função rng_levy48 na forma da primeira função de Levy do script dos aneis
double rng_levy48(double alpha, double R0){
 //double ee, phi;
 double mu=alpha;
 double f;
 //double mu1=mu-1;
 //double xmu=1/mu;
 //double xmu1=xmu-1;
 //phi=(drand48()-0.5)*PI;
 //ee=-log(drand48());
 f = R0*exp(log(drand48())*(1/(1-mu)));
 //int condicao = f>10000 || f<-10000;
 if(f>99999999 || f<-99999999){
 //printf("Deu fora do alcance, f = %lf\n", f);
 do{
 f = R0*exp(log(drand48())*(1/(1-mu)));
 }while(f>99999999 || f<-99999999);
 //printf("Novo valor de f = %lf\n", f);
 }
 //printf("valor f=%lf\n", f);
 return f;
//  return rr*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}
*/



//////////////////////////
/*Função de Levy dependente de beta*/

double rng_levy48(double alpha, double rr, double beta){
double ee, phi;
double mu=alpha;
double mu1=mu-1;
double xmu=1/mu;
double xmu1=xmu-1;
double delta=0.000000000000001;
double ksi;
phi=(drand48()-0.5)*PI;
ee=-log(drand48());
double zeta=-beta*tan(PI*mu/2.0);
if((1.0-delta)<mu && mu<(1.0+delta)){
	ksi=PI/2;
	return (1.0/ksi)*((PI/2.0+beta*phi)*tan(phi)-beta*log((PI/2.0)*ee*cos(phi)/(PI/2.0+beta*phi)));
}else{
	ksi=(1/mu)*atan(-zeta);
	return pow((1+SQR(zeta)),0.5*xmu)*sin(mu*(phi+ksi))/pow(cos(phi),xmu)*pow(cos(phi-mu*(phi+ksi))/ee, xmu1);
}
}


//////////////////////////

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
}*/


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
	double rry;

/* Gerando todos os valores de ell e mu*/
	FILE *arq1;
	arq1 = fopen("dados_levy.csv", "w+");

	for(index=0.1;index<2.1;index+=INDEX_INCREMENT){
		for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
			if((1.9999999999999<index && index<2.0000000000001) && (9.99999999999999<ell && ell<10.00000000000001)){
				fprintf(arq1, "(alpha=%.2lf x R0=%.2lf)\n", index, ell);
			}
			else{
			fprintf(arq1, "(alpha=%.2lf x R0=%.2lf),", index, ell);
			}
		}

	}

	for(values=10;values<10000;values++){
	for(index=0.1; index<2.1; index+=INDEX_INCREMENT){
	for(ell=0.1;ell<=10;ell/=ELL_INCREMENT){
		if((1.999999999<index && index<2.000000001) && (9.9999999999<ell && ell<10.0000000001)){
//			do{
			rry = rng_levy48(index, ell, 1.00);
//			}while(rry<0 || rry>LARGESTFLIGHT);

			fprintf(arq1, "%lf \n", rry);
		}
		else{

//			do{
			rry = rng_levy48(index, ell, 1.00);
//			}while(rry<0 || rry>LARGESTFLIGHT);
			fprintf(arq1, "%lf,", rry);
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
