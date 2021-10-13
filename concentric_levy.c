#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>



/*
 * Notas
 * Gandhi Viswanathan11:48
2 novas colunas
target interior
target exterior
Gandhi Viswanathan11:55
4 colunas
interior/total
exterior/total
Gandhi Viswanathan11:57
alfa=0.1 (mu=1.1) tem que ter 50% no anel interno*/
static double sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


//Defining some parameters
#define PI 3.14159265358979323846

#define RV 1    // Internal annulus radius
#define L 50    // External annulus radius
#define R0 1    // Scale factor
#define LC 1.0001  // Searcher start point
#define PRECISION_INTERVAL 0.0000000000000001

#define X_OUT 12345.0
#define TOTALDISTANCE 10000000  // total distance before stoping
#define LARGESTFLIGHT (L*1000)  // maximum levy step size


#define ALPHA_INC 0.1  // step for incrementing alpha
#define MAX_ALPHA_ENTRIES 100  //this should be >3/ALPHA_INC

// defining the function drand48()
double drand48();
void exit(int status);

// defining some variables
static double travel; // travel distance
static double x,y; // searcher position
static double alpha; // levy index alpha

static int tt; // dummy for easy type casting

static double t1=X_OUT;
static double t2=X_OUT;
static double t3=X_OUT;
static double t4=X_OUT;

static double distance_histogram[MAX_ALPHA_ENTRIES];
static long target_histogram[MAX_ALPHA_ENTRIES];
static long flight_histogram[MAX_ALPHA_ENTRIES];
static long inside_histogram[MAX_ALPHA_ENTRIES];
static long outside_histogram[MAX_ALPHA_ENTRIES];

/* Defining the levy alpha function
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
	ksi=PI/2.0;
	return rr*(1.0/ksi)*((PI/2.0+beta*phi)*tan(phi)-beta*log((PI/2.0)*ee*cos(phi)/(PI/2.0+beta*phi)))+(2.0/PI)*beta*rr*log(rr);
}else{
	ksi=(1/mu)*atan(-zeta);
	return rr*pow((1+SQR(zeta)),0.5*xmu)*sin(mu*(phi+ksi))/pow(cos(phi),xmu)*pow(cos(phi-mu*(phi+ksi))/ee, xmu1);
}
}

    
double pickmin(double n1, double n2, double n3, double n4){
    double nn1, nn2, nn3, nn4;
    nn1=n1;nn2=n2;nn3=n3;nn4=n4;

    if(nn1<0.0) nn1=X_OUT; if(nn2<0.0) nn2=X_OUT; if(nn3<0.0) nn3=X_OUT; if(nn4<0.0) nn4=X_OUT;
    double minimum_n = fmin(fmin(nn1, nn2), fmin(nn3, nn4));

    return minimum_n;
}


void interval_solution(double x, double y, double xnew, double ynew, double rv, double l){
    //long double a, b, c_in, c_out, discriminant_in, discriminant_out, delta_in, delta_out;
	double a, b, c_in, c_out, discriminant_in, discriminant_out, delta_in, delta_out;
    //long double cx, cy;
    double cx, cy;

    cx=xnew-x;
    cy=ynew-y;
    a=SQR(cx)+SQR(cy);

    b=2*cx*x+2*cy*y;

    c_in = SQR(x)+SQR(y)-SQR(RV);
    c_out = SQR(x)+SQR(y)-SQR(L);

    discriminant_in=SQR(b)-4*a*c_in;
    discriminant_out=SQR(b)-4*a*c_out;

    if(discriminant_in>0){
        delta_in = sqrt(discriminant_in);
        t1=(-b+delta_in)/(2*a);
        t2=(-b-delta_in)/(2*a);
    }

    if(discriminant_out>0){
        delta_out = sqrt(discriminant_out);
        t3=(-b+delta_out)/(2*a);
        t4=(-b-delta_out)/(2*a);
    }
}

void annulus_intersect(double tt1, double tt2, double tt3, double tt4, double tt_min){
	//double el=0.0000000000000001;
	short bool_1 = (tt1-PRECISION_INTERVAL)<tt_min && tt_min<(tt1+PRECISION_INTERVAL);
	short bool_2 = (tt2-PRECISION_INTERVAL)<tt_min && tt_min<(tt2+PRECISION_INTERVAL);
	short bool_3 = (tt3-PRECISION_INTERVAL)<tt_min && tt_min<(tt3+PRECISION_INTERVAL);
	short bool_4 = (tt4-PRECISION_INTERVAL)<tt_min && tt_min<(tt4+PRECISION_INTERVAL);
	if(bool_1 || bool_2){
		inside_histogram[tt=alpha/ALPHA_INC]++;
	}else if(bool_3 || bool_4){
		outside_histogram[tt=alpha/ALPHA_INC]++;
	}else{
		printf("annulus intersect Error.\n");
		exit(0);
	}

}



void initialize_search(){
    x=LC;
    y=0.0;
    travel=0.0;
}


void find_target(){
    double rrx, rry;                // Random number generation
    double vx, vy;                  // Velocity unit vector components
    double ell;                     // Levy walk jump size
    double theta;                     // Velocity angle
    short targetnotfound;           // Boolean for while loop
    double xnew, ynew;            // farthest away possible new searcher position
    double t_min;
    int contador=0;
    double mu=alpha+1;

    // setting targetnotfound as false
    targetnotfound=1;


    // the following loop only stops upon finding a target
    ell=0.0;
    
    while(targetnotfound){
        //rry = LARGESTFLIGHT+1;

		rry= LARGESTFLIGHT+1;
		while (rry>LARGESTFLIGHT)
		{
		rrx=drand48();
		rry=R0*exp(log(rrx)*(1/(1-mu)));
		}

  //      do{
//        while(rry<0){
            //rrx=drand48();
         //rry=rng_levy48(alpha, R0, 1);

         //if(rry<0) rry=-rry;

         //if(rry>LARGESTFLIGHT) rry=LARGESTFLIGHT;
//        }while(rry<0 || rry>LARGESTFLIGHT);

        flight_histogram[tt=alpha/ALPHA_INC]++;
        ell=rry;
        theta = 2.0*drand48()*PI;
        vx=cos(theta); vy=sin(theta);

        xnew=x+ell*vx;
        ynew=y+ell*vy;

    interval_solution(x, y, xnew, ynew, RV, L);

    t_min = pickmin(t1, t2, t3, t4);

    if (0.0<t_min && t_min<1.0 + PRECISION_INTERVAL){
        annulus_intersect(t1, t2, t3, t4, t_min);
        targetnotfound=0;
        travel+=ell*t_min; // t is the fraction traversed

    }


//    if(flight_histogram[tt=alpha/ALPHA_INC]>11409230){
//    	printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf, alpha=%lf \n" ,x,y,xnew,ynew,t_min, alpha);
//    	printf("ell=%lf\n", ell);
//    }
    //printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf, alpha=%lf \n" ,x,y,xnew,ynew,t_min, alpha);
    //printf("ell=%lf\n", ell);
    if(SQR(x)+SQR(y)>SQR(L)){
    	printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf, alpha=%lf \n" ,x,y,xnew,ynew,t_min, alpha);
    	printf("ell=%lf\n", ell);
    	contador++;
    	if(contador>10){
    		printf("Parou!\n");
    		printf("%ld", flight_histogram[tt=alpha/ALPHA_INC]);
    		exit(0);
    	}
    }
    x=xnew;
    y=ynew;


    if(targetnotfound) travel+=ell;


    	//exit(0);
    }
}



void main(){

  // function to evaluate the time elapsed
  clock_t tic = clock();
  time_t curr_time;
  struct tm *info;
  time(&curr_time);
  info = localtime(&curr_time);
  

  //initialize result array
  for (alpha=0.1;alpha<2.1;alpha+=ALPHA_INC) {
    distance_histogram[ tt=alpha/ALPHA_INC]=0;
    target_histogram[ tt=alpha/ALPHA_INC]=0;
    flight_histogram[ tt=alpha/ALPHA_INC]=0;
    inside_histogram[ tt=alpha/ALPHA_INC]=0;
    outside_histogram[ tt=alpha/ALPHA_INC]=0;
  }
  
  

  
  
  for (alpha=0.1;alpha<2.1;alpha+=ALPHA_INC){

    while (distance_histogram[tt=alpha/ALPHA_INC]<TOTALDISTANCE)
    {
	initialize_search(); //put searcher in the right position
	find_target();  // search until target found
	distance_histogram[tt=alpha/ALPHA_INC]+=travel; // sum the distances and store
	target_histogram[tt=alpha/ALPHA_INC]++;
    }
  }
  

  printf("\n alpha, eta, distance, targets, number-of-flights, inside, outside, inside-percent, outside-percent\n");
  
  /*print result array*/
  for (alpha=0.1;alpha<2.1;alpha+=ALPHA_INC){
    printf("%lf %lf %lf %ld %ld %ld %ld, %lf, %lf\n",
    	      alpha,
    	      target_histogram[tt=alpha/ALPHA_INC] / distance_histogram[tt=alpha/ALPHA_INC]*(SQR(L)/RV),
    	      distance_histogram[tt=alpha/ALPHA_INC],
    	      target_histogram[tt=alpha/ALPHA_INC],
    	      flight_histogram[tt=alpha/ALPHA_INC],
    	      inside_histogram[tt=alpha/ALPHA_INC],
    	      outside_histogram[tt=alpha/ALPHA_INC],
    		  (double)inside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC],
    		  (double)outside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC]);
    fflush(stdout);
  }



  /*
  for (alpha=0.1;alpha<2.1;alpha+=ALPHA_INC){
    printf("%lf %lg %lf %ld %ld\n",
	   alpha,
	   target_histogram[ tt=alpha/ALPHA_INC]/distance_histogram[ tt=alpha/ALPHA_INC]*(SQR(L)/RV) ,
	   distance_histogram[tt=alpha/ALPHA_INC],
	   target_histogram[tt=alpha/ALPHA_INC],
	   flight_histogram[tt=alpha/ALPHA_INC]);
    fflush(stdout);
  }*/

  /*Creating the csv file and storage the data in it*/
  FILE * arq;
  arq = fopen("concentric_levy.csv", "w+");
  fprintf(arq, "alpha,eta,distance,targets,number-of-flights,inside,outside, inside-percent, outside-percent\n");
  for (alpha = 1.1; alpha < 3.1; alpha += ALPHA_INC) {
    fprintf(arq, "%lf,%lf,%lf,%ld,%ld,%ld,%ld, %lf, %lf\n",
      alpha,
      target_histogram[tt=alpha/ALPHA_INC] / distance_histogram[tt=alpha/ALPHA_INC]*(SQR(L)/RV),
      distance_histogram[tt=alpha/ALPHA_INC],
      target_histogram[tt=alpha/ALPHA_INC],
      flight_histogram[tt=alpha/ALPHA_INC],
      inside_histogram[tt=alpha/ALPHA_INC],
      outside_histogram[tt=alpha/ALPHA_INC],
	  (double)inside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC],
	  (double)outside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC]);
  }
  fclose(arq);



  /*
  FILE * arq;
  arq = fopen("concentric_levy.csv", "w+");
  fprintf(arq, "mu,eta,distance,targets,number-of-flights\n");
  for (alpha = 0.1; alpha < 2.1; alpha += ALPHA_INC) {
    fprintf(arq, "%lf,%lg,%lf,%ld,%ld\n",
      alpha,
      target_histogram[tt=alpha/ALPHA_INC] / distance_histogram[tt=alpha/ALPHA_INC]*(SQR(L)/RV),
      distance_histogram[tt=alpha/ALPHA_INC],
      target_histogram[tt=alpha/ALPHA_INC],
      flight_histogram[tt=alpha/ALPHA_INC]);
  }
  fclose(arq);*/

  clock_t toc = clock();
  printf("Elapsed: %f seconds, TD=%d\n", (double)(toc - tic)/CLOCKS_PER_SEC, (int)TOTALDISTANCE);


}
