/* Levy Searcher Simulation*/


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>


// Macro to obtain the square numbers
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


//Defining some parameters
#define PI 3.14159265358979323846
#define RV 1    // Internal annulus radius
//#define L 50    // External annulus radius
//#define R0 1    // Scale factor
//#define LC 1.0001  // Searcher start point
#define PRECISION_INTERVAL 0.000000000001 // used to compare float numbers

// input values
  // rho, delta, sigma parameters
//#define sigma 0.0001
//#define rho 0.0001
//#define Del 0.0001

//static double L = 1 / sqrt(PI * rho); // external radius
//static double LC = RV * (Del + 1); // start point
//static double R0 = sigma * RV;


#define X_OUT 12345.0 // default value to numbers outside the interval [0,1]
#define TOTALDISTANCE 100000  // total distance before stoping
#define LARGESTFLIGHT 9999999999  // maximum levy step size


#define ALPHA_INC 0.1  // step for incrementing alpha
#define MAX_ALPHA_ENTRIES 100  //this should be >3/ALPHA_INC


// calling some functions
double drand48();
void exit(int status);


// defining some variables
static double travel; // travel distance
static double x,y; // searcher position
static double alpha; // levy index alpha
static int tt; // dummy for easy type casting

// solutions for the intersection with the inner ring
static double t1=X_OUT;
static double t2=X_OUT;

// solution for hte intersectio with the outer ring
static double t3=X_OUT;
static double t4=X_OUT;


// call for the results array
static double distance_histogram[MAX_ALPHA_ENTRIES];
static long target_histogram[MAX_ALPHA_ENTRIES];
static long flight_histogram[MAX_ALPHA_ENTRIES];
static long inside_histogram[MAX_ALPHA_ENTRIES];
static long outside_histogram[MAX_ALPHA_ENTRIES];


/* Original stable levy function proposed by Prof. Sergey*/
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



/* Stable levy function with beta dependency
* Reference: Chambers et al. Journal of Statistical Association, 1976
* alpha-> is the alpha index in the stable distribution, with interval defined as (0,2]
* rr-> is the scale factor in the stable distribution
* beta-> is the skewness parameter, with interval defined as [-1, 1]

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
*/

/* Select the smallest positive value defined
*  among four double variables */
double pickmin(double n1, double n2, double n3, double n4){
    double nn1, nn2, nn3, nn4;
    nn1=n1;nn2=n2;nn3=n3;nn4=n4;

    if(nn1<0.0) nn1=X_OUT; if(nn2<0.0) nn2=X_OUT; if(nn3<0.0) nn3=X_OUT; if(nn4<0.0) nn4=X_OUT;
    double minimum_n = fmin(fmin(nn1, nn2), fmin(nn3, nn4));

    return minimum_n;
}


/*Given the initial position (x, y), the final position (xnew, ynew),
* the radius of the inner (rv) and outer ring (l)
* this function calculates if there will real solutions and
* store the real solutions in t1, t2, t3 and t4*/
void interval_solution(double x, double y, double xnew, double ynew, double ll){
    // Defining variables that does not depend on the inner or outer radius
	double a, b, c_in, c_out, discriminant_in, discriminant_out, delta_in, delta_out;
    double cx, cy;

    cx=xnew-x;
    cy=ynew-y;

    a=SQR(cx)+SQR(cy);
    b=2*cx*x+2*cy*y;

    // Variables that depends on the inner radius
    c_in = SQR(x)+SQR(y)-SQR(RV);
    discriminant_in=SQR(b)-4*a*c_in;

    // Variables that depends on the outer radius
    c_out = SQR(x)+SQR(y)-SQR(ll);
    discriminant_out=SQR(b)-4*a*c_out;

    // Loop to certify if there is a real solution for the inner ring
    if(discriminant_in>0){
        // If it does, then obtain t1 and t2
        delta_in = sqrt(discriminant_in);
        t1=(-b+delta_in)/(2*a);
        t2=(-b-delta_in)/(2*a);
    }

    // Loop to certify if there is a real solution for the outer ring
    if(discriminant_out>0){
        // If it does, then obtain t3 and t4
        delta_out = sqrt(discriminant_out);
        t3=(-b+delta_out)/(2*a);
        t4=(-b-delta_out)/(2*a);
    }
}

/*If the searcher finds a target, this functions returns what ring the searcher has found*/
void annulus_intersect(double tt1, double tt2, double tt3, double tt4, double tt_min){
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



void initialize_search(double delta){
	double LC = RV * (delta + 1); // start point
    x=LC;
    y=0.0;
    travel=0.0;
}


void find_target(double sigma, double delta, double rho, double total_distance){
    double rrx, rry;                // Random number generation
    double vx, vy;                  // Velocity unit vector components
    double ell;                     // Levy walk jump size
    double theta;                     // Velocity angle
    short targetnotfound;           // Boolean for while loop
    double xnew, ynew;            // farthest away possible new searcher position
    double t_min;
    int count=0;

    double L = 1 / sqrt(PI * rho); // external radius
    double LC = RV * (delta + 1); // start point
    double R0 = sigma * RV;



    // setting targetnotfound as false
    targetnotfound=1;


    // the following loop only stops upon finding a target
    ell=0.0;

    while(targetnotfound){

    // Starts the power-law generator
      // Adaptating the power-law generator to receive alpha index values
    	//double mu=alpha+1;

    	//rry= LARGESTFLIGHT+1;



		//do{
		//rrx=drand48();
		//rry=R0*exp(log(rrx)*(1/(1-mu)));
		rry=rng_levy48(alpha, R0);

		//rry=rng_levy48(alpha, R0, 1);
		if(rry<0) rry = -rry;
		if(rry>LARGESTFLIGHT) rry = LARGESTFLIGHT;
//		}while (rry<0 || rry>LARGESTFLIGHT);



    flight_histogram[tt=alpha/ALPHA_INC]++;

    ell=rry;
    theta = 2.0*drand48()*PI;
    vx=cos(theta); vy=sin(theta);

    xnew=x+ell*vx;
    ynew=y+ell*vy;

    interval_solution(x, y, xnew, ynew, L);

    t_min = pickmin(t1, t2, t3, t4);

    if (0.0<t_min && t_min<1.0 + PRECISION_INTERVAL){
        annulus_intersect(t1, t2, t3, t4, t_min);
        targetnotfound=0;
        travel+=ell*t_min; // t is the fraction traversed

    }

    // Verifying if there is a compare precision bug
    if(SQR(x)+SQR(y)>SQR(L)){
      // The searcher is outside the outer ring without the intersection being correctly identified
    	printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf, alpha=%lf \n" ,x,y,xnew,ynew,t_min, alpha);
    	printf("ell=%lf\n", ell);
    	count++;
    	if(count>10){
    		printf("Stopped!\n");
    		printf("%ld", flight_histogram[tt=alpha/ALPHA_INC]);
    		exit(0);
    	}
    }
    x=xnew;
    y=ynew;

	if((distance_histogram[ tt=alpha/ALPHA_INC]+travel)>total_distance){
		targetnotfound=0;
	}
    if(targetnotfound) travel+=ell;
    }
}



void starting_simulation(float alpha_min, float alpha_max){

  for (alpha=alpha_min;alpha<alpha_max;alpha+=ALPHA_INC) {
	distance_histogram[ tt=alpha/ALPHA_INC]=0;
	target_histogram[ tt=alpha/ALPHA_INC]=0;
	flight_histogram[ tt=alpha/ALPHA_INC]=0;
	inside_histogram[ tt=alpha/ALPHA_INC]=0;
	outside_histogram[ tt=alpha/ALPHA_INC]=0;
  }
}

void simulation(double alpha_min, double alpha_max, double sigma, double delta, double rho, double total_distance){

	double L = 1 / sqrt(PI * rho); // external radius
	double LC = RV * (delta + 1); // start point
	double R0 = sigma * RV;

  for (alpha=alpha_min;alpha<alpha_max;alpha+=ALPHA_INC){

	while (distance_histogram[tt=alpha/ALPHA_INC]<total_distance)
	{
	  initialize_search(delta); //put searcher in the right position
	  find_target(sigma, delta, rho, total_distance);  // search until target found
	  //if(alpha>0.9){
		//printf("Aqui...");
	  //}
	  distance_histogram[tt=alpha/ALPHA_INC]+=travel; // sum the distances and store
	  target_histogram[tt=alpha/ALPHA_INC]++;
	}
  }


}


void result_save(double alpha_min, double alpha_max, double sigma, double delta, double rho, int n ){
	 FILE * arq;
	 char filename[5000];
	 sprintf(filename, "/home/lucas/eclipse-workspace/Estudos_Doutorado/contas_doutorado/PRL_Results/levy_n=%d.csv", n);
	  arq = fopen(filename, "w+");
	  fprintf(arq, "alpha,eta,distance,targets,number-of-flights,inside,outside,inside-percent,outside-percent\n");
	  for (alpha = alpha_min; alpha < alpha_max; alpha += ALPHA_INC) {
	    fprintf(arq, "%lf,%lf,%lf,%ld,%ld,%ld,%ld,%lf,%lf\n",
	      alpha,
	      target_histogram[tt=alpha/ALPHA_INC] / distance_histogram[tt=alpha/ALPHA_INC]*(1/(rho*RV)),
	      distance_histogram[tt=alpha/ALPHA_INC],
	      target_histogram[tt=alpha/ALPHA_INC],
	      flight_histogram[tt=alpha/ALPHA_INC],
	      inside_histogram[tt=alpha/ALPHA_INC],
	      outside_histogram[tt=alpha/ALPHA_INC],
		  (double)inside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC],
		  (double)outside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC]);
	  }



	  fprintf(arq, "\n d=%.1e_s=%.1e_r=%.1e", sigma, delta, rho);

	  fclose(arq);

}



void main(){

  // function to evaluate the time elapsed
  clock_t tic = clock();
  time_t curr_time;
  struct tm *info;
  time(&curr_time);
  info = localtime(&curr_time);
  
	double sigma;
	double rho;
	double Del;


	float alpha_min = 0.1;
	float alpha_max = 1.8;
	float tt_d = 1000.0;


	int i;
	for(sigma=0.01, Del=0.01,rho=0.01,i=1;sigma>=0.0001, Del>=0.0001,rho>=0.0001,i<=27;sigma*=0.1, Del*=0.1, rho*=0.1,i++){
	starting_simulation(alpha_min, alpha_max);

	simulation(alpha_min, alpha_max, sigma, Del, rho, tt_d);

	result_save(alpha_min, alpha_max, sigma, Del, rho, i);
	}
  /*
  printf("\n alpha, eta, distance, targets, number-of-flights, inside, outside, inside-percent, outside-percent\n");
  
  Console output result table
  for (alpha=alpha_min;alpha<alpha_max;alpha+=ALPHA_INC){
    printf("%lf %lf %lf %ld %ld %ld %ld, %lf, %lf\n",
    	      alpha,
    	      target_histogram[tt=alpha/ALPHA_INC] / distance_histogram[tt=alpha/ALPHA_INC]*(1/(rho*RV)),
    	      distance_histogram[tt=alpha/ALPHA_INC],
    	      target_histogram[tt=alpha/ALPHA_INC],
    	      flight_histogram[tt=alpha/ALPHA_INC],
    	      inside_histogram[tt=alpha/ALPHA_INC],
    	      outside_histogram[tt=alpha/ALPHA_INC],
    		  (double)inside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC],
    		  (double)outside_histogram[tt=alpha/ALPHA_INC]/(double)target_histogram[tt=alpha/ALPHA_INC]);
    fflush(stdout);
  }
  */



  /*estimating the total time elapsed*/
  clock_t toc = clock();
  printf("Elapsed: %f seconds, TD=%d\n", (double)(toc - tic)/CLOCKS_PER_SEC, (int)TOTALDISTANCE);


}
