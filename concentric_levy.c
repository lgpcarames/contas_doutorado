#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

static float sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


//Defining some parameters
#define PI 3.14159265358979323846

#define RV 1    // Internal annulus radius
#define L 50    // External annulus radius
#define R0 1    // Scale factor
#define LC 1.0001   // Searcher start point

#define X_OUT 2.00
#define TOTALDISTANCE 1000000000 // total distance before stoping
#define LARGESTFLIGHT (L*1000)  // maximum levy step size


#define ALPHA_INC 0.1  // step for incrementing alpha
#define MAX_ALPHA_ENTRIES 100  //this should be >3/ALPHA_INC

// defining the function drand48()
double drand48();

// defining some variables
static double travel; // travel distance
static double x,y; // searcher position
static double alpha; // levy index alpha
static int tt; // dummy for easy type casting

static double t1;
static double t2;
static double t3;
static double t4;

static double distance_histogram[MAX_ALPHA_ENTRIES];
static long target_histogram[MAX_ALPHA_ENTRIES];
static long flight_histogram[MAX_ALPHA_ENTRIES];
static long inside_histogram[MAX_ALPHA_ENTRIES];
static long outside_histogram[MAX_ALPHA_ENTRIES];

// Defining the levy alpha function
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

    
double pickmin(double n1, double n2, double n3, double n4){
    double minimum_n = fmin(fmin(n1, n2), fmin(n3, n4));

    return minimum_n;
}


void interval_solution(double x, double y, double xnew, double ynew, double cx, double cy, double rv, double l){
    double a, b, c_in, c_out, discriminant_in, discriminant_out, delta_in, delta_out;
    t1=X_OUT;
    t2=X_OUT;
    t3=X_OUT;
    t4=X_OUT;


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
        if(t1<0){t1=X_OUT;}
        if(t2<0){t2=X_OUT;}
    }

    if(discriminant_out>0){
        delta_out = sqrt(discriminant_out);
        t3=(-b+delta_out)/(2*a);
        t4=(-b-delta_out)/(2*a);
        if(t3<0){t3=X_OUT;}
        if(t4<0){t4=X_OUT;}
    }
}

//this function will verify what annulus the searcher have found
void annulus_intersect(double n1, double n2, double n3, double n4){
    double tmin = pickmin(n1, n2, n3, n4);

    if(n1 == tmin || n2 == tmin){
        inside_histogram[tt==alpha/ALPHA_INC];
    }
    if(n3 == tmin || n4 == tmin){
        outside_histogram[tt==alpha/ALPHA_INC];
    }
}

// Defining the initialize_search() function

void initialize_search(){
    x=LC;
    y=0;
    travel=0;
}


void find_target(){
    double rrx, rry;                // Random number generation
    double vx, vy;                  // Velocity unit vector components
    double ell;                     // Levy walk jump size
    double phi;                     // Velocity angle
    short targetnotfound;           // Boolean for while loop
    double xnew, ynew, cx, cy;              // farthest away possible new searcher position
    double t_min;


    // setting targetnotfound as false
    targetnotfound=1;


    // the following loop only stops upon finding a target
    ell=0;
    
    while(targetnotfound){
        rry = LARGESTFLIGHT+1;
        while(rry>LARGESTFLIGHT){
            rrx=drand48();
            rry=rng_levy48(alpha, R0);
        }

        ell=rry;
        flight_histogram[tt=alpha/ALPHA_INC];

        phi = drand48()*PI;
        vx=cos(phi); vy=sin(phi);

        xnew=x+ell*vx;
        ynew=y+ell*vy;

        cx=xnew-x;
        cy=ynew-y;

    interval_solution(x, y, xnew, ynew, cx, cy, RV, L);

    t_min = pickmin(t1, t2, t3, t4);

    if (0<t_min<1){
        annulus_intersect(t1, t2, t3, t4);
        targetnotfound=0;
        travel+=ell*t_min; // t is the fraction traversed

    }


    x=xnew;
    y=ynew;
    if(targetnotfound) travel+=ell;

    }



}



void main(){
  int i;
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
  
  

  
  
  for (alpha=1.1;alpha<=3.1;alpha+=ALPHA_INC){

    while (distance_histogram[tt=alpha/ALPHA_INC]<TOTALDISTANCE)
    {
	initialize_search(); //put searcher in the right position
	find_target();  // search until target found
	distance_histogram[tt=alpha/ALPHA_INC]+=travel; // sum the distances and store
	target_histogram[tt=alpha/ALPHA_INC]++;
    }
  }
  

  printf("\n#% mu, eta, distance, targets, number-of-flights, inside outside\n");
  
  //print result array
  for (alpha=1.1;alpha<=3.1;alpha+=ALPHA_INC){
    printf("%lf %lg %lf %ld %ld %ld %ld\n",
	   alpha,
	   target_histogram[ tt=alpha/ALPHA_INC]/distance_histogram[ tt=alpha/ALPHA_INC]*(SQR(L)/RV) ,
	   distance_histogram[tt=alpha/ALPHA_INC],
	   target_histogram[tt=alpha/ALPHA_INC],
	   flight_histogram[tt=alpha/ALPHA_INC],
	   inside_histogram[tt=alpha/ALPHA_INC],
	   outside_histogram[tt=alpha/ALPHA_INC]);
    fflush(stdout);
  }


  //Creating the csv file and storage the data in it
  FILE * arq;
  arq = fopen("concentric_levy.csv", "w+");
  fprintf(arq, "mu,eta,distance,targets,number-of-flights,inside,outside\n");
  for (alpha = 1.1; alpha <= 3.1; alpha += ALPHA_INC) {
    fprintf(arq, "%lf,%lg,%lf,%ld,%ld,%ld,%ld\n",
      alpha,
      target_histogram[tt=alpha/ALPHA_INC] / distance_histogram[tt=alpha/ALPHA_INC]*(SQR(L)/RV),
      distance_histogram[tt=alpha/ALPHA_INC],
      target_histogram[tt=alpha/ALPHA_INC],
      flight_histogram[tt=alpha/ALPHA_INC],
      inside_histogram[tt=alpha/ALPHA_INC],
      outside_histogram[tt=alpha/ALPHA_INC]);
  }
  fclose(arq);

  clock_t toc = clock();
  printf("Elapsed: %f seconds, TD=%d\n", (double)(toc - tic)/CLOCKS_PER_SEC, (int)TOTALDISTANCE);


}
