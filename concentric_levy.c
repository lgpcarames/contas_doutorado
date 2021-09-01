#include <stdio.h>
#include <math.h>
#include <time.h>

static float sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


//Defining some parameters
#define PI 3.14159265358979323846

#define RV 1    // Internal annulus radius
#define L 50    // External annulus radius
#define R0 1    // Scale factor
#define LC 1.0001   // Searcher start point

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


// Defining the initialize_search() function

void initialize_search(){
    x=LC;
    y=0;
    travel=0;
}


// Defining the find_target() function
void find_target(){
    double rrx, rry;                // Random number generation
    double vx, vy;                  // Velocity unit vector components
    double ell;                     // Levy walk jump size
    double phi;                     // Velocity angle
    short targetnotfound;           // Boolean for while loop
    double xnew, ynew;              // farthest away possible new searcher position
    double cx, cy,t,t2;             // parametrization variables, see below
	double delta,discriminant;      // for the quadratic formula
	double a,b,c;                   // quadratic formula...

    // setting targetnotfound as false
    targetnotfound=1;


    // the following loop only stops upon finding a target
    ell=0;
    
    while(targetnotfound){

        // filter while loop

        rry = LARGESTFLIGHT+1;
        while(rry<LARGESTFLIGHT || rry<0){
            //rrx=drand48()
            rry=rng_levy48(alpha, R0);
        }
    
        ell=rry;
        flight_histogram[tt=alpha/ALPHA_INC]++;

        phi=drand48()*PI;
        vx=cos(phi);vy=sin(phi);

        xnew=x+ell*vx;
        ynew=y+ell*vy;

        // parametrize in t according to
        // x(t) = x + t (xnew-x) = x + cx t
        // y(t) = y + t (ynew-x) = t + cy t
        // this means that t=0 gives the original position and
        // t=1 gives the final position.

        cx=xnew-x;
        cy=ynew-y;

        // Now let us write (x(t))^2 + (y(t))^2= radius^2 and
        // find the values of t.
        // if t is complex then there is no intersection
        // it t is real and inside  [0,1] then there is intersection
        // with the circle of the given radius of the line segment

        // So x^2+ cx^2 t^2 + 2 cx x t + y^2 + cy^2 t^2 + 2 cy y t -R^2=0
        // t^2(cx^2 + cy^2) + t(2 cx x + 2 cy y) + x^2 + y^2 - R^2 =0
        // So in the quadratic formula we will have
        a=SQR(cx)+SQR(cy);
        b= 2*cx*x+2*cy*y;
        // c= x^2 + y^2 - R^2

        // Searching if target is found in internal radius

        c=SQR(x)+SQR(y)-SQR(RV);

        discriminant=SQR(b)-4*a*c;
        if(discriminant>=0){
            delta=sqrt(discriminant);
            t= ( - b + delta )/(2*a);
            t2=( - b - delta )/(2*a);
            if ((0<t2)&&(t2<t)&&(t<1)) t=t2;

                if ((0<=t) && (t<=1))
                    // target found
                    {
                    targetnotfound=0;
                    travel+=ell*t; // t is the fraction traversed
                    inside_histogram[ tt=alpha/ALPHA_INC]++;
                }
        }
        else{
            //Searching if the target is found in external radius
            c=SQR(x)+SQR(y)-SQR(L);

            discriminant=SQR(b)-4*a*c;

            if(discriminant>=0){
                delta=sqrt(discriminant);
                t= ( - b + delta )/(2*a);
        
                if ((-0.001<=t)&&(t<=1.001)){
                    targetnotfound=0;
                    travel+=ell*t; // t is the fraction traversed
                    outside_histogram[ tt=alpha/ALPHA_INC]++;
                }
            }
            else{
                t= ( - b - delta )/(2*a);
                if ((-0.001<=t)&&(t<=1.001)){
                    targetnotfound=0;
                    travel+=ell*t; // t is the fraction traversed
                    outside_histogram[ tt=alpha/ALPHA_INC]++;
                }
            }
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
