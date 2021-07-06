/* levy walk search inside an annulus */ 
#include <stdio.h>
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr


// #define TARGETS 1000                   //ensemble average over TARGETS
#define PI 3.14159265358979323846
#define RV 1                         //radius of vision (smaller annulus radius)
// rho, delta, sigma parameters

#define rho 0.0001
#define sigma 0.0001
#define Delta 2


//#define L 50                      //linear size (outer radius)
#define L  (sqrt(1/PI*rho))
//#define LC 5                      //distance from closest target (from center)

#define LC (RV*(Delta+1))

#define scale (sigma*RV)

#define TOTALRANGE 999999999999999999
#define TOTALDISTANCE 100000000         // total distance traveled before stopping
//#define TOTALDISTANCE 2000000         // total distance traveled before stopping
#define LARGESTFLIGHT (L*1000)       // maximum levy step size 

#define R0 1                         //smallest step


// the search space is a 2D torus of size LxL
// the periodic boundary conditions (BC) are chosen so that
// the coordinates x and y can vary from -L/2 to +L/2.
// The target is always located at the origin (0,0)
// This choice of target location makes the code run faster


#define MU_INC 0.1                   //step for incrementing mu
#define MAX_MU_ENTRIES 100           //this should be > 3/MU_INC


// library functions 
double drand48();                   
void exit(int status);


long n;
static double travel;                // travel distance 
static double x,y;                   //searcher coordinate 
static double mu;                    //levy index alpha + 1 
static int tt;                       // dummy for easier type casting

static double
distance_histogram[MAX_MU_ENTRIES];  // table of results
static long
target_histogram[MAX_MU_ENTRIES];    // table of results
static long
flight_histogram[MAX_MU_ENTRIES];    // debugging
static long
inside_histogram[MAX_MU_ENTRIES];    // debugging
static long
outside_histogram[MAX_MU_ENTRIES];    // debugging
 


/*modified function
double rng_levy48(double alpha, double rr){
  double ee, phi;
  double mu=alpha-1;
  double mu1=mu-1;
  double xmu=1/mu;
  double xmu1=xmu-1;
  phi=(drand48()-0.5)*PI;
  ee=-log(drand48());
  //double f;
  return rr*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
  }
*/

/*double valor(void){
	printf("%lf", drand48());
}*/

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

if(f>TOTALRANGE || f<0){
//printf("Deu fora do alcance, f = %lf\n", f);
do{
phi=(drand48()-0.5)*PI;
ee=-log(drand48());
f = rr*sin(mu*phi)/pow(cos(phi),xmu)*pow(cos(phi*mu1)/ee,xmu1);
}while(f>TOTALRANGE || f<0);
//printf("Novo valor de f = %lf\n", f);
}
//printf("valor f=%lf\n", f);
return f;
}


void  initialize_search(){
  double phi;                        // random angle
  //Boundary condition 1
  //x=LC;
  
  //Boundary condition 2
  x=LC;
  y=0;

  // walker always starts at this position.
  
  travel=0;   // initialize travel distance 
  
  //printf("%lf %lf \n",x,y);
}


void find_target(){
  
  //unsigned long here,this;
  //double rrx,rry;                        // random number generation
  double rry;
  double vx,vy;                      // velocity unit vector components 
  double ell;                        // levy walk jump size
  double phi;                        // velocity  angle
  short targetnotfound;              // boolean for while loop
  double xnew,ynew;                  // farthest away possible new searcher position
  double cx, cy,t,t2;                // parametrization variables, see below
  double delta,discriminant;         // for the quadratic formula
  double a,b,c;                      // quadratic formula... 


  targetnotfound=1;


  //  drand48(); // this call is due to a peculiarity of drand48()...

  
  // the following loop only stops upon finding a target
  ell=0;
  while (targetnotfound){

    rry= LARGESTFLIGHT+1;
    while (rry>LARGESTFLIGHT)
      {
    //rrx=drand48();

    //rry=R0*exp(log(rrx)*(1/(1-mu)));
	rry=rng_levy48(mu, scale);
      }
        //printf("mu=%lf flight=%lf  \n",mu,rry);
    
    ell=rry;
    flight_histogram[tt=mu/MU_INC]++;    

    
    phi=drand48()*PI;
    vx=cos(phi);vy=sin(phi);
    
    xnew=x+ell*vx; // 1.000174347
    ynew=y+ell*vy; // 0.0000668769
    
    // parametrize in t according to 
    // x(t) = x + t (xnew-x) = x + cx t
    // y(t) = y + t (ynew-x) = t + cy t
    // this means that t=0 gives the original position and
    // t=1 gives the final position. 

    cx=xnew-x; // 0.000074347
    cy=ynew-y; // 0.0000668769

    // Now let us write (x(t))^2 + (y(t))^2= radius^2 and
    // find the values of t.
    // if t is complex then there is no intersection
    // it t is real and inside  [0,1] then there is intersection
    // with the circle of the given radius of the line segment

    // So x^2+ cx^2 t^2 + 2 cx x t + y^2 + cy^2 t^2 + 2 cy y t -R^2=0
    // t^2(cx^2 + cy^2) + t(2 cx x + 2 cy y) + x^2 + y^2 - R^2 =0
    // So in the quadratic formula we will have
    a=  SQR(cx)  + SQR(cy); // 0.00000000999999616261
    b= 2 *cx*x + 2*cy*y; // 0.0001487088694
    // c= x^2 + y^2 - R^2
    //printf(" (%lf,%lf) --> (%lf,%lf), increment = %lf\n" ,x,y,xnew,ynew, ell);
        
    // start with outer radius L

    if (SQR(xnew)+SQR(ynew)>=SQR(L))
      { //target found
	c= SQR(x) + SQR(y) - SQR(L);
	discriminant = SQR(b) - 4*a*c;// 4*a*c=0.00000003999542465218984984
	//printf("discriminant = %lf\n a = %lf\n cx=%lf\n, cy=%lf\n", discriminant, a, cx, cy);
	if (discriminant<0) {
		printf("y=%lf, x=%lf, ynew=%lf, xnew=%lf, phi=%lf, vx=%lf, vy=%lf\n", y, x, ynew, xnew, phi, vx, vy);
		printf("a=%lf, b=%lf, c=%lf\n", a, b, c);
		printf("rho=%lf, sigma=%lf, delta=%lf, L=%lf, ell=%lf, LC=%lf, rry=%lf", rho, sigma, Delta, L, scale, LC, rry);
		printf("\n Serious discriminant error for outer radius %lf \n", discriminant); exit(0);}
	delta=sqrt(discriminant);
	
	t= ( - b + delta )/(2*a);
	//printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf \n" ,x,y,xnew,ynew,t); 
	if ((-0.001<=t)&&(t<=1.001))
	  {
	    targetnotfound=0;
	    travel+=ell*t; // t is the fraction traversed
	    outside_histogram[ tt=mu/MU_INC]++;
	  }
	else
	  {
	    printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf \n" ,x,y,xnew,ynew,t); 
	    t= ( - b - delta )/(2*a);
	    printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf \n" ,x,y,xnew,ynew,t); 

	    printf("%lf", ell);
	    if ((-0.001<=t)&&(t<=1.001))
	      {
		targetnotfound=0;
		travel+=ell*t; // t is the fraction traversed
		outside_histogram[ tt=mu/MU_INC]++;
	      }
	    else {printf("\n Serious error regarding outer radius\n"); exit(0);}
	  }
      }
    else
      // so walker is inside larger radius
      // in this case either there is zero or two intersections
      // with a very small probability  of a tanget intersection
      {
	c= SQR(x) + SQR(y) - SQR(RV);
	discriminant = SQR(b) - 4*a*c;
	if (discriminant>=0)
	  // this means there is an intersection
	  {
	    delta=sqrt(discriminant);
	    t= ( - b + delta )/(2*a);
	    t2=( - b - delta )/(2*a);
	    if ((0<t2)&&(t2<t)&&(t<1)) t=t2;
	    // the above is to choose the entry rather than exit intersection

	    if ((0<=t) && (t<=1))
	      // target found 
	      {
//				printf("inside!!  (%lf,%lf) --> (%lf,%lf) target at t= %lf \n" ,x,y,xnew,ynew,t);
		targetnotfound=0;
		travel+=ell*t; // t is the fraction traversed
		inside_histogram[ tt=mu/MU_INC]++;
	      }
	    // else the the intersection occurs outside t\in[0,1]
	    // so that it is not physical
	  }
	// else the negative discriminant means there no intersection
      }
    {
      //      printf(" (%lf,%lf) --> (%lf,%lf) target at t= %lf \n" ,x,y,xnew,ynew,t);
      //      printf(" At tt= %lf the walker is at (%lf,%lf) \n" ,t,x+cx*t,y+cy*t);
      //      exit(0);
	}
    x=xnew;
    y=ynew;
    if (targetnotfound) travel+=ell;
  }
}
  
void main(){
//int i=0;

  //initialize result array
  for (mu=1.1;mu<=3.1;mu+=MU_INC) {
    distance_histogram[ tt=mu/MU_INC]=0;
    target_histogram[ tt=mu/MU_INC]=0;
    flight_histogram[ tt=mu/MU_INC]=0;
    inside_histogram[ tt=mu/MU_INC]=0;
    outside_histogram[ tt=mu/MU_INC]=0;
  }
  
  

  
  
  for (mu=1.1;mu<=3.1;mu+=MU_INC){

	  while (distance_histogram[tt=mu/MU_INC]<TOTALDISTANCE){

	initialize_search(RV); //put searcher in the right position
	find_target(0.1);  // search until target found
	distance_histogram[tt=mu/MU_INC]+=travel; // sum the distances and store
	target_histogram[tt=mu/MU_INC]++;
      }
  }
  fflush(stdout);
  

  /*printf(" mu  eta         distance     targets      numer-of-flights inside outside\n");
  
  print result array
  for (mu=1.1;mu<=3.1;mu+=MU_INC){
    printf("%lf %lg %lf %ld %ld %ld %ld\n",
	   mu,
	   target_histogram[ tt=mu/MU_INC]/distance_histogram[ tt=mu/MU_INC]*(SQR(L)/RV) ,
	   distance_histogram[ tt=mu/MU_INC],
	   target_histogram[ tt=mu/MU_INC],
	   flight_histogram[ tt=mu/MU_INC],
	   inside_histogram[ tt=mu/MU_INC],
	   outside_histogram[ tt=mu/MU_INC]);

    fflush(stdout);
  }*/

//ell = Scale parameter
//RE = External Radius
//BC1 = first boundary condition defined as x=LC
//BC2 = second boundary condition defined as x=R0*1.0001

   FILE *arq;
    arq = fopen("concentric_levy.csv", "w+");
    fprintf(arq, "mu,eta,distance,targets,number-of-flights,inside,outside\n");
    for (mu=1.1;mu<=3.1;mu+=MU_INC){
    	fprintf(arq, "%lf,%lg,%lf,%ld,%ld,%ld,%ld\n",
	   mu,
	   target_histogram[ tt=mu/MU_INC]/distance_histogram[ tt=mu/MU_INC]*(SQR(L)/RV) ,
	   distance_histogram[ tt=mu/MU_INC],
	   target_histogram[ tt=mu/MU_INC],
	   flight_histogram[ tt=mu/MU_INC],
	   inside_histogram[ tt=mu/MU_INC],
	   outside_histogram[ tt=mu/MU_INC]);
    }
    fclose(arq);

}
  








