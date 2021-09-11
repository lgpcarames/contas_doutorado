#include <stdio.h>
#include <math.h>


#define X_OUT 1234.0
#define L 50
#define RV 1

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double t1, t2, t3, t4;
double drand48();
#define PI 3.14159265358979323846
#define ALPHA_INC 0.1
#define MAX_ALPHA_ENTRIES 100
#define MAXVALUES 1000
#define R0 1


void interval_solution(double x, double y, double xnew, double ynew, double rv, double l){
    double a, b, c_in, c_out, discriminant_in, discriminant_out, delta_in, delta_out;
    t1=X_OUT;
    t2=X_OUT;
    t3=X_OUT;
    t4=X_OUT;
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



int main(){
	interval_solution(3.32, 0.79, 66.06, 22.49, 1, 50);
	printf("t1=%lf, t2=%lf, t3=%lf, t4=%lf\n", t1, t2, t3, t4);

	return 0;
}
