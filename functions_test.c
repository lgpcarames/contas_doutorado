#include <stdio.h>
#include <math.h>
#include <stdlib.h>
static float sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define X_OUT 2.00
static double t1;
static double t2;
static double t3;
static double t4;
    
double theless(double n1, double n2, double n3, double n4){
    double minimum_n = fmin(fmin(n1, n2), fmin(n3, n4));

    return minimum_n;
}


void interval_solution(double x, double y, double xnew, double ynew, double cx, double cy, double RV, double L){
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


int main(){

    double L=5.74;
    double RV=1.02;
    double x=3.53;
    double y=-2.9;
    double xnovo=-5.72319;
    double ynovo=0.46105;
    double cx=xnovo-x;
    double cy=ynovo-y;


    interval_solution(x, y, xnovo, ynovo, cx, cy, RV, L);

    printf("t1=%.2f, t2=%.2f, t3=%.2f, t4=%.2f\n", t1, t2, t3, t4);

    printf("%lf", theless(t1, t2, t3, t4));

    return 0;
}