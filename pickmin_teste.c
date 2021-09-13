#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define RAND_MAX 2147483647.0


#define X_OUT 1234.0
#define L 50
#define RV 1

#define PI 3.14159265358979323846

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double t1, t2, t3, t4;
double drand48();
#define PI 3.14159265358979323846
#define ALPHA_INC 0.1
#define MAX_ALPHA_ENTRIES 100
#define MAXVALUES 1000
#define R0 1


double pickmin(double n1, double n2, double n3, double n4){
    double nn1, nn2, nn3, nn4;
    nn1=n1;nn2=n2;nn3=n3;nn4=n4;

    if(nn1<0) nn1=2; if(nn2<0) nn2=2; if(nn3<0) nn3=2; if(nn4<0) nn4=2;
    double minimum_n = fmin(fmin(nn1, nn2), fmin(nn3, nn4));

    return minimum_n;
}


int main(){

    printf("%lf", pickmin(-0.000011, -0.000034, 0.000547, -0.000592));

    return 0;
}