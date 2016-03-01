#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <stdio.h>
#include <gsl/gsl_rng.h>

double v[3] = {1, 2, 3};
double q[4] = {1, 2, 3, 4};



double ***r_body;


double q_v_prod (double quaternion[4], double vector[3]) {
    double product[4] = {0, 0, 0, 0};
    for (int i = 0; i < 3; i++) {
        product[0] -=  (quaternion[i+1] * vector[i]);
        product[i+1] = quaternion[0] * vector[i+1] ;
    }
    product[1] += quaternion[2] * vector[2] - quaternion[3] * vector[1];
    product[2] += quaternion[3] * vector[0] - quaternion[1] * vector[2];
    product[3] += quaternion[1] * vector[1] - quaternion[2] * vector[0];
    return product[0];
}


int main() {

    std::cout << q_v_prod (q, v) << std::endl;

    return 0;
};
