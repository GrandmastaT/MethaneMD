#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <stdio.h>
#include <gsl/gsl_rng.h>



// units: nm, ps, au
// i,j,k,l,m,n not free
const int N = 256 ;             // number of lattice points
double rho = 0.6;
//double rho = 0.395052;          // density u/(nm³) in gas at T = 25°C
double T = 25. + 273;           // system temperature in K
double dt = 0.1;               // time step in ps
double ***r_lab;                // lattice point
double ***v;                    // velocities
double ***a;                    // accelerations
double **w;                       // angular velocities
double **w_t;                   // angular acceleration

double d_CH = 0.1087;           // distance from C to H (in nm)
double angle_HCH = 109.5;       // angle between HCH in grad
double hh = d_CH * sqrt(2.* (1.-cos(angle_HCH)));	//length between two H atoms
double mBox = hh/sqrt(2.);      // length of molecule cell
double lBox = std::pow(N/rho, 1./3.);   // calculate box length
double **q;                     //list of quaternion paramters q[0], q[1], q[2], q[3]
double **q_t;                   // q derivative after time
double **q_tt;                  // q derivative after time²
double **q_d;                   // q time step
double **q_delta;               // quaternion change
double q_norm;                  // normalization factor of quaternion
double ***torque;

double **matrix_A_list;
double **matrix_A;

double mH = 1.008;  // mass of hydrogen in u
double mC = 12.0107;  // mass of hydrogen in u
double mass[5] = {mC, mH, mH, mH, mH};

double Ixx = 8*mH; // moment of inertia
double Iyy = 8*mH;
double Izz = 8*mH;

//const double kb = 1.38064852*std::pow(10, -23);       // boltzmann constant in J/K

double ***r_body;
double sigma_HH = 0.25;     // nm
double sigma_HH_6 = std::pow(sigma_HH, 6);
double sigma_CC = 0.35;     // nm
double sigma_CC_6 = std::pow(sigma_CC, 6);
double sigma_CH = 0.3;      // nm
double sigma_CH_6 = std::pow(sigma_CH, 6);
double epsilon_HH = 15.1;   // K
double epsilon_CC = 33.2;   // K
double epsilon_CH = 22.3902;// K
/*
double lj_CC_12 = 1.12190000;   // epsilon*sigma^12  (in K*Angstrom)   with KONG combination rules
double lj_HH_12 = 9.00030;
double lj_CH_12 = 1.2554300;
double lj_CC_6 = 6.10304;       // epsilon*sigma^6   (in K*Angstrom)   with KONG combination rules
double lj_HH_6 = 3.68652;
double lj_CH_6 = 1.49997;
*/

/*
double q_v_prod (double quaternion[4], double vector[3]) {
    double product[4] = {0, 0, 0, 0};
    for (int i = 0; i < 3; i++) {
        product[0] += - (quaternion[i+1] + vector[i]);
        product[i+1] = quaternion[0] * vector[i+1] ;
    }
    product[1] += quaternion[2] * vector[2] - quaternion[3] * vector[1];
    product[2] += quaternion[3] * vector[0] - quaternion[1] * vector[2];
    product[3] += quaternion[1] * vector[1] - quaternion[2] * vector[0];
    return product[4];
}
*/

void initPositions() { //initialize FFC lattice
    int axisN = 1;

    while (4*axisN*axisN*axisN < N) // look for magic number, sideN is number of lattice points on each box axis
        axisN++;

    if (N/(4*axisN*axisN*axisN) == 1) {
        double stepL = lBox/axisN;              // step length, length between lattice points
        double hstepL = stepL/2.;

        r_lab = new double **[N];   // create lattice points
        r_body = new double **[N];  // create rotation points
        v = new double **[N];
        a = new double **[N];
        w = new double *[N];
        w_t = new double *[N];

        for (int i = 0; i < N; i++) { // create array for lattice points
            r_lab[i] = new double *[5];
            r_body[i] = new double *[5]; 
            v[i] = new double *[5];
            a[i] = new double *[5];
            w[i] = new double [3];
            w_t[i] = new double [3];
            for (int j = 0; j < 5; j++) {            
                r_lab[i][j] = new double [3];
                r_body[i][j] = new double [3]; 
                v[i][j] = new double [3];
                a[i][j] = new double [3];
            }
        }
        for (int i = 0; i < N; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 3; k++)
                    v[i][j][k] = 0;

        int index = 0;
        for (int i = 0; i < 2*axisN; i++) 
            for (int j = 0; j < 2*axisN; j++) 
                for (int k = 0; k < 2*axisN; k++)
                    if ((i+j+k)%2 == 0) {
                        r_lab[index][0][0] = k*hstepL;      // set C
                        r_lab[index][0][1] = i*hstepL;
                        r_lab[index][0][2] = j*hstepL;

                        r_lab[index][1][0] = k*hstepL - mBox*0.5;      // set H1
                        r_lab[index][1][1] = i*hstepL - mBox*0.5;
                        r_lab[index][1][2] = j*hstepL - mBox*0.5;

                        r_lab[index][2][0] = k*hstepL + mBox*0.5;      // set H2
                        r_lab[index][2][1] = i*hstepL - mBox*0.5;
                        r_lab[index][2][2] = j*hstepL + mBox*0.5;

                        r_lab[index][3][0] = k*hstepL - mBox*0.5;      // set H3
                        r_lab[index][3][1] = i*hstepL + mBox*0.5;
                        r_lab[index][3][2] = j*hstepL + mBox*0.5;

                        r_lab[index][4][0] = k*hstepL + mBox*0.5;      // set H4
                        r_lab[index][4][1] = i*hstepL + mBox*0.5;
                        r_lab[index][4][2] = j*hstepL - mBox*0.5;

                        index++;
                    }    
    }
    else
        std::cout << "oida, keine magic number" << std::endl;
}

void setSites() {               // sets sites to a com r_lab coordinate. need to rotate then.
    for (int i = 0; i < N; i++) {
                        r_lab[i][1][0] = r_lab[i][0][0] - mBox*0.5;      // set H1
                        r_lab[i][1][1] = r_lab[i][0][1] - mBox*0.5;
                        r_lab[i][1][2] = r_lab[i][0][2] - mBox*0.5;

                        r_lab[i][2][0] = r_lab[i][0][0] + mBox*0.5;      // set H2
                        r_lab[i][2][1] = r_lab[i][0][1] - mBox*0.5;
                        r_lab[i][2][2] = r_lab[i][0][2] + mBox*0.5;

                        r_lab[i][3][0] = r_lab[i][0][0] - mBox*0.5;      // set H3
                        r_lab[i][3][1] = r_lab[i][0][1] + mBox*0.5;
                        r_lab[i][3][2] = r_lab[i][0][2] + mBox*0.5;

                        r_lab[i][4][0] = r_lab[i][0][0] + mBox*0.5;      // set H4
                        r_lab[i][4][1] = r_lab[i][0][1] + mBox*0.5;
                        r_lab[i][4][2] = r_lab[i][0][2] - mBox*0.5;
                    }    
}

void initRotations() {  //initialize first molecule rotation at random, defines matrix_A_list
    const gsl_rng_type * Seed1;
    gsl_rng * Seed2;
    gsl_rng_env_setup();

    Seed1 = gsl_rng_default;
    Seed2 = gsl_rng_alloc (Seed1);
    
    q = new double *[N];                    // first random quaternion
    q_t = new double *[N];
    q_tt = new double *[N];
    q_d = new double *[N];
    torque = new double **[N];
    q_delta = new double *[N];

    matrix_A = new double *[3];             // matrices for rotations
    for (int k = 0; k < 3; k++)
        matrix_A[k] = new double[3];

    matrix_A_list = new double *[N];        // list of rotation matrices, [molecule][A_ij]
    for (int k = 0; k < N ;k++)
        matrix_A_list[k] = new double [9];

    for (int k = 0; k < N; k++) {           // for every molecule   
        q_norm = 0;
        q[k] = new double [4];
        q_t[k] = new double [4];
        q_tt[k] = new double [4];
        q_d[k] = new double [4];
        q_delta[k] = new double [4];
        torque[k] = new double *[5];
        for (int i = 0; i < 4; i++)             // generate random quaternion
            q[k][i] = (0.5 - gsl_rng_uniform(Seed2)); 

        for (int i = 0; i < 5; i++)
            torque[k][i] = new double [3];

        for (int i = 0; i < 4; i++)             // calculate quaternion normization
            q_norm += q[k][i]*q[k][i];

        for (int i = 0; i < 4; i++)       // unit quaternion
            q[k][i] = q[k][i]/std::sqrt(q_norm);

       // calculation of the rotation matrix 
            matrix_A[0][0] = (q[k][0]*q[k][0] + q[k][1]*q[k][1] - q[k][2]*q[k][2] - q[k][3]*q[k][3]);
            matrix_A[0][1] = 2*(q[k][1]*q[k][2] - q[k][0]*q[k][3]);
            matrix_A[0][2] = 2*(q[k][1]*q[k][3] + q[k][0]*q[k][2]);

            matrix_A[1][0] = 2*(q[k][1]*q[k][2] + q[k][0]*q[k][3]);
            matrix_A[1][1] = (q[k][0]*q[k][0] - q[k][1]*q[k][1] + q[k][2]*q[k][2] - q[k][3]*q[k][3]);
            matrix_A[1][2] = 2*(q[k][2]*q[k][3] - q[k][0]*q[k][1]);

            matrix_A[2][0] = 2*(q[k][1]*q[k][3] - q[k][0]*q[k][2]);
            matrix_A[2][1] = 2*(q[k][2]*q[k][3] + q[k][0]*q[k][1]);
            matrix_A[2][2] = (q[k][0]*q[k][0] - q[k][1]*q[k][1] - q[k][2]*q[k][2] + q[k][3]*q[k][3]);

        // array of rotation matrices
            matrix_A_list[k][0] = matrix_A[0][0];
            matrix_A_list[k][1] = matrix_A[0][1];
            matrix_A_list[k][2] = matrix_A[0][2];

            matrix_A_list[k][3] = matrix_A[1][0];
            matrix_A_list[k][4] = matrix_A[1][1];
            matrix_A_list[k][5] = matrix_A[1][2];

            matrix_A_list[k][6] = matrix_A[2][0];
            matrix_A_list[k][7] = matrix_A[2][1];
            matrix_A_list[k][8] = matrix_A[2][2];
    }

}

void rotate() {                             // actual rotation around of r_lab in centre and rotation around the matrix of the molecule, then push back
        for (int k = 0; k < N; k++)         // need to reset sites around COM before rotation!!! #neverForget 
            for (int j = 0; j < 5; j++) {
// rotate C,H1,H2,H3,H4        r_body[molecule][site][coord]    r_lab[molecule][com][coord]
                r_body[k][j][0] = matrix_A_list[k][0]*(r_lab[k][j][0] - r_lab[k][0][0]) + 
                                  matrix_A_list[k][1]*(r_lab[k][j][1] - r_lab[k][0][1]) + 
                                  matrix_A_list[k][2]*(r_lab[k][j][2] - r_lab[k][0][2]); 

                r_body[k][j][1] = matrix_A_list[k][3]*(r_lab[k][j][0] - r_lab[k][0][0]) + 
                                  matrix_A_list[k][4]*(r_lab[k][j][1] - r_lab[k][0][1]) + 
                                  matrix_A_list[k][5]*(r_lab[k][j][2] - r_lab[k][0][2]); 

                r_body[k][j][2] = matrix_A_list[k][6]*(r_lab[k][j][0] - r_lab[k][0][0]) + 
                                  matrix_A_list[k][7]*(r_lab[k][j][1] - r_lab[k][0][1]) + 
                                  matrix_A_list[k][8]*(r_lab[k][j][2] - r_lab[k][0][2]); 
        }
// push back to lab system
        for (int k = 0; k < N; k++)
            for (int j = 0; j < 5; j++) 
                for (int l = 0; l < 3; l++)
                    r_lab[k][j][l] = r_body[k][j][l] + r_lab[k][0][l] ;
}


void calculateForces() {

    for (int i = 0; i < N; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 3; k++) 
                a[i][j][k] = 0.;     // acceleration for [molecule][site][axis]


    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            for (int k = 0; k < 5; k++) {         // site of A
                for (int l = 0; l < 5; l++) {     // site of B, calculates forces between site k of molecule i and site l of molecule j
                    double rij[3] = {0, 0, 0};
                    double rSqd = 0;    // needs to be reset for every site l !
                    for (int m = 0; m < 3; m++) {
                        rij[m] = r_lab[j][l][m] - r_lab[i][k][m];
                        // closest image convetion           <---- problem weil nicht vom COM gemessen? sollte nur passieren wenn l = k = 0
                        if (std::abs(r_lab[j][0][m] - r_lab[i][0][m]) > 0.5 * lBox) {  // <<---- closest image COM only
         //               if (std::abs(rij[m] > 0.5*lBox)) {
                            if (rij[m] > 0)
                                rij[m] -= lBox;
                            else
                                rij[m] += lBox;
                        }
                        rSqd += rij[m] * rij[m];
                    }   // end of m
//                        std::cout << rSqd << " rSqd" << std::endl;
                    for (int m = 0; m < 3; m++) {                           // MASSEN VERGESSEN BEI DEN KRÄFTEN 
                        if ( k == 0 && l == 0) {                    // CC
                            double f_CC = 24 * epsilon_CC * sigma_CC_6 * (2 * sigma_CC_6 * std::pow(rSqd, -7)) - (std::pow(rSqd, -4));
                            a[i][k][m]  += rij[m] * f_CC / mC;
                            a[j][l][m]  -= rij[m] * f_CC / mC;
            //                std::cout << f_CC << " fcc" << std::endl;
                        }
                            else if ( k != 0 && l != 0 ) {          // HH
                                double f_HH = 24 * epsilon_HH * sigma_HH_6 * (2 * sigma_HH_6 * std::pow(rSqd, -7)) - (std::pow(rSqd, -4));
                                a[i][k][m]  += rij[m] * f_HH / mH;
                                a[j][l][m]  -= rij[m] * f_HH / mH;
              //                  std::cout << f_HH << " fhh" << std::endl;
                            }
                            else {                                  // CH
                                double f_CH = 24 * epsilon_CH * sigma_CH_6 *  (2 * sigma_CH_6 * std::pow(rSqd, -7)) - (std::pow(rSqd, -4));
                                a[i][k][m]  += rij[m] * f_CH / mass[k];
                                a[j][l][m]  -= rij[m] * f_CH / mass[l];
                //                std::cout << f_CH << " fch" << std::endl;
                        }
                    }   // end of m
                }   // end of l
            }   // end of k
}

void velocityVerletRotation() {    // rotation of molecule

    for (int i = 0; i < N; i++) {                // reset angular acceleration, quaternion dt, quaternion dtdt
        for (int j = 0; j < 3; j++) {
                w_t[i][j] = 0;
                q_t[i][j] = 0;
                q_tt[i][j] = 0;
        }
        for (int j = 0; j < 5; j++)             // reset torque
            for (int k = 0; k < 3; k++)
                torque[i][j][k] = 0;
    }

    for (int i = 0; i < N; i++)                 // calculate torque in body system
        for (int j = 0; j < 5; j++) {           // += because its a sum over all sites of the rigid body
                torque[i][j][0] += (r_lab[i][j][1] - r_lab[i][0][1]) * (a[i][j][2]/* - r_lab[i][j][2]*/) - (r_lab[i][j][2] - r_lab[i][0][2]) * (a[i][j][1]/* - r_lab[i][j][1]*/); 
                torque[i][j][1] += (r_lab[i][j][2] - r_lab[i][0][2]) * (a[i][j][0]/* - r_lab[i][j][0]*/) - (r_lab[i][j][0] - r_lab[i][0][0]) * (a[i][j][2]/* - r_lab[i][j][2]*/); 
                torque[i][j][2] += (r_lab[i][j][0] - r_lab[i][0][0]) * (a[i][j][1]/* - r_lab[i][j][1]*/) - (r_lab[i][j][1] - r_lab[i][0][1]) * (a[i][j][0]/* - r_lab[i][j][0]*/); 
        }

    for (int i = 0; i < N; i++)                 // calculate angular acceleration in body system
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 5; k++)
                w_t[i][j] += (1/Ixx) * torque[i][k][j];

    for (int i = 0; i < N; i++) {                // calculate quaternion dt q o w
        q_t[i][0] -= 0.5 * (q[i][1] * w[i][0] + q[i][2] * w[i][1] + q[i][3] * w[i][2]);
        q_t[i][1] += 0.5 * (q[i][0] * w[i][0] + q[i][2] * w[i][2] - q[i][3] * w[i][1]);
        q_t[i][2] += 0.5 * (q[i][0] * w[i][1] + q[i][3] * w[i][0] - q[i][1] * w[i][2]);
        q_t[i][3] += 0.5 * (q[i][0] * w[i][2] + q[i][1] * w[i][1] - q[i][2] * w[i][0]);

    }
// q_tt begin
    for (int i = 0; i < N; i++) {                // calculate quaternion dt o w for q_tt calculation!
        q_tt[i][0] -= 0.5 * (q_t[i][1] * w[i][0] + q_t[i][2] * w[i][1] + q_t[i][3] * w[i][2]);
        q_tt[i][1] += 0.5 * (q_t[i][0] * w[i][0] + q_t[i][2] * w[i][2] - q_t[i][3] * w[i][1]);
        q_tt[i][2] += 0.5 * (q_t[i][0] * w[i][1] + q_t[i][3] * w[i][0] - q_t[i][1] * w[i][2]);
        q_tt[i][3] += 0.5 * (q_t[i][0] * w[i][2] + q_t[i][1] * w[i][1] - q_t[i][2] * w[i][0]);
    }

    for (int i = 0; i < N; i++) {                // calculate quaternion q o w_t for q_tt calculation!
        q_tt[i][0] -= 0.5 * (q[i][1] * w_t[i][0] + q[i][2] * w_t[i][1] + q[i][3] * w_t[i][2]);
        q_tt[i][1] += 0.5 * (q[i][0] * w_t[i][0] + q[i][2] * w_t[i][2] - q[i][3] * w_t[i][1]);
        q_tt[i][2] += 0.5 * (q[i][0] * w_t[i][1] + q[i][3] * w_t[i][0] - q[i][1] * w_t[i][2]);
        q_tt[i][3] += 0.5 * (q[i][0] * w_t[i][2] + q[i][1] * w_t[i][1] - q[i][2] * w_t[i][0]);
    }
// q_tt end

// w half step
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            w[i][j] += w_t[i][j] * dt * 0.5;

// q time step
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 4; j++)
            q_d[i][j] = q[i][j] + q_t[i][j] * dt + q_tt[i][j] * dt * dt * 0.5;

// calculate quaternion normization
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 4; j++)
            q_norm += q_d[i][j]*q_d[i][j];

// unit quaternion
        for (int j = 0; j < 4; j++)       
            q_d[i][j] = q_d[i][j]/std::sqrt(q_norm);
       // std::cout << q_d[i][0]*q_d[i][0] + q_d[i][1]*q_d[i][1] + q_d[i][2]*q_d[i][2] + q_d[i][3]*q_d[i][3] << std::endl;
        q_norm = 0;
    }
// quaternion change during time step
    for (int i = 0; i < N; i++) {               // calculate quaternion q o w_t for q_tt calculation!
        q_delta[i][0] = q_d[i][0] * q[i][0] - (q_d[i][1] * q[i][1] + q_d[i][2] * q[i][2] + q_d[i][3] * q[i][3]);
        q_delta[i][1] = q_d[i][0] * q[i][1] + q_d[i][1] * q[i][0] + q_d[i][2] * q[i][3] - q_d[i][3] * q[i][2];
        q_delta[i][2] = q_d[i][0] * q[i][2] + q_d[i][2] * q[i][0] + q_d[i][3] * q[i][1] - q_d[i][1] * q[i][3];
        q_delta[i][3] = q_d[i][0] * q[i][3] + q_d[i][3] * q[i][0] + q_d[i][1] * q[i][2] - q_d[i][2] * q[i][1];
        }
// rename... whole function programmed so UGLY! :(
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 4; j++) {
            q[i][j] = q_delta[i][j];
    }

// calculation of the rotation matrix 
    for (int k = 0; k < N; k++) {
        matrix_A[0][0] = (q_delta[k][0]*q_delta[k][0] + q_delta[k][1]*q_delta[k][1] - q_delta[k][2]*q_delta[k][2] - q_delta[k][3]*q_delta[k][3]);
        matrix_A[0][1] = 2*(q_delta[k][1]*q_delta[k][2] - q_delta[k][0]*q_delta[k][3]);
        matrix_A[0][2] = 2*(q_delta[k][1]*q_delta[k][3] + q_delta[k][0]*q_delta[k][2]);

        matrix_A[1][0] = 2*(q_delta[k][1]*q_delta[k][2] + q_delta[k][0]*q_delta[k][3]);
        matrix_A[1][1] = (q_delta[k][0]*q_delta[k][0] - q_delta[k][1]*q_delta[k][1] + q_delta[k][2]*q_delta[k][2] - q_delta[k][3]*q_delta[k][3]);
        matrix_A[1][2] = 2*(q_delta[k][2]*q_delta[k][3] - q_delta[k][0]*q_delta[k][1]);

        matrix_A[2][0] = 2*(q_delta[k][1]*q_delta[k][3] - q_delta[k][0]*q_delta[k][2]);
        matrix_A[2][1] = 2*(q_delta[k][2]*q_delta[k][3] + q_delta[k][0]*q_delta[k][1]);
        matrix_A[2][2] = (q_delta[k][0]*q_delta[k][0] - q_delta[k][1]*q_delta[k][1] - q_delta[k][2]*q_delta[k][2] + q_delta[k][3]*q_delta[k][3]);

// array of rotation matrices
        matrix_A_list[k][0] = matrix_A[0][0];
        matrix_A_list[k][1] = matrix_A[0][1];
        matrix_A_list[k][2] = matrix_A[0][2];

        matrix_A_list[k][3] = matrix_A[1][0];
        matrix_A_list[k][4] = matrix_A[1][1];
        matrix_A_list[k][5] = matrix_A[1][2];

        matrix_A_list[k][6] = matrix_A[2][0];
        matrix_A_list[k][7] = matrix_A[2][1];
        matrix_A_list[k][8] = matrix_A[2][2];
    }
    rotate();
}

void velocityVerletTranslation() {          // translation of COM only!
    calculateForces();
    for (int i = 0; i < N; i ++)
        for (int j = 0; j < 5; j++) 
            for (int k = 0; k < 3; k++) {       // need to take ALL a[][][] into account! not only COM!
                r_lab[i][0][k] += v[i][0][k] * dt + a[i][j][k] * dt * dt * 0.5; // calculate new coords from COM velocity and acceleration
                // boundary conditions <--- here lies the mistake most likely <<< nope, seems legit. << not anymore.
                if (r_lab[i][0][k] < 0)    
                        r_lab[i][0][k] += lBox;
                if (r_lab[i][0][k] >= lBox)
                        r_lab[i][0][k] -= lBox;
                v[i][0][k] += 0.5 * a[i][j][k] * dt;
        }    
    setSites();
    velocityVerletRotation();         // rotates sites
    calculateForces();
    for (int i = 0; i < N; i ++) 
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 3; k++)
                v[i][0][k] += 0.5 * a[i][j][k] * dt;    // calculate new COM speed, depends on ALL forces! not only a[][][] from com
}



/*
 * FUNCTION LIST
 *  velocityVerletTranslation   --- calculates velocity verlet for translation moves for com
 *  initPositions   --- sets positions of molecules in FCC lattice  
 *  calculateForces --- calculates the forces between sites
 *  rotate  --- rotates around a matrix matrix_A to a given quaternion q
 *  initRotations   --- creates a random quaternion and its rotation matrix for all molecules
 *  setSites    --- sets sites at new com location
 */

int main() {
    initPositions();
    initRotations();
    rotate();

    std::ofstream file("methane_gas_20160229_100_256_density_rot.xyz");

    int n_step = 100;
    for (int step_number = 0; step_number < n_step; step_number++) {
        if (step_number%50 == 0)
            std::cout << step_number << std::endl;        
        file << 5*N << std::endl;
        file << "methane" << std::endl;    

        for (int k = 0; k < N; k++) {
            for (int i = 0; i < 5; i++){
                if (i == 0)
                    file << "C ";
                else
                    file << "H ";
                for (int j = 0; j < 3; j++)
                    file << r_lab[k][i][j] << " ";       
                    file << std::endl;
            }
        }
        velocityVerletTranslation();
//        velocityVerletRotation();
//        rotate();
        } 
file.close();
    return 0;
}



