#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <gsl/gsl_rng.h>

// units: nm, ps, au
// i,j,k,l,m,n not free
const int N = 256 ;             // number of lattice points
double rho = 0.395052;          // density u/(nm³) in gas ad T = 25°C
double T = 25. + 273;           // system temperature in K
double dt = 0.01;               // time step in ps
double ***r_lab;                // lattice point
double ***v;                    // velocities
double ***a;                    // accelerations

double d_CH = 0.1087;           // distance from C to H (in nm)
double angle_HCH = 109.5;       // angle between HCH in grad
double hh = d_CH * sqrt(2.* (1.-cos(angle_HCH)));	//length between two H atoms
double mBox = hh/sqrt(2.);      // length of molecule cell
double lBox = std::pow(N/rho, 1./3.);   // calculate box length
double **q;                     //list of quaternion paramters q[0], q[1], q[2], q[3]
double q_norm;                  // normalization factor of quaternion

double **matrix_A_list;
double **matrix_A;

double ***r_body;
double sigma_HH = 0.25;     // nm
double sigma_CC = 0.35;     // nm
double sigma_CH = 0.3;      // nm
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

        for (int i = 0; i < N; i++) { // create array for lattice points
            r_lab[i] = new double *[5];
            r_body[i] = new double *[5]; 
            v[i] = new double *[5];
            a[i] = new double *[5];
            for (int j = 0; j < 5; j++) {            
                r_lab[i][j] = new double [3];
                r_body[i][j] = new double [3]; 
                v[i][j] = new double [3];
                a[i][j] = new double [3];
            }

        }

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

    matrix_A = new double *[3];             // matrices for rotations
    for (int k = 0; k < 3; k++)
        matrix_A[k] = new double[3];

    matrix_A_list = new double *[N];        // list of rotation matrices, [molecule][A_ij]
    for (int k = 0; k < N ;k++)
        matrix_A_list[k] = new double [9];

    for (int k = 0; k < N; k++) {           // for every molecule   
        q_norm = 0;
        q[k] = new double [4];
        for (int i = 0; i < 4; i++)             // generate random quaternion
            q[k][i] = (0.5 - gsl_rng_uniform(Seed2)); 

        for (int i = 0; i < 4; i++)             // calculate quaternion randomization
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
                a[i][j][k] = 0;     // acceleration for [molecule][site][axis]

    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            for (int k = 0; k < 5; k++) {         // site of A
                double rij[5][3];
                for (int l = 0; l < 5; l++) {     // site of B, calculates forces between site k of molecule i and site l of molecule j
                    double rSqd = 0;    // needs to be reset for every site l !
                    for (int m = 0; m < 3; m++) {
                        rij[l][m] = r_lab[i][k][m] - r_lab[j][l][m];
                        if (std::abs(rij[l][m] > 0.5*lBox))
                            if (rij[l][m] > 0)
                                rij[l][m] -= lBox;
                            else
                                rij[l][m] += lBox;

                        rSqd += rij[l][m] * rij[l][m];
                    }   // end of m
                    for (int m = 0; m < 3; m++) {   
                        if ( k == 0 && l == 0) {                    // CC
                            double f_CC = 24 * epsilon_CC * (2 * (std::pow(sigma_CC, 12) / std::pow(rSqd, -14)) - std::pow(sigma_CC, 6) / (std::pow(rSqd, -8)));
                            a[i][k][m]  += rij[l][m] * f_CC;
                            a[j][k][m]  -= rij[l][m] * f_CC;
                        }
                            else if ( k != 0 && l != 0 ) {          // HH
                                double f_HH = 24 * epsilon_HH * (2 * (std::pow(sigma_HH, 12) / std::pow(rSqd, -14)) - std::pow(sigma_HH, 6) / (std::pow(rSqd, -8)));
                                a[i][k][m]  += rij[l][m] * f_HH;
                                a[j][k][m]  -= rij[l][m] * f_HH;
                            }
                            else {                                  // CH
                                double f_CH = 24 * epsilon_CH * (2 * (std::pow(sigma_CH, 12) / std::pow(rSqd, -14)) - std::pow(sigma_CH, 6) / (std::pow(rSqd, -8)));
                                a[i][k][m]  += rij[l][m] * f_CH;
                                a[j][k][m]  -= rij[l][m] * f_CH;
                        }
                    }   // end of m
                }   // end of l
            }   // end of k
}

/*
void calculateForces() {            // calculate LJ forces between all sites
    for (int i = 0; i < N; i++)
        for(int j = 0; j < 5; j++)
            for (int k = 0; k < 3; k++)
                a[i][j][k] = 0;
    for (int i = 0; i < N-1; i++)                   // molecule i        
        for (int j = i+1; j < N; j++) {             // molecule j
            double rij[5][3];                       // rij[sites][coords]
            for (int k = 0; k < 5; k++)             // sites of i 
                for (int l = 0; l < 5; l++) {       // sites of j
                    double rSqd = 0.;               //set rSqd = 0 for site k on site l
                    for (int m = 0; m < 3; m++) {   // coords
                        rij[k][m] = r_lab[i][k][m] - r_lab[j][l][m];    //distance between sites
                        if (abs(rij[0][m]) > 0.5 * lBox) {              // closest image convention
                            if (rij[k][m] > 0)
                                for (int var = 0; var < 5; var++)
                                    rij[var][m] -= lBox;
                            else
                                for (int var = 0; var < 5; var++)
                                    rij[var][m] += lBox;
                        }
                        rSqd += rij[k][m] * rij[k][m];
                    }
                    for (int m = 0; m < 3; m++) {   // coords
                        if ( k == 0 && l == 0) {                    // CC
                            double f_CC = 24 * ((lj_CC_12 * std::pow(rSqd, -7)) - (lj_CC_6 * std::pow(rSqd, -4)));
                            a[i][k][m]  += rij[k][m] * f_CC;
                            a[j][k][m]  -= rij[k][m] * f_CC;
                        }
                            else if ( k != 0 && l != 0 ) {          // HH
                                double f_HH = 24 * ((lj_HH_12 * std::pow(rSqd, -7)) - (lj_HH_6 * std::pow(rSqd, -4)));
                                a[i][k][m]  += rij[k][m] * f_HH;
                                a[j][k][m]  -= rij[k][m] * f_HH;
                            }
                            else {                                  // CH
                                double f_CH = 24 * ((lj_CH_12 * std::pow(rSqd, -7)) - (lj_CH_6 * std::pow(rSqd, -4)));
                                a[i][k][m]  += rij[k][m] * f_CH;
                                a[j][k][m]  -= rij[k][m] * f_CH;
                        }
                    }
            }
        }   
}
*/

void velocityVerletTranslation() {          // translation of COM only!
    calculateForces();
    for (int i = 0; i < N; i ++) 
        for (int k = 0; k < 3; k++) {
            r_lab[i][0][k] += v[i][0][k] * dt + a[i][0][k] * dt * dt * 0.5; // calculate new coords from COM velocity and acceleration
            if (r_lab[i][0][k] < 0)     // boundary conditions <--- here lies the mistake most likely
                    r_lab[i][0][k] += lBox;
            if (r_lab[i][0][k] >= lBox)
                    r_lab[i][0][k] -= lBox;
            v[i][0][k] += 0.5 * a[i][0][k] * dt;
        }    
    setSites();
    calculateForces();
    for (int i = 0; i < N; i ++) 
        for (int k = 0; k < 3; k++)
            v[i][0][k] += 0.5 * a[i][0][k] * dt;    // calculate new COM speed
}

void velocityVerletRotation() {             // rotation of molecule
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

    std::ofstream file("methane_gas_20160225_200_256.xyz");
    std::ofstream file_test("control_file.txt");
    int n_step = 5;
    for (int step_number = 0; step_number < n_step; step_number++) {
        if (step_number%50 == 0)
            std::cout << step_number << std::endl;        
/*
        for  (int k = 0; k < N; k++){ 
            double norm_test = 0;
            for (int i = 0; i < 4; i++) {
                std::cout << q[k][i] << std::endl;
                norm_test += q[k][i]*q[k][i];
            }           
            std::cout << norm_test << " norm"  << std::endl;
        }
*/

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
        for (int k = 0; k < N; k++){
            file_test << step_number << std::endl;
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 3; j++) {
                    file_test << r_lab[k][i][j] << " r ";
                    file_test << a[k][i][j] << " a ";
                }
                file_test << i << " site " << std::endl;
            }
        }
        velocityVerletTranslation();

        for (int k = 0; k < N; k++) {
            file_test << " after velocity verlet " << std::endl;
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 3; j++) {
                    file_test << r_lab[k][i][j] << " r " ;
                    file_test << a[k][i][j] << " a " ;
                }
                file_test << i << " site " << std::endl;
            }
        }
        rotate();
        
        for (int k = 0; k < N; k++) {
            file_test << " after rotate " << std::endl;
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 3; j++) {
                    file_test << r_lab[k][i][j] << " r ";
                    file_test << a[k][i][j] << " a ";
                }
                file_test << i << " site " << std::endl;
            }
        }
    }
    file.close();
//    std::cout << "test " << 0%3 << std::endl;
    return 0;
}



