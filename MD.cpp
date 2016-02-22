#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <gsl/gsl_rng.h>


// i,j,k,l,m,n not free
const int N = 256 ;             // number of lattice points
double rho = 0.456;             // density
double T = 1.0;                 // system temperature
double dt = 0.01;               // time step
double **r;                     // lattice point
double **v;                     // velocities
double **a;                     // accelerations
double d_CH = 2.*0.1087;        // distance from C to H (in pm)
double angle_HCH = 109.5;       // angle between HCH in grad
double hh = sqrt(2.*(d_CH*d_CH) * (1.-cos(angle_HCH)));	//length between two H atoms
double mBox = hh/sqrt(2.);      // length of molecule cell
double lBox = std::pow(N/rho, 1./3.);   // calculate box length
double ang_v;
double ang_L;
double **q;                     //list of quaternion paramters q[0], q[1], q[2], q[3]
double q_norm;                  // normalization factor of quaternion
double **r_rot;

class methane {                // class & its declarations
    public:  
    double com[3];
    double H1[3];
    double H2[3];
    double H3[3];
    double H4[3];
    double ** meth_coord;

    double set_coord(double rx, double ry, double rz) {

    meth_coord = new double* [5];
    for (int i = 0; i < 5; i++)
        meth_coord[i] = new double [3];

    meth_coord[0][0] = rx ;       // update values for com
    meth_coord[0][1] = ry ;
    meth_coord[0][2] = rz ;

    meth_coord[1][0] = rx - mBox/2.;       // update values for H1
    meth_coord[1][1] = ry - mBox/2.;
    meth_coord[1][2] = rz - mBox/2.;

    meth_coord[2][0] = rx + mBox/2.;        // update values for H2
    meth_coord[2][1] = ry - mBox/2.;
    meth_coord[2][2] = rz + mBox/2.;

    meth_coord[3][0] = rx - mBox/2.;         // update values for H3
    meth_coord[3][1] = ry + mBox/2.;
    meth_coord[3][2] = rz + mBox/2.;

    meth_coord[4][0] = rx + mBox/2.;         // update values for H4
    meth_coord[4][1] = ry + mBox/2.;
    meth_coord[4][2] = rz - mBox/2.;

    }
};

methane molecules[N];

/*
class rotation_q {
    double matrix_A[3][3];

    double set_q (double quat[4]) {    
        matrix_A[0][0] = (quat[0]*quat[0] + quat[1]*quat[1] - quat[2]*quat[2] - quat[3]*quat[3]);
        matrix_A[0][1] = 2*(quat[1]*quat[2] - quat[0]*quat[3]);
        matrix_A[0][2] = 2*(quat[1]*quat[3] + quat[0]*quat[2]);

        matrix_A[1][1] = 2*(quat[1]*quat[2] + quat[0]*quat[3]);
        matrix_A[1][0] = (quat[0]*quat[0] - quat[1]*quat[1] + quat[2]*quat[2] - quat[3]*quat[3]);
        matrix_A[1][2] = 2*(quat[2]*quat[3] - quat[0]*quat[1]);

        matrix_A[2][1] = 2*(quat[1]*quat[3] - quat[0]*quat[2]);
        matrix_A[2][2] = 2*(quat[2]*quat[3] + quat[0]*quat[1]);
        matrix_A[2][0] = (quat[0]*quat[0] - quat[1]*quat[1] - quat[2]*quat[2] + quat[3]*quat[3]);
    }

};
*/

void initPositions() { //initialize FFC lattice
    int axisN = 1;

    while (4*axisN*axisN*axisN < N) // look for magic number, sideN is number of lattice points on each box axis
        axisN++;

    if (N/(4*axisN*axisN*axisN) == 1) {
        double stepL = lBox/axisN;              // step length, length between lattice points
        double hstepL = stepL/2.;

        r = new double *[N];    // create lattice points
        r_rot = new double *[N]; // create rotation points

        for (int i = 0; i < N; i++) // create array for lattice points
            r[i] = new double [3];

        for (int i = 0; i < N; i++)    // create array for rotated points
            r_rot[i] = new double [12];  // 0-2 H1, 3-5 H2, 6-8 H3, 9-11 H4

        int index = 0;
        for (int i = 0; i < 2*axisN; i++) 
            for (int j = 0; j < 2*axisN; j++) 
                for (int k = 0; k < 2*axisN; k++)
                    if ((i+j+k)%2 == 0) {
                        r[index][0] = k*hstepL;
                        r[index][1] = i*hstepL;
                        r[index][2] = j*hstepL;
                        index++;
                    }    
    }
    else
        std::cout << "oida, keine magic number" << std::endl;
}


void initRotations() {  //initialize molecule rotation
    const gsl_rng_type * Seed1;
    gsl_rng * Seed2;
    gsl_rng_env_setup();

    Seed1 = gsl_rng_default;
    Seed2 = gsl_rng_alloc (Seed1);
    
    
    q = new double *[N];
    double matrix_A[3][3];
    double matrix_A_list[N][9];

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
        for (int k = 0; k < N; k++) {
// rotate H1
            r_rot[k][0] = matrix_A_list[k][0]*(molecules[k].meth_coord[1][0] - r[k][0]) 
                        + matrix_A_list[k][1]*(molecules[k].meth_coord[1][1] - r[k][1])  
                        + matrix_A_list[k][2]*(molecules[k].meth_coord[1][2] - r[k][2]);            
            r_rot[k][1] = matrix_A_list[k][3]*(molecules[k].meth_coord[1][0] - r[k][0]) 
                        + matrix_A_list[k][4]*(molecules[k].meth_coord[1][1] - r[k][1]) 
                        + matrix_A_list[k][5]*(molecules[k].meth_coord[1][2] - r[k][2]);
            r_rot[k][2] = matrix_A_list[k][6]*(molecules[k].meth_coord[1][0] - r[k][0]) 
                        + matrix_A_list[k][7]*(molecules[k].meth_coord[1][1] - r[k][1]) 
                        + matrix_A_list[k][8]*(molecules[k].meth_coord[1][2] - r[k][2]);
// rotate H2
            r_rot[k][3] = matrix_A_list[k][0]*(molecules[k].meth_coord[2][0] - r[k][0]) 
                        + matrix_A_list[k][1]*(molecules[k].meth_coord[2][1] - r[k][1]) 
                        + matrix_A_list[k][2]*(molecules[k].meth_coord[2][2] - r[k][2]);            
            r_rot[k][4] = matrix_A_list[k][3]*(molecules[k].meth_coord[2][0] - r[k][0]) 
                        + matrix_A_list[k][4]*(molecules[k].meth_coord[2][1] - r[k][1]) 
                        + matrix_A_list[k][5]*(molecules[k].meth_coord[2][2] - r[k][2]);
            r_rot[k][5] = matrix_A_list[k][6]*(molecules[k].meth_coord[2][0] - r[k][0]) 
                        + matrix_A_list[k][7]*(molecules[k].meth_coord[2][1] - r[k][1]) 
                        + matrix_A_list[k][8]*(molecules[k].meth_coord[2][2] - r[k][2]);
// rotate H3
            r_rot[k][6] = matrix_A_list[k][0]*(molecules[k].meth_coord[3][0] - r[k][0]) 
                        + matrix_A_list[k][1]*(molecules[k].meth_coord[3][1] - r[k][1]) 
                        + matrix_A_list[k][2]*(molecules[k].meth_coord[3][2] - r[k][2]);
            r_rot[k][7] = matrix_A_list[k][3]*(molecules[k].meth_coord[3][0] - r[k][0]) 
                        + matrix_A_list[k][4]*(molecules[k].meth_coord[3][1] - r[k][1]) 
                        + matrix_A_list[k][5]*(molecules[k].meth_coord[3][2] - r[k][2]);
            r_rot[k][8] = matrix_A_list[k][6]*(molecules[k].meth_coord[3][0] - r[k][0]) 
                        + matrix_A_list[k][7]*(molecules[k].meth_coord[3][1] - r[k][1])  
                        + matrix_A_list[k][8]*(molecules[k].meth_coord[3][2] - r[k][2]);
// rotate H4
            r_rot[k][9] = matrix_A_list[k][0]*(molecules[k].meth_coord[4][0] - r[k][0]) 
                        + matrix_A_list[k][1]*(molecules[k].meth_coord[4][1] - r[k][1]) 
                        + matrix_A_list[k][2]*(molecules[k].meth_coord[4][2] - r[k][2]);            
            r_rot[k][10] = matrix_A_list[k][3]*(molecules[k].meth_coord[4][0] - r[k][0]) 
                        + matrix_A_list[k][4]*(molecules[k].meth_coord[4][1] - r[k][1]) 
                        + matrix_A_list[k][5]*(molecules[k].meth_coord[4][2] - r[k][2]);
            r_rot[k][11] = matrix_A_list[k][6]*(molecules[k].meth_coord[4][0] - r[k][0]) 
                        + matrix_A_list[k][7]*(molecules[k].meth_coord[4][1] - r[k][1]) 
                        + matrix_A_list[k][8]*(molecules[k].meth_coord[4][2] - r[k][2]);
        }


        for (int k = 0; k < N; k++)
            for (int j = 0; j < 12; j=j+3) {
                r_rot[k][j] += r[k][0] ;
                r_rot[k][j+1] += r[k][1] ;
                r_rot[k][j+2] += r[k][2] ;
            }
}

int main() {
    initPositions();

    for (int k = 0; k < N; k++)
        molecules[k].set_coord(r[k][0], r[k][1], r[k][2]);

    initRotations();

    for  (int k = 0; k < N; k++){ 
        double norm_test = 0;
        for (int i = 0; i < 4; i++) {
            std::cout << q[k][i] << std::endl;
            norm_test += q[k][i]*q[k][i];
        }           
        std::cout << norm_test << " norm"  << std::endl;
    }
    std::ofstream file_r("position_test.xyz");
    std::ofstream file("methane_20160220_1_256.xyz");

    file << 5*N << std::endl;
    file << "methane" << std::endl;    

    file_r << 5*N << std::endl;
    file_r << "test for lattice structure" << std::endl;

    for (int k = 0; k < N; k++) {
        for (int j = 0; j < 5; j++){
            if (k == 0){
                file_r << "C ";
                for (int l = 0; l < 3; l++)
                    file_r << " " << molecules[k].meth_coord[j][l];
            }
            else
                file_r << "H ";
                for (int l = 0; l < 3; l++)
                    file_r << " " << molecules[k].meth_coord[j][l];
            file_r << std::endl;
        }


        file << "O " << molecules[k].meth_coord[0][0] << " " << molecules[k].meth_coord[0][1] << " " << molecules[k].meth_coord[0][2]; 
        for (int j = 0; j < 12; j++) {
            if (j%3 == 0) {
                file << std::endl;
                file << "H ";
                file << r_rot[k][j] << " ";
            }
            else
                file << r_rot[k][j] << " ";
        }
        file << std::endl;
    }
    file.close();
    file_r.close();
//    std::cout << "test " << 0%3 << std::endl;
    return 0;
}



