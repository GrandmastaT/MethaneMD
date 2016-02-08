#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

int N = 256 ;           // number of lattice points
double rho = 0.456;        // density
double T = 1.0;         // system temperature
double dt = 0.01;
double **r;             // lattice point
double **v;             // velocities
double **a;             // accelerations
double sigma = 2.*0.1087;   // distance from C to H (in pm)
double alpha = 109.5;   // angle between HCH in grad
double hh = sqrt(2.*(sigma*sigma) * (1.-cos(alpha)));	//length between two H atoms
double mbox = hh/sqrt(2.);  // length of molecule cell
double boxL = std::pow(N/rho, 1./3.);   // calculate box length

void initialize();
void initPositions();
void initVelocities();
void rescaleVelocities();
double gasdev();


void initPositions() { //initialize FFC lattice
    int sideN = 1;

    while (4*sideN*sideN*sideN < N) // look for magic number, sideN is number of lattice points on each box axis
        sideN++;

    if (N/(4*sideN*sideN*sideN) == 1) {
        double stepL = boxL/sideN;              // step length, length between lattice points
        double hstepL = stepL/2.;

        std::cout << sideN << " number of sites" << std::endl;
        std::cout << boxL << " box length" << std::endl;
        std::cout << stepL << " step between lattice sites" << std::endl;
        std::cout << mbox << " molecule length" << std::endl;

        r = new double *[N];    // create lattice points
        for (int i = 0; i < N; i++)
            r[i] = new double [3];

        int index = 0;
        for (int i = 0; i < 2*sideN; i++) 
            for (int j = 0; j < 2*sideN; j++) 
                for (int k = 0; k < 2*sideN; k++)
                    if ((i+j+k )%2 == 0) {
                        r[index][0] = k*hstepL;
                        r[index][1] = i*hstepL;
                        r[index][2] = j*hstepL;
                        index++;
                    }    
    }
    else
        std::cout << "oida, keine magic number" << std::endl;
}

double gasdev () {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = std::sqrt(-2.0 * std::log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        return v2*fac;
    } else {
        available = false;
        return gset;
    }
}

void rescaleVelocities() {
    double vSqdSum = 0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vSqdSum += v[n][i] * v[n][i];
    double lambda = std::sqrt (3 * (N-1) * T / vSqdSum);
        for (int n = 0; n < N; n++)
            for (int i = 0; i < 3; i++)
                v[n][i] *= lambda;
}

void initVelocities() {

    a = new double *[N];   
    v = new double *[N];
    for (int i = 0; i < N; i++) {
        v[i] = new double [3];
        a[i] = new double [3];
    }
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = gasdev();

    double vCM[3] = {0, 0, 0};
    for (int i = 0; i < N; i++)     // keep COM of cube fixed
        for (int j = 0; j < 3; j++)
            vCM[j] += v[i][j];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++) 
            v[i][j] -= vCM[j];
        
     rescaleVelocities();
}


void computeAccelerations() {           // compute starting velocities
    for (int i = 0; i < N; i++)         // set all accelerations to zero
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;

    for (int i = 0; i < N-1; i++)       // loop over all distinct pairs
        for (int j = i+1; j < N; j++) {
            double rij[3];              // position of i relative to j
            double rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k];
                                        // closest image convention
                if (std::abs(rij[k]) > 0.5 * boxL) {
                    if (rij[k] > 0)
                        rij[k] -= boxL;
                    else
                        rij[k] += boxL;
                }
                rSqd += rij[k] * rij[k];
             }
            double f = 24 * (2 * std::pow(rSqd, -7) - std::pow(rSqd, -4));
            for (int k = 0; k < 3; k++) {
                a[i][k] += rij[k] * f;
                a[j][k] -= rij[k] * f;
            }
        }
}

void velocityVerlet(double dt) {
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 *a[i][k] * dt * dt;
                                    // periodic boundary conditions
            if (r[i][k] < 0)
                r[i][k] += boxL;
            if (r[i][k] >= boxL)
                r[i][k] -= boxL;
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            v[i][k] += 0.5 * a[i][k] * dt;
}

double instantTemp() {
    double sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}

class methane {                // class & its declarations
    public:  
    double com[3];
    double H1[3];
    double H2[3];
    double H3[3];
    double H4[3];
    double set_coord(double rx, double ry, double rz) {
       
        com[0] = rx ;       // update values for com
        com[1] = ry ;
        com[2] = rz ;

        H1[0] = rx - mbox/2.;       // update values for H1
        H1[1] = ry - mbox/2.;
        H1[2] = rz - mbox/2.;

        H2[0] = rx + mbox/2.;        // update values for H2
        H2[1] = ry - mbox/2.;
        H2[2] = rz + mbox/2.;

        H3[0] = rx - mbox/2.;         // update values for H3
        H3[1] = ry + mbox/2.;
        H3[2] = rz + mbox/2.;

        H4[0] = rx + mbox/2.;         // update values for H4
        H4[1] = ry + mbox/2.;
        H4[2] = rz - mbox/2.;

    }

    methane() // constructor! gnarhrhrharhrhraharharhrahararh 
    {
    }
};


int main()
{

    int Nrun = 500;
    methane methr[N];
    initPositions();   
    initVelocities();

    for (int i = 0; i < N; i++)     // creates array of class methane objects
        methr[i];
    std::ofstream file("test5.xyz");   
    std::ofstream fileT("temp.txt");

    for (int k = 0; k < Nrun; k++) {
        std::cout << k << std::endl;
        fileT << instantTemp() << std::endl;
        file << 5*N << std::endl;
        file << "blah" << std::endl;

        if (k % 25 == 0) 
            rescaleVelocities();

        for (int j = 0; j < N; j++)     // set methane on lattice point
            methr[j].set_coord(r[j][0], r[j][1], r[j][2]);
                                        //        methane meth;
        for (int j = 0; j < N; j++) {
            file << "C " << methr[j].com[0] << " " << methr[j].com[1] << " " << methr[j].com[2] << std::endl;
            file << "H1 " << methr[j].H1[0] << " " << methr[j].H1[1] << " " << methr[j].H1[2] << std::endl;
            file << "H2 " << methr[j].H2[0] << " " << methr[j].H2[1] << " " << methr[j].H2[2] << std::endl;
            file << "H3 " << methr[j].H3[0] << " " << methr[j].H3[1] << " " << methr[j].H3[2] << std::endl;
            file << "H4 " << methr[j].H4[0] << " " << methr[j].H4[1] << " " << methr[j].H4[2] << std::endl;
        }
        
        velocityVerlet(dt);
    }
    file.close();
    fileT.close();

    return 0;
}
