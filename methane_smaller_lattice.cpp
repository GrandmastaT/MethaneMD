#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

int N = 256 ;           // number of lattice points
double rho = 0.456;        // density
double **r;             // lattice point
double sigma = 1.087;   // distance from C to H (in pm)
double alpha = 109.5;   // angle between HCH in grad
double a = sqrt(2.*(sigma*sigma) * (1.-cos(alpha)));	//length between two H atoms
double k = a/sqrt(2.);  // length of molecule cell

void initPositions() { //initialize FFC lattice
    int sideN = 0;

    while (4*sideN*sideN*sideN < N) // look for magic number, sideN is number of lattice points on each box axis
        sideN++;

    if (N/(4*sideN*sideN*sideN) == 1) {
        double boxL = std::pow(N/rho, 1./3.);   // calculate box length
        double stepL = boxL/sideN;              // step length, length between lattice points
        double hstepL = stepL/2.;
        
        r = new double *[N];
        for (int i = 0; i < N; i++)
            r[i] = new double [3];

        int index = 0;
        for (int i = 0; i < (2*sideN)-1; i++) 
            for (int j = 0; j < (2*sideN)-1; j++) 
                for (int k = 0; k < (2*sideN)-1; k++)
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

class methane {                // class & its declarations
    public:
    
    double com[3];
    double H1[3];
    double H2[3];
    double H3[3];
    double H4[3];

    methane() // constructor! gnarhrhrharhrhraharharhrahararh 
   /* : com ({0., 0., 0.}), H1 {-dhh,d,0.}, H2 {dhh,d,0.}, H3 {0.,-d,dhh}, H4 {0.,-d,-dhh}*/
    {
        for (int i = 0; i < 3; i++)
            com[i] = 0;     // set starting values for C (com)

        H1[0] = -k/2.;       // set starting values for H1
        H1[1] = -k/2.;
        H1[2] = -k/2.;

        H2[0] = k/2.;        // set starting values for H2
        H2[1] = -k/2.;
        H2[2] = k/2.;

        H3[0] = -k/2.;         // set starting values for H3
        H3[1] = k/2.;
        H3[2] = k/2.;

        H4[0] = k/2;         // set starting values for H4
        H4[1] = k/2;
        H4[2] = -k/2;

    }

};


int main()
{
    initPositions();   
    methane methr[N];
    for (int i = 0; i < N; i++)     // creates array of class methane objects
        methr[i];
    
    for (int i = 0; i < N; i++){    // set methane on lattice point

        methr[i].com[0] += r[i][0];
        methr[i].com[1] += r[i][1];
        methr[i].com[2] += r[i][2];

        methr[i].H1[0] += r[i][0];
        methr[i].H1[1] += r[i][1];
        methr[i].H1[2] += r[i][2];

        methr[i].H2[0] += r[i][0];
        methr[i].H2[1] += r[i][1];
        methr[i].H2[2] += r[i][2];

        methr[i].H3[0] += r[i][0];
        methr[i].H3[1] += r[i][1];
        methr[i].H3[2] += r[i][2];

        methr[i].H4[0] += r[i][0];
        methr[i].H4[1] += r[i][1];
        methr[i].H4[2] += r[i][2];

        }


    std::ofstream file("test_box.xyz");
    file << 5*N << std::endl;
    file << "blah" << std::endl;

    methane meth;
    std::cout << meth.H1[0] << " H10" << std::endl;
    std::cout << meth.H1[1] << " H11" << std::endl;
    std::cout << meth.com[0] << " com0" << std::endl;
    std::cout << meth.com[1] << " com1" << std::endl;
    for (int i = 0; i < N; i++) {
        file << "C " << methr[i].com[0] << " " << methr[i].com[1] << " " << methr[i].com[2] << std::endl;
        file << "H1 " << methr[i].H1[0] << " " << methr[i].H1[1] << " " << methr[i].H1[2] << std::endl;
        file << "H2 " << methr[i].H2[0] << " " << methr[i].H2[1] << " " << methr[i].H2[2] << std::endl;
        file << "H3 " << methr[i].H3[0] << " " << methr[i].H3[1] << " " << methr[i].H3[2] << std::endl;
        file << "H4 " << methr[i].H4[0] << " " << methr[i].H4[1] << " " << methr[i].H4[2] << std::endl;
/*
        file << "H2H3 " << sqrt( (methr[i].H2[0]-methr[i].H3[0])*(methr[i].H2[0]-methr[i].H3[0]) + (methr[i].H2[1]-methr[i].H3[1])*(methr[i].H2[1]-methr[i].H3[1]) + (methr[i].H2[2]-methr[i].H3[2])*(methr[i].H2[2]-methr[i].H3[2])) << std::endl;

        file << "H2H1 " << sqrt( (methr[i].H2[0]-methr[i].H1[0])*(methr[i].H2[0]-methr[i].H1[0]) + (methr[i].H2[1]-methr[i].H1[1])*(methr[i].H2[1]-methr[i].H1[1]) + (methr[i].H2[2]-methr[i].H1[2])*(methr[i].H2[2]-methr[i].H1[2])) << std::endl;

        file << "H2H4 " << sqrt( (methr[i].H2[0]-methr[i].H4[0])*(methr[i].H2[0]-methr[i].H4[0]) + (methr[i].H2[1]-methr[i].H4[1])*(methr[i].H2[1]-methr[i].H4[1]) + (methr[i].H2[2]-methr[i].H4[2])*(methr[i].H2[2]-methr[i].H4[2])) << std::endl;
*/
        }

    file.close();

    return 0;
}
