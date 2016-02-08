#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

int N = 256 ;           // number of lattice points
double rho = 1.;        // density
double **r;             // lattice point
double sigma = 108.7;   // distance from C to H (in pm)
double alpha = 109.5;   // angle between HCH in grad
double d = sigma * cos(alpha/2);
double dhh = sigma * sin(alpha/2);

void initPositions() { //initialize FFC lattice
    int sideN = 0;

    while (4*sideN*sideN*sideN < N) // look for magic number, sideN is number of lattice points on each box axis
        sideN++;

    if (N/(4*sideN*sideN*sideN) == 1) {
        double boxL = std::pow(N/rho, 1./3.);   // calculate box length
        double stepL = boxL/sideN;              // step length, length between lattice points
        double hstepL = stepL/2;
        
        r = new double *[N];
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

class methane {                // class & its declarations
    public:
    
    double com[3];
    double H1[3];
    double H2[3];
    double H3[3];
    double H4[3];
        
    void addr ();   // function declaration

    methane() : com {0., 0., 0.}, H1 {-dhh,d,0.}, H2 {dhh,d,0.}, H3 {0.,-d,dhh}, H4 {0.,-d,-dhh}
    {

    }

};


int main()
{
    std::cout << d << " " << dhh << std::endl;
    initPositions();   
    std::ofstream file("test.txt");
    file << N << " N" << std::endl;
    file << "blah" << std::endl;
    methane meth;
    std::cout << meth.H1[0] << " H10" << std::endl;
    std::cout << meth.H1[1] << " H11" << std::endl;
    std::cout << meth.com[0] << " com0" << std::endl;
    std::cout << meth.com[1] << " com1" << std::endl;
    for (int i = 0; i < N; i++) 
        file << "C " << r[i][0] << " " << r[i][1] << " " << r[i][2] << std::endl;
    
    file.close();

    return 0;
}
