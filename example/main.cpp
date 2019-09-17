/*
 * IPsat dipole amplitude fitted to HERA data
 *
 * Reference: H. Mäntysaari, P. Zurita,
 * Phys.Rev. D98 (2018) 036002 (arXiv:1804.05311)
 *
 * This code solves DGLAP using the fitted parametrization
 * for the initial condition, supporting both non-linear IPsat
 * and linearized IPnonsat (fitted separately)
 *
 * Questions and comments: Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include <dipoleamplitude.hpp>
#include <iostream>
#include <iomanip>
using namespace std;

// Example
int main(int argc, char* argv[])
{
    MZ_ipsat::DipoleAmplitude amplitude(MZ_ipsat::MZ_IPSAT); // use MZ_IPNONSAT for linearized
	// Tabulate DGLAP evolution, faster evaluation of the dipole
	amplitude.EnableLookupTable();    
    
    
    double minr=1.1e-6;
    double maxr=100;
    int points=50;
	cout <<" # r [1/GeV]  N(r, x=0.01, b=0 [1/GeV])   N(r, x=0.001, b=1 [1/GeV])" << endl;
    for (int i=0; i<points; i++)
    {
        double r =  minr * pow((maxr/minr), ((double)i)/((double)points));
        cout <<  std::scientific << std::setprecision(9) << r << " " << amplitude.N(r, 0.01, 0) << " " << amplitude.N(r, 0.001, 1) << endl;
    }
    
    
    
    return 0;
}
