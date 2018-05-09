#include <chrono>
#include "constants.h"

namespace constants
{
    //global constants
    extern const double pi(3.14159265359);
    extern const double ln2(0.69314718056);
    
    /*################################################################*/
    //Simulation parameters
    extern const int Nt(4);					//lattice size in temporal direction
    extern const int Ns(4);					//lattice siye in spatial direction
		extern const int V(Nt*Ns);			//# of lattice sites, i.e. the volume
		extern const int max_bessel(50); //max order for Bessel functions
		extern const double m(0.0);			//bare mass
		extern const int nmeas(1e4);		//m
		extern const int nequi(1e4);
		extern const int nskip(2e1);
    extern const double theta(0.0);
    //extern double theta_bar(theta/2.0/pi);
    extern double theta_bar(0.0);
    extern double beta(0.1);
    extern double eta(0.0);
    extern double eta_bar(0.0);
    
    /*################################################################*/
    //Random number seed
    extern const unsigned seed(std::chrono::system_clock::now().time_since_epoch().count());
}
