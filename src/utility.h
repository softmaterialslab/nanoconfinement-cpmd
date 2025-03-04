// This is header file for the UTILITY class.
// This file includes standard library files and gsl functions that are utilized in the code
// This file also has useful constant parameters for the problem at hand

#ifndef _UTILITY_H
#define _UTILITY_H

#include<iostream>
#include<iomanip>
#include<fstream>
#include<math.h>
#include<vector>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

//OPENMP
#include <omp.h>
//BOOST MPI
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/program_options.hpp>

using namespace boost::program_options;
namespace mpi = boost::mpi;

//#define CHUNKSIZE 50				// for parallel implementation
//#define THREADSIZE 16

using namespace std;

extern mpi::environment env;
extern mpi::communicator world;

//MPI boundary parameters
extern unsigned int lowerBoundIons;
extern unsigned int upperBoundIons;
extern unsigned int sizFVecIons;
extern unsigned int lowerBoundMesh;
extern unsigned int upperBoundMesh;
extern unsigned int sizFVecMesh;


// constants
const double pi = 3.141593;							// Pi
const double dcut = 1.122462048;						// Cutoff distance in Lennard-Jones = 2 ^ (1/6)
const double dcut2 = 1.259921049;						// Cutoff distance in Lennard-Jones = 2 ^ (1/6)
const double lB_water = 0.714;							// Bjerrum length in water, in nanometers
const double epsilon_water = 80.1;						// Dielectric constant of water
const double room_temperature = 298;						// Room temperature in Kelvin
const double unitlength = 0.5*lB_water;						// Unit of length is this much nanometers
const double unitenergy = 1.3807 * pow(10.0,-16) * room_temperature;		// Unit of energy (thermal energy at room temoperature in CGS)
const double unitmass = 23 * 1.67 * pow(10.0, -24);				// Unit of mass (mass of sodium ion in CGS)
const double unittime = sqrt(unitmass * unitlength * pow(10.0,-7) * unitlength / unitenergy);	// Unit of time (length expressed in cms), result in seconds
const double kB = 1;								// Boltzmann constant in reduced units
const double scalefactor = epsilon_water * lB_water / unitlength;		// Reduced units lead to this scale factor for Coloumb interaction

class UTILITY 
{
  public:
    unsigned long int seed;
    const gsl_rng_type * T;
    gsl_rng * r;
    
    UTILITY() 
    {
      T = gsl_rng_default;
      r = gsl_rng_alloc(T);
      
      gsl_rng_env_setup();
//       srand((time(0)));                	// srand & time are built-in
//       unsigned long int s = random();  	// gsl_rng_uniform will eventually
//       gsl_rng_set(r,s); 		// seed the random number generator;
    }
    
    ~UTILITY() 
    {
      gsl_rng_free(r);
    }
};

#endif
