// Last update: April 5, 2013
// This is main.
// This is MD simulation of an electrolyte confined within planar walls
// Problem : Compute density profiles of ions trapped within planar walls
/* Useful studies :	
		     1. Role of dielectric contrast
		     2. Role of valency of ions
		     3. Role of varying salt concentration
// To implement:
		     1. Charged walls
		     
*/
/* lessons learned :
		     1. do not use exit(1) to exit the program. this leads to memory leaks.

*/

#include <boost/program_options.hpp>
#include <time.h>
#include "utility.h"
#include "interface.h"
#include "particle.h"
#include "vertex.h"
#include "BIN.h"
#include "control.h"
#include "functions.h"
#include "precalculations.h"
#include "thermostat.h"
#include "fmd.h"

//MPI boundary parameters
unsigned int lowerBoundIons;
unsigned int upperBoundIons;
unsigned int sizFVecIons;
unsigned int lowerBoundMesh;
unsigned int upperBoundMesh;
unsigned int sizFVecMesh;

mpi::environment env;
mpi::communicator world;

void
cpmd(vector<PARTICLE> &, INTERFACE &, vector<THERMOSTAT> &, vector<THERMOSTAT> &, vector<BIN> &, CONTROL &, CONTROL &);
// removed particle colloid, fake thermostat and fmdremote

int main(int argc, char *argv[]) {

    // Electrostatic system variables
    double bx, by, bz;        // lengths of the box
    double ein;            // permittivity of inside medium
    double e_left;                // permittivity of leftside medium
    double e_right;               // permittivity of rightside medium
    int pz_in;            // positive valency of ions inside
    int nz_in;            // negative valency of ions inside
    double salt_conc_in;        // salt concentration outside	(enter in M)
    double saltion_diameter_in;    // inside salt ion diameter	(positive and negative ions assumed to have same diameter at this point)
    double T, fake_T;            // temperature at which the system of ions is

    // Simulation related variables
    double fraction_diameter;        // fraction that multiplies the diameter to generate the discretization width for the interface
    double Q, fake_Q;                // thermostat mass required to generate canonical ensemble
    unsigned int chain_length_real, chain_length_fake;
    double bin_width;            // width of the bins used to compute density profiles
    CONTROL fmdremote, cpmdremote;            // remote control for cpmd

    // Different parts of the system
    vector<PARTICLE> saltion_in;        // salt ions inside
    vector<PARTICLE> ion;        // all ions in the system
    INTERFACE box;            // interface, z planes hard walls; rest periodic boundaries

    // Analysis
    vector<BIN> bin;            // bins

    // Fake system
    double fakemass;

    // Get input values from the user
    options_description desc("Usage:\nrandom_mesh <options>");
    desc.add_options()
            ("help,h", "print usage message")
            ("bx,X", value<double>(&bx)->default_value(22.848),
             "box length in x direction in nanometers")                    // enter in nanometers
            ("by,Y", value<double>(&by)->default_value(22.848),
             "box length in y direction in nanometers")                    // enter in nanometers
            ("bz,Z", value<double>(&bz)->default_value(3),
             "box length in z direction in nanometers")                    // enter in nanometers
            ("epsilon_in,e", value<double>(&ein)->default_value(80.1), "dielectric const inside")
            ("epsilon_left,l", value<double>(&e_left)->default_value(20), "dielectric const leftside")
            ("epsilon_right,r", value<double>(&e_right)->default_value(80.1), "dielectric const rightside")
            ("pz_in,p", value<int>(&pz_in)->default_value(1), "positive valency inside")
            ("nz_in,n", value<int>(&nz_in)->default_value(-1), "negative valency inside")
            ("salt_conc_in,c", value<double>(&salt_conc_in)->default_value(0.1), "salt concentration inside")
            ("saltion_diameter_in,d", value<double>(&saltion_diameter_in)->default_value(0.714),
             "salt ion diameter inside")        // enter in nanometers
            ("fraction_diameter,g", value<double>(&fraction_diameter)->default_value(1.4),
             "for interface discretization width")
            ("thermostat_mass,Q", value<double>(&Q)->default_value(1.0), "thermostat mass")
            ("chain_length_real,L", value<unsigned int>(&chain_length_real)->default_value(5),
             "chain length for real system: enter L+1 if you want L thermostats")
            ("fake_temperature,k", value<double>(&fake_T)->default_value(0.002), "fake temperature")
            ("fake_thermostat_mass,q", value<double>(&fake_Q)->default_value(1.0), "fake thermostat mass")
            ("chain_length_fake,D", value<unsigned int>(&chain_length_fake)->default_value(5),
             "chain length for fake system: enter L+1 if you want L thermostats")
            ("fakemass,f", value<double>(&fakemass)->default_value(100.0), "fake mass")
            ("cpmd_fake_mass,M", value<double>(&cpmdremote.fakemass)->default_value(1.0), "cpmd fake mass")
            ("bin_width,B", value<double>(&bin_width)->default_value(0.1), "bin width")
            ("cpmd_timestep,T", value<double>(&cpmdremote.timestep)->default_value(0.0005), "time step used in cpmd")
            ("fmd_steps,s", value<int>(&fmdremote.steps)->default_value(200), "steps used in fmd")
            ("cpmd_steps,S", value<int>(&cpmdremote.steps)->default_value(1000), "steps used in cpmd")
            ("cpmd_eqm,P", value<int>(&cpmdremote.hiteqm)->default_value(10), "production begin (cpmd)")
            ("cpmd_freq,F", value<int>(&cpmdremote.freq)->default_value(10), "sample frequency (cpmd)")
            ("cpmd_verify,y", value<int>(&cpmdremote.verify)->default_value(1000), "verify (cpmd)")
            ("cpmd_extra_compute,x", value<int>(&cpmdremote.extra_compute)->default_value(10),
             "compute additional (cpmd)")
            ("cpmd_writedensity,w", value<int>(&cpmdremote.writedensity)->default_value(100), "write density files")
            ("cpmd_movie_freq,m", value<int>(&cpmdremote.moviefreq)->default_value(100), "compute additional (cpmd)");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    // Set up the system
    T = 1;        // set temperature
    box.ein = ein;
    box.e_left = e_left;
    box.e_right = e_right;
    box.set_up(salt_conc_in, 0, 0, pz_in, 0, 0, bx / unitlength, by / unitlength, bz / unitlength);
    box.put_saltions_inside(saltion_in, pz_in, nz_in, salt_conc_in, saltion_diameter_in, ion);
    // box.put_one_saltion_in_center(saltion_in, pz_in, nz_in, saltion_diameter_in, ion);						//put one ion in center

    make_bins(bin, box, bin_width);    // set up bins to be used for computing density profiles
    vector<double> initial_density;
    bin_ions(ion, box, initial_density, bin);    // bin the ions to get initial density profile

//   box.number_of_vertices = total_gridpoints;
    box.discretize(saltion_diameter_in / unitlength, fraction_diameter);

    fmdremote.fakemass = fakemass;

    // output to screen the parameters of the problem
    if (world.rank() == 0) {
        cout << "\nProgram starts\n";
        cout << "Dielectric constant of water is " << epsilon_water << endl;
        cout << "Reduced units: scalefactor entering in Coloumb interaction is " << scalefactor << endl;
        cout << "Box dimensions x | y | z " << setw(15) << box.lx << setw(15) << box.ly << setw(15) << box.lz << endl;
        cout << "Permittivity inside " << box.ein << endl;
        cout << "Permittivity leftside " << box.e_left << endl;
        cout << "Permittivity rightside " << box.e_right << endl;
        cout << "Contrast strength left " << 2 * (box.e_left - box.ein) / (box.e_left + box.ein) << endl;
        cout << "Contrast strength right " << 2 * (box.e_right - box.ein) / (box.e_right + box.ein) << endl;
        cout << "Positive ion valency inside " << pz_in << endl;
        cout << "Negative ion valency inside " << nz_in << endl;
        cout << "Salt ion diameter inside " << saltion_diameter_in / unitlength << endl;
        cout << "Salt concentration inside " << salt_conc_in << endl;
        cout << "Debye length inside " << box.inv_kappa_in << endl;
        cout << "Mean separation inside " << box.mean_sep_in << endl;
        cout << "Number of salt ions inside " << saltion_in.size() << endl;
        cout << "Temperature " << T << endl;
        cout << "Binning width (uniform) " << bin[0].width << endl;
        cout << "Number of bins " << bin.size() << endl;
        cout << "Number of points discretizing the left and right z planes " << box.leftplane.size() << setw(10)
             << box.rightplane.size() << endl;
    }

    int numOfNodes = world.size();
    if (world.rank() == 0) {
#pragma omp parallel default(shared)
        {
            if (omp_get_thread_num() == 0) {
                printf("The app comes with MPI and OpenMP (Hybrid) parallelization)\n");
                printf("Number of MPI processes used %d\n", numOfNodes);
                printf("Number of OpenMP threads per MPI process %d\n", omp_get_num_threads());
                printf("Make sure that number of grid points / ions is greater than %d\n",
                       omp_get_num_threads() * numOfNodes);
            }
        }
    }

    // write to files

    // initial configuration
    if (world.rank() == 0) {
        ofstream initial_configuration("outfiles/initialconfig.dat");
        for (unsigned int i = 0; i < ion.size(); i++)
            initial_configuration << "ion" << setw(5) << ion[i].id << setw(15) << "charge" << setw(5) << ion[i].q
                                  << setw(15) << "position" << setw(15) << ion[i].posvec << endl;
        initial_configuration.close();

        // initial density
        ofstream density_profile("outfiles/initial_density_profile.dat", ios::out);
        for (unsigned int b = 0; b < initial_density.size(); b++)
            density_profile << bin[b].lower << setw(15) << initial_density.at(b) << endl;
        density_profile.close();
    }
    // check point
    double totalions = 0;
    for (unsigned int b = 0; b < initial_density.size(); b++)
        totalions += initial_density.at(b) * bin[b].volume;
    if (world.rank() == 0)
        cout << "total ions " << totalions << "  should be " << ion.size() << endl;
    int totalpions = 0, totalnions = 0;
    for (unsigned int i = 0; i < ion.size(); i++) {
        if (ion[i].valency > 0)
            totalpions += 1;
        else if (ion[i].valency < 0)
            totalnions += 1;
    }

    // some calculations before simulation begins
    double charge = box.total_charge_inside(ion);

    if (world.rank() == 0) {
        cout << "Total positive ions " << totalpions << endl;
        cout << "Total negative ions " << totalnions << endl;
        cout << "Total charge inside" << charge << endl;
    }
    // prepare for cpmd : make real baths
    vector<THERMOSTAT> real_bath;
    if (chain_length_real == 1)
        real_bath.push_back(THERMOSTAT(0, T, 3 * ion.size(), 0.0, 0, 0));
    else {
        real_bath.push_back(THERMOSTAT(Q, T, 3 * ion.size(), 0, 0, 0));
        while (real_bath.size() != chain_length_real - 1)
            real_bath.push_back(THERMOSTAT(Q / (3 * ion.size()), T, 1, 0, 0, 0));
        real_bath.push_back(THERMOSTAT(0, T, 3 * ion.size(), 0.0, 0, 0));
        // final bath is dummy bath (dummy bath always has zero mass)
    }

    vector<THERMOSTAT> fake_bath;
    if (chain_length_fake == 1)
        fake_bath.push_back(THERMOSTAT(0, fake_T, box.leftplane.size(), 0.0, 0, 0));
    else {
        fake_bath.push_back(THERMOSTAT(fake_Q, fake_T, box.leftplane.size(), 0, 0, 0));
        while (fake_bath.size() != chain_length_fake - 1)
            fake_bath.push_back(THERMOSTAT(fake_Q / box.leftplane.size(), fake_T, 1, 0, 0, 0));
        fake_bath.push_back(THERMOSTAT(0, fake_T, box.leftplane.size(), 0.0, 0,
                                       0));            // finally, the coding trick: dummy bath (dummy bath always has zero mass)
    }

    // Induced charge system
    // fake positions initialized
    for (unsigned int i = 0; i < box.leftplane.size(); i++)
        box.leftplane[i].w = 0.0;
    double total_ind_charge_left = 0;
    for (unsigned int i = 0; i < box.leftplane.size(); i++)
        total_ind_charge_left += box.leftplane[i].w * box.leftplane[i].a;

    for (unsigned int i = 0; i < box.rightplane.size(); i++)
        box.rightplane[i].w = 0.0;
    double total_ind_charge_right = 0;
    for (unsigned int i = 0; i < box.rightplane.size(); i++)
        total_ind_charge_right += box.rightplane[i].w * box.rightplane[i].a;

    if (world.rank() == 0) {
        cout << "total induced left interface" << total_ind_charge_left << endl;
        cout << "total induced right interface" << total_ind_charge_right << endl;
    }

    clock_t t;
    t = clock();
    precalculate(box); // calculate the Greens function within interface or between interface
    t = clock() - t;
    //exact incuded charge density calculation
    //vector<VERTEX> exact_induced_density_left;
    //vector<VERTEX> exact_induced_density_right;

    //exact_induced_density_left = box.compute_exact_induced_density(ion,'l');
    //exact_induced_density_right = box.compute_exact_induced_density(ion,'r');
    //double total_energy = box.exact_total_energy(exact_induced_density_left, exact_induced_density_right, ion);
    //cout << "total energy should be " << total_energy << endl;
    if (world.rank() == 0)
        cout << "pre_calculation finished,  time taken: " << ((double) t) / CLOCKS_PER_SEC << " seconds" << endl;

    // Fictitious Molecular Dynamics

    t = clock();
    fmdremote.timestep = 0.0005;
    fmdremote.verify = 0;
    fmdremote.anneal = 'y';
    fmd(ion, box, fmdremote, cpmdremote);
    t = clock() - t;
    if (world.rank() == 0)
        cout << "fmd calculation finished,  time taken: " << ((double) t) / CLOCKS_PER_SEC << " seconds" << endl;
/*
  ofstream induced_density_left ("outfiles_windows/induced_density_left.dat");
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    induced_density_left << k+1 <<setw(15) << box.leftplane[k].posvec << setw(15) << 
    sqrt(box.leftplane[k].posvec.x * box.leftplane[k].posvec.x + box.leftplane[k].posvec.y * box.leftplane[k].posvec.y) << setw(15) << 
    box.leftplane[k].w << setw(15) << box.leftplane[k].wmean << setw(15)<< exact_induced_density_left[k].w << endl;		// NOTE plotting also the cartesian coordinates
  induced_density_left.close();
  ofstream induced_density_right("outfiles_windows/induced_density_right.dat");
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    induced_density_right << k+1 <<setw(15) << box.rightplane[k].posvec << setw(15) <<
    sqrt(box.rightplane[k].posvec.x * box.rightplane[k].posvec.x + box.rightplane[k].posvec.y * box.rightplane[k].posvec.y) << setw(15) <<
    box.rightplane[k].w << setw(15) << box.rightplane[k].wmean << setw(15)<< exact_induced_density_right[k].w << endl;		// NOTE plotting also the cartesian coordinates
  induced_density_right.close();
 */



    //MPI Boundary calculations for ions
    unsigned int rangeIons = (ion.size() + world.size() - 1) / (1.0 * world.size());
    lowerBoundIons = world.rank() * rangeIons;
    upperBoundIons = (world.rank() + 1) * rangeIons - 1;
    sizFVecIons = rangeIons;
    if (world.rank() == world.size() - 1) {
        upperBoundIons = ion.size() - 1;
    }

    //MPI Boundary calculations for meshPoints
    unsigned int rangeMesh = (box.leftplane.size() + world.size() - 1) / (1.0 * world.size());
    lowerBoundMesh = world.rank() * rangeMesh;
    upperBoundMesh = (world.rank() + 1) * rangeMesh - 1;
    sizFVecMesh = rangeMesh;
    if (world.rank() == world.size() - 1) {
        upperBoundMesh = box.leftplane.size() - 1;
    }

    // Car-Parrinello Molecular Dynamics
    cpmd(ion, box, real_bath, fake_bath, bin, fmdremote, cpmdremote);

//   // Post simulation analysis (useful for short runs, but performed otherwise too)
    //   auto_correlation_function();
    if (world.rank() == 0) {
        cout << "MD trust factor R (should be < 0.05) is " << compute_MD_trust_factor_R(cpmdremote.hiteqm) << endl;
        cout << "Program ends \n\n";
    }

    return 0;
}
// End of main

// useful codes

/* random number generation code from gsl
int N = 10;  // number of random numbers to generate

// random nmbr declarations

const gsl_rng_type * T;
gsl_rng * r;

T = gsl_rng_default;
r = gsl_rng_alloc (T);


// prepare for seeding of gsl random nmbr generator:

//   srand((time(0)));                // srand & time are built-in
//   unsigned long int s = random();  // gsl_rng_uniform will eventually
unsigned long int s = 23897897;  // gsl_rng_uniform will eventually
                                 // want a non-negative "long" integer


gsl_rng_env_setup();
gsl_rng_set(r,s); // seed the random number generator;

// generate N random numbers:

for(int i = 0; i < N; ++i) {
  cout << " random number = " << gsl_rng_uniform (r) << endl;;
  }

gsl_rng_free (r);

return 0;
*/

