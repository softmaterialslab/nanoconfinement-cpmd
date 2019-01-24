// This is Car-Parrinello molecular dynamics (CPMD)

#include "particle.h"
#include "vertex.h"
#include "interface.h"
#include "thermostat.h"
#include "control.h"
#include "forces.h"
#include "energies.h"
#include "functions.h"

void cpmd(vector<PARTICLE>& ion, INTERFACE& box, vector<THERMOSTAT>& real_bath, vector<THERMOSTAT>& fake_bath, vector<BIN>& bin, CONTROL& fmdremote,  CONTROL& cpmdremote)
{
  
  // Part I: Initialize
  // initialize fake positions and velocities
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
  {
    box.leftplane[k].mu = cpmdremote.fakemass * box.leftplane[k].a * box.leftplane[k].a;
    box.leftplane[k].w = box.leftplane[k].wmean;
  }
  initialize_fake_velocities(box.leftplane, fake_bath);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
  {
    box.rightplane[k].mu = cpmdremote.fakemass * box.rightplane[k].a * box.rightplane[k].a;
    box.rightplane[k].w = box.rightplane[k].wmean;
  }
  initialize_fake_velocities(box.rightplane, fake_bath);
  // satisfy constraints
  long double sigma_left = constraint_left(ion, box);
  for (unsigned int k = 0; k <  box.leftplane.size(); k++)
     box.leftplane[k].w =  box.leftplane[k].w - sigma_left / ( box.leftplane[k].a *( box.leftplane.size()) );		
  long double sigma_right = constraint_right( ion, box);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].w = box.rightplane[k].w - sigma_right / (box.rightplane[k].a *( box.rightplane.size()) );	
  long double sigmadot_left = dotconstraint_left(box);
  for (unsigned int k = 0; k <  box.leftplane.size(); k++)
     box.leftplane[k].vw = box.leftplane[k].vw - sigmadot_left / ( box.leftplane[k].a * ( box.leftplane.size()));		
  long double sigmadot_right = dotconstraint_right(box);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
     box.rightplane[k].vw =box.rightplane[k].vw - sigmadot_right / (box.rightplane[k].a * (box.rightplane.size()));
  long double fake_ke = fake_kinetic_energy(box);
  // initialize real positions and velocities
  initialize_particle_velocities(ion, real_bath);		// particle velocities initialized
//  compute_force_fake_dof( ion, box, 'n');			// calculate forces using the new positions
  compute_force_real_dof(ion, box, 'y');			// force on particles and fake degrees initialized
  long double particle_ke = particle_kinetic_energy(ion);	// compute initial kinetic energy
  long double potential_energy; 
  potential_energy= energy_functional(ion, box);	// compute initial potential energy
  
  // Output cpmd essentials
  cout << "\n";
  cout << "C P M D" << " on " << endl;
  cout << "Mass assigned to the fake degrees " << box.leftplane[0].mu << setw(10) << box.rightplane[0].mu << endl;
  cout << "Total induced charge on leftplane " << box.total_induced_charge_left() << endl;
  cout << "Total induced charge on rightplane " << box.total_induced_charge_right() << endl;
  cout << "Constraint on leftplane " << constraint_left(ion, box) << endl;
  cout << "Constraint on rightplane " << constraint_right(ion, box) << endl;
  cout << "Dot Constraint on leftplane " << dotconstraint_left(box) << endl;
  cout << "Dot Constraint on rightplane " << dotconstraint_right(box) << endl;
  cout << "Initial force on fake degree at vertex 0 on left plane " << box.leftplane[0].fw << endl;
  cout << "Initial force on fake degree at vertex 0 on right plane " << box.rightplane[0].fw << endl;
  cout << "Initial force on ion 0 " << ion[0].forvec << endl;
  cout << "Initial force on ion 1 " << ion[1].forvec << endl;
  cout << "Initial force on fake degree at vertex 0 on right plane " << box.rightplane[0].fw << endl;
  cout << "Initial ion kinetic energy " << particle_ke << endl; 
  cout << "Initial fake kinetic energy " << fake_ke << endl; 
  cout << "Inital potential energy " << potential_energy << endl;
  cout << "Initial system energy " << particle_ke + fake_ke + potential_energy << endl;
  cout << "Real (ion) chain length (L+1) implementation " << real_bath.size() << endl;
  cout << "Real thermostat temperature " << real_bath[0].T << endl;
  cout << "Real thermostat mass " << real_bath[0].Q << endl;
  cout << "Fake chain length (L+1) implementation " << fake_bath.size() << endl;
  cout << "Fake thermostat temperature " << fake_bath[0].T << endl;
  cout << "Fake thermostat mass " << fake_bath[0].Q << endl;
  cout << "Number of bins used for computing density profiles " << bin.size() << endl;
  cout << "Time step " << cpmdremote.timestep << endl;
  cout << "Number of steps " << cpmdremote.steps << endl;
  cout << "Production begins at " << cpmdremote.hiteqm << endl;
  cout << "Sampling frequency " << cpmdremote.freq << endl;
  cout << "Extra computation every " << cpmdremote.extra_compute <<  " steps" << endl;
  cout << "Verification every " << cpmdremote.verify << " steps" << endl;
  cout << "Write density profile every " << cpmdremote.writedensity << endl;
  
//   for (unsigned int i = 0; i < ion.size(); i++)
//     cout << "ion   " << i << " start position  " << ion[i].posvec << endl;
  
  // for movie
  int moviestart = 1;					// starting point of the movie
  
  // for energy
  double energy_samples = 0;
  
  // for verification with fmd (exact results)
  double verification_samples = 0;
  double average_functional_deviation = 0;

  // for density profile
  vector<double> mean_positiveion_density;				// average density profile
  vector<double> mean_negativeion_density;				// average density profile
  vector<double> mean_sq_positiveion_density;			// average of square of density
  vector<double> mean_sq_negativeion_density;			// average of square of density
  for (unsigned int b = 0; b < bin.size(); b++)
  {
    mean_positiveion_density.push_back(0.0);
    mean_negativeion_density.push_back(0.0);
    mean_sq_positiveion_density.push_back(0.0);
    mean_sq_negativeion_density.push_back(0.0);
  }
  double density_profile_samples = 0;			// number of samples used to estimate density profile
  
  long double expfac_real, expfac_fake;				// exponential factors useful in velocity Verlet routine
  
  /*
  // print out initial position, velocity, force
  cout << "x y z  " <<  ion[0].posvec << endl;
  cout << "vx vy vz  " <<  ion[0].velvec << endl;
  cout << "fx fy fz  " <<  ion[0].forvec << endl;
  */
  
  // TEMPORARY!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING 
//   ofstream list_tic ("outfiles/total_induced_charge.dat", ios::app);  
  
  // Part II : Propagate
  for (int num = 1; num <= cpmdremote.steps; num++)
  {
    //list_tic << num << setw(15) << box.total_induced_charge() << endl; // WARNING WARNING
    
    //if (num % 100 == 0) cout << "what step is this " << num << endl;
    // INTEGRATOR
    //! begins
    // reverse update of Nose-Hoover chain
    // update xi ionic system
    for (int j = real_bath.size() - 1; j > -1; j--)
      update_chain_xi(j, real_bath, cpmdremote.timestep, particle_ke);			
    // update xi fake system
    for (int j = fake_bath.size() - 1; j > -1; j--)
      update_chain_xi(j, fake_bath, cpmdremote.timestep, fake_ke);	
    // update eta ionic system
    for (unsigned int j = 0; j < real_bath.size(); j++)
      real_bath[j].update_eta(cpmdremote.timestep);
    // update eta fake system
    for (unsigned int j = 0; j < fake_bath.size(); j++)
      fake_bath[j].update_eta(cpmdremote.timestep);	
    
    expfac_real = exp(-0.5 * cpmdremote.timestep * real_bath[0].xi);
    expfac_fake = exp(-0.5 * cpmdremote.timestep * fake_bath[0].xi);
    
    // core loop: velocity --> position --> force --> velocity
    for (unsigned int i = 0; i < ion.size(); i++)
      ion[i].new_update_velocity(cpmdremote.timestep, real_bath[0], expfac_real);
    for (unsigned int k = 0; k < box.leftplane.size(); k++)
      box.leftplane[k].new_update_velocity(cpmdremote.timestep, fake_bath[0], expfac_fake);	
    for (unsigned int k = 0; k < box.rightplane.size(); k++)
      box.rightplane[k].new_update_velocity(cpmdremote.timestep, fake_bath[0], expfac_fake);	
    
    for (unsigned int i = 0; i < ion.size(); i++)
      ion[i].update_position(cpmdremote.timestep);
    for (unsigned int k = 0; k < box.leftplane.size(); k++)	
      box.leftplane[k].update_position(cpmdremote.timestep);
    for (unsigned int k = 0; k < box.rightplane.size(); k++)	
      box.rightplane[k].update_position(cpmdremote.timestep);
    
    SHAKE_left(ion, box, cpmdremote);			// shake to ensure constraint
    SHAKE_right(ion, box, cpmdremote);			// shake to ensure constraint
    
  //  compute_force_fake_dof(ion, box, 'n');	// calculate forces using the new positions
    compute_force_real_dof(ion, box, 'y');

    for (unsigned int k = 0; k < box.leftplane.size(); k++)
      box.leftplane[k].new_update_velocity(cpmdremote.timestep, fake_bath[0], expfac_fake);
    for (unsigned int k = 0; k < box.rightplane.size(); k++)
      box.rightplane[k].new_update_velocity(cpmdremote.timestep, fake_bath[0], expfac_fake);
    for (unsigned int i = 0; i < ion.size(); i++)
      ion[i].new_update_velocity(cpmdremote.timestep, real_bath[0], expfac_real);
    
    RATTLE_left(box);						// rattle to ensure time derivative of the constraint
    RATTLE_right(box);						// rattle to ensure time derivative of the constraint  
    
    // kinetic energies needed to set canonical ensemble
    fake_ke = fake_kinetic_energy(box);
    particle_ke = particle_kinetic_energy(ion);
    
    // forward update of Nose-Hoover chain
    for (unsigned int j = 0; j < fake_bath.size(); j++)
      fake_bath[j].update_eta(cpmdremote.timestep);	
    for (unsigned int j = 0; j < real_bath.size(); j++)
      real_bath[j].update_eta(cpmdremote.timestep);
    for (unsigned int j = 0; j < fake_bath.size(); j++)
      update_chain_xi(j, fake_bath, cpmdremote.timestep, fake_ke);
    for (unsigned int j = 0; j < real_bath.size(); j++)
      update_chain_xi(j, real_bath, cpmdremote.timestep, particle_ke);
    //! ends
    
    // extra computations
    if (num == 1 || num % cpmdremote.extra_compute == 0) 
    {
      energy_samples++;
      compute_n_write_useful_data(num, ion, real_bath, fake_bath, box);
      write_basic_files(num, ion, real_bath, box);
    }
    // verify with F M D
    if (num % cpmdremote.verify == 0)
    {
      double functional_deviation = verify_with_FMD(num, ion, box, fmdremote, cpmdremote);
      verification_samples++;
      average_functional_deviation += functional_deviation;
    }
    // make a movie
    if (num >= moviestart && num % cpmdremote.moviefreq == 0) 
      make_movie(num, ion, box); 
//     
//     // compute density profile
     if (num >= cpmdremote.hiteqm && (num % cpmdremote.freq == 0))
     {
       density_profile_samples++;
       compute_density_profile(num, density_profile_samples, mean_positiveion_density, mean_sq_positiveion_density, mean_negativeion_density, mean_sq_negativeion_density, ion, box, bin, cpmdremote);
     }
  }

   // Part III : Analysis
   // 1. density profile
   vector<double> positiveion_density_profile;
   vector<double> negativeion_density_profile;
   for (unsigned int b = 0; b < mean_positiveion_density.size(); b++)
     positiveion_density_profile.push_back(mean_positiveion_density.at(b) / density_profile_samples);
   for (unsigned int b = 0; b < mean_negativeion_density.size(); b++)
     negativeion_density_profile.push_back(mean_negativeion_density.at(b) / density_profile_samples);
   
   // 2. error bars
   vector<double> p_error_bar;
   vector<double> n_error_bar;
   for (unsigned int b = 0; b < positiveion_density_profile.size(); b++)
     p_error_bar.push_back( sqrt(1.0/density_profile_samples) * sqrt( mean_sq_positiveion_density.at(b)/density_profile_samples - positiveion_density_profile.at(b)*positiveion_density_profile.at(b) ) );
   for (unsigned int b = 0; b < negativeion_density_profile.size(); b++)
     n_error_bar.push_back( sqrt(1.0/density_profile_samples) * sqrt( mean_sq_negativeion_density.at(b)/density_profile_samples - negativeion_density_profile.at(b)*negativeion_density_profile.at(b) ) );
   
   // 3. write results
   ofstream list_p_profile ("outfiles/p_denisty_profile.dat", ios::out);
   ofstream list_n_profile ("outfiles/n_denisty_profile.dat", ios::out);
   for (unsigned int b = 0; b < positiveion_density_profile.size(); b++)
     list_p_profile << (-0.5*box.lz + b * bin[b].width) * unitlength << setw(15) << positiveion_density_profile.at(b) << setw(15) << p_error_bar.at(b) << endl; // change in the z coordinate, counted from leftwall
   ofstream list_profile ("outfiles/n_denisty_profile.dat", ios::out);
   for (unsigned int b = 0; b < negativeion_density_profile.size(); b++)
     list_n_profile << (-0.5*box.lz + b * bin[b].width) * unitlength << setw(15) << negativeion_density_profile.at(b) << setw(15) << n_error_bar.at(b) << endl; // change in the z coordinate, counted from leftwall
   
   ofstream final_configuration ("outfiles/final_configuration.dat");
   for (unsigned int i = 0; i < ion.size(); i++)
     final_configuration << ion[i].posvec << endl;
   
 //   ofstream final_induced_density ("outfiles/final_induced_density.dat");
 //   for (unsigned int k = 0; k < s.size(); k++)
 //     final_induced_density << k+1 << setw(15) << s[k].theta << setw(15) << s[k].phi << setw(15) << s[k].w << endl;
   
   cout << "Number of samples used to compute energy" << setw(10) << energy_samples << endl;
   cout << "Number of samples used to get density profile" << setw(10) << density_profile_samples << endl;
  
//   cout << "Number of samples used to verify on the fly results" << setw(10) << verification_samples << endl;
//   cout << "Average deviation of the functional from the BO surface" << setw(15) << average_functional_deviation / verification_samples << endl;

  return ;
}
// End of CPMD routine

// old routine
/*
    //! begins
    // reverse update of Nose-Hoover chain
    for (int j = real_bath.size() - 1; j > -1; j--)
      update_chain_xi(j, real_bath, cpmdremote.timestep, particle_ke);			// update xi for real baths in reverse order
    for (int j = fake_bath.size() - 1; j > -1; j--)
      update_chain_xi(j, fake_bath, cpmdremote.timestep, fake_ke);			// update xi for fake baths in reverse order
    for (unsigned int j = 0; j < real_bath.size(); j++)
      real_bath[j].update_eta(cpmdremote.timestep);					// update eta for real baths
    for (unsigned int j = 0; j < fake_bath.size(); j++)
      fake_bath[j].update_eta(cpmdremote.timestep);					// update eta for fake baths
    // pre-compute expfac_real and expfac_fake to be used in velocity Verlet
    expfac_real = exp(-0.5 * cpmdremote.timestep * real_bath[0].xi);
    expfac_fake = exp(-0.5 * cpmdremote.timestep * fake_bath[0].xi);
    // Modified velocity Verlet (with thermostat effects)
    for (unsigned int i = 0; i < ion.size(); i++)
      ion[i].new_update_velocity(cpmdremote.timestep, real_bath[0], expfac_real);	// update particle velocity half time step
    for (unsigned int k = 0; k < s.size(); k++)
      s[k].new_update_velocity(cpmdremote.timestep, fake_bath[0], expfac_fake);		// update fake velocity half time step
    for (unsigned int i = 0; i < ion.size(); i++)
      ion[i].update_position(cpmdremote.timestep);					// update particle position full time step
    for (unsigned int k = 0; k < s.size(); k++)						// update fake position full time step
      s[k].update_position(cpmdremote.timestep);
    SHAKE(s, ion, dsphere, cpmdremote);							// shake to ensure constraint is satisfied
    for_cpmd_calculate_force(s, ion, dsphere, colloid, 'y');					// calculate forces on ion and fake degree
    for (unsigned int k = 0; k < s.size(); k++)
      s[k].new_update_velocity(cpmdremote.timestep, fake_bath[0], expfac_fake);		// update fake velocity half time step  
    for (unsigned int i = 0; i < ion.size(); i++)
      ion[i].new_update_velocity(cpmdremote.timestep, real_bath[0], expfac_real);	// update particle velocity half time step
    RATTLE(s);										// rattle to ensure time derivative of the constraint is satisfied
    // kinetic energies needed to set canonical ensemble
    fake_ke = fake_kinetic_energy(s);
    particle_ke = particle_kinetic_energy(ion);
    // Forward Nose-Hoover chain
    for (unsigned int j = 0; j < fake_bath.size(); j++)
      fake_bath[j].update_eta(cpmdremote.timestep);					// update eta for fake baths
    for (unsigned int j = 0; j < real_bath.size(); j++)
      real_bath[j].update_eta(cpmdremote.timestep);					// update eta for real baths
    for (unsigned int j = 0; j < fake_bath.size(); j++)
      update_chain_xi(j, fake_bath, cpmdremote.timestep, fake_ke);			// update xi for fake baths in forward order
    for (unsigned int j = 0; j < real_bath.size(); j++)
      update_chain_xi(j, real_bath, cpmdremote.timestep, particle_ke);			// update xi for real baths in forward order
    //! ends
*/

