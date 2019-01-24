// This is fictitious molecular dynamics
// This program is used to estimate the correct w(k)'s on the interface

#include "fmd.h"

void fmd( vector<PARTICLE>& ion, INTERFACE& box, CONTROL& fmdremote, CONTROL& cpmdremote)
{

  // Part I : Initialize
  for (unsigned int k = 0; k <  box.leftplane.size(); k++)
  {
       box.leftplane[k].mu = fmdremote.fakemass * box.leftplane[k].a *  box.leftplane[k].a;						// Assign mass to the fake degree
   // s[k].w = 0.0;								// Initialize fake degree value		(unconstrained)
       box.leftplane[k].vw = 0.0;								// Initialize fake degree velocity	(unconstrained)
  }
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
  {
      box.rightplane[k].mu = fmdremote.fakemass *box.rightplane[k].a * box.rightplane[k].a;						// Assign mass to the fake degree
   // s[k].w = 0.0;								// Initialize fake degree value		(unconstrained)
      box.rightplane[k].vw = 0.0;								// Initialize fake degree velocity	(unconstrained)
  }

/*
  long double sigma = constraint( ion, box);
  for (unsigned int k = 0; k <  box.leftplane.size(); k++)
     box.leftplane[k].w =  box.leftplane[k].w - sigma / ( box.leftplane[k].a *( box.leftplane.size()+box.rightplane.size()) );				// Force to satisfy constraint
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].w = box.rightplane[k].w - sigma / (box.rightplane[k].a *( box.leftplane.size()+box.rightplane.size()) );				// Satisfy constraint

 long double sigmadot = dotconstraint(box);
  for (unsigned int k = 0; k <  box.leftplane.size(); k++)
     box.leftplane[k].vw = box.leftplane[k].vw - sigmadot / ( box.leftplane[k].a * ( box.leftplane.size()+box.rightplane.size()));				// Force to atisfy time derivative of the constraint
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
     box.rightplane[k].vw =box.rightplane[k].vw - sigmadot / (box.rightplane[k].a * ( box.leftplane.size()+box.rightplane.size()));				// Satisfy time derivative of the constraint
*/

  long double sigma_left = constraint_left(ion, box);
  for (unsigned int k = 0; k <  box.leftplane.size(); k++)
     box.leftplane[k].w =  box.leftplane[k].w - sigma_left / ( box.leftplane[k].a *( box.leftplane.size()) );				// Force to satisfy constraint
  long double sigma_right = constraint_right( ion, box);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].w = box.rightplane[k].w - sigma_right / (box.rightplane[k].a *( box.rightplane.size()) );				// Satisfy constraint

 long double sigmadot_left = dotconstraint_left(box);
  for (unsigned int k = 0; k <  box.leftplane.size(); k++)
     box.leftplane[k].vw = box.leftplane[k].vw - sigmadot_left / ( box.leftplane[k].a * ( box.leftplane.size()));				// Force to atisfy time derivative of the constraint
  long double sigmadot_right = dotconstraint_right(box);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
     box.rightplane[k].vw =box.rightplane[k].vw - sigmadot_right / (box.rightplane[k].a * (box.rightplane.size()));				// Satisfy time derivative of the constraint

  compute_force_fake_dof(ion, box, 'n');				// Compute initial force

  long double kinetic_energy = 0.0;
  kinetic_energy = fake_kinetic_energy(box);				// Compute initial kinetic energy
  double potential_energy = energy_functional(ion, box);	// Compute initial potential energy
  cout << "potential energy :  " << potential_energy<<endl;
  fmdremote.annealfreq = 10;
  fmdremote.hiteqm = 100;
  fmdremote.freq =10;

  // Output fmd essentials
  cout << "\n";
  cout << "F M D" << " at " << fmdremote.verify << endl;
  cout << "Mass assigned to the fake degrees " << box.leftplane[0].mu << endl;
  cout << "Total induced charge on the interface " << box.total_induced_charge() << endl;
  cout << "Constraint is (zero if satisfied) " << constraint(ion, box) << endl;
  cout << "Time derivative of the constraint is " << dotconstraint(box) << endl;
  cout << "Initial force on fake degree at vertex 0 at left interface " << box.leftplane[0].fw << endl;
  cout << "Initial force on fake degree at vertex 0 at right interface" << box.rightplane[0].fw << endl;
  cout << "Initial fake kinetic energy " << kinetic_energy << endl;
  cout << "Inital potential energy " << potential_energy << endl;
  cout << "Initial total energy " << kinetic_energy + potential_energy << endl;
  cout << "Time step " << fmdremote.timestep << endl;
  cout << "Number of steps " << fmdremote.steps << endl;


  // create fmd output files
/*  ofstream fmdv, fmde, fmdtic;
  char data[200];
  sprintf(data, "outfiles/fmdv_%.06d.dat", fmdremote.verify);
  fmdv.open(data);
  sprintf(data, "outfiles/fmde_%.06d.dat", fmdremote.verify);
  fmde.open(data);
  sprintf(data, "outfiles/fmdtic_%.06d.dat", fmdremote.verify);
  fmdtic.open(data);
*/
  //cout << "Time step " << fmdremote.timestep << endl;

  double average_total_induced_charge = 0.0;					// initialize useful averages
  double average_total_induced_charge_left = 0.0;					// initialize useful averages
  double average_total_induced_charge_right= 0.0;					// initialize useful averages
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    box.leftplane[k].wmean = 0;
  for (unsigned int l = 0; l < box.rightplane.size(); l++)
    box.rightplane[l].wmean = 0;
  int samples = 0;								// number of samples

  fmdremote.extra_compute = 1000;						// new addition, if one wants to have longer fmd since the startind density is not coming out right
  //cout << "extra_compute " <<  fmdremote.extra_compute << endl;

  // PART II : Propagate
  /*.......................................Fictitious molecular dynamics.......................................*/

  for (int num = 1; num <= fmdremote.steps; num++)
  {
 //   if (num % 50 == 0) cout << "what step is this " << num << endl;
//     cout << num <<endl;
//     fmdv << num << "  " << box.leftplane[0].w << "  " << box.rightplane[0].w << "  " << box.leftplane[0].fw << "  " << box.rightplane[0].fw <<  "  " << box.leftplane[0].vw << "  " << box.rightplane[0].vw <<endl;

    // total induced charge
//     fmdtic << num << "  " << box.total_induced_charge_left() << "  " << box.total_induced_charge_right() << endl;

    // INTEGRATOR
    //! begins
    for (unsigned int k = 0; k < box.leftplane.size(); k++)
      box.leftplane[k].update_velocity(fmdremote.timestep);		// update velocity half time step
    for (unsigned int k = 0; k < box.rightplane.size(); k++)
      box.rightplane[k].update_velocity(fmdremote.timestep);		// update velocity half time step
    for (unsigned int k = 0; k <box.leftplane.size(); k++)
      box.leftplane[k].update_position(fmdremote.timestep);		// update position full time step
    for (unsigned int k = 0; k <box.rightplane.size(); k++)
      box.rightplane[k].update_position(fmdremote.timestep);		// update position full time step
    //SHAKE_left(ion, box, fmdremote);			// shake to ensure constraint
    //SHAKE_right(ion, box, fmdremote);			// shake to ensure constraint
    compute_force_fake_dof( ion, box, 'n');	// calculate forces using the new positions
    for (unsigned int k = 0; k < box.leftplane.size(); k++)
      box.leftplane[k].update_velocity(fmdremote.timestep);		// update velocity half time step
    for (unsigned int k = 0; k < box.rightplane.size(); k++)
      box.rightplane[k].update_velocity(fmdremote.timestep);		// update velocity half time step
    //RATTLE_left(box);						// rattle to ensure time derivative of the constraint
    //RATTLE_right(box);						// rattle to ensure time derivative of the constraint

    //! ends

    // anneal (only if necessary)
    if (fmdremote.anneal == 'y' && num >= fmdremote.hiteqm &&  num % fmdremote.annealfreq == 0)
    {
      for (unsigned int k = 0; k < box.leftplane.size(); k++)
        box.leftplane[k].vw = 0.9 * box.leftplane[k].vw;
      for (unsigned int k = 0; k < box.rightplane.size(); k++)
        box.rightplane[k].vw = 0.9 * box.rightplane[k].vw;
    }

    // sample collection for averages
    if (num > fmdremote.hiteqm && (num % fmdremote.freq) == 0)
    {
      samples = samples + 1;
      for (unsigned int k = 0; k < box.leftplane.size(); k++)
        box.leftplane[k].wmean += box.leftplane[k].w;
      for (unsigned int k = 0; k < box.rightplane.size(); k++)
        box.rightplane[k].wmean += box.rightplane[k].w;
      average_total_induced_charge += box.total_induced_charge();
      average_total_induced_charge_left += box.total_induced_charge_left();
      average_total_induced_charge_right += box.total_induced_charge_right();
    }


/*    // additional computations (turn off for faster simulations)
    if (num % fmdremote.extra_compute == 0)
    {
      kinetic_energy =0.0;
      kinetic_energy = fake_kinetic_energy(box);				// Compute  kinetic energy
      double potential_energy = energy_functional(ion, box);
      double extended_energy = kinetic_energy + potential_energy;
      fmde << num << "  " << extended_energy << "  " << kinetic_energy << "  " <<  potential_energy << endl;
    }
*/
  }
 /*.................................................fmd ends..................................................*/


  // Part III: Compute averages
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    box.leftplane[k].wmean /= samples;				// average induced charge at each vertex
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].wmean /= samples;				// average induced charge at each vertex
  average_total_induced_charge /= samples;		// average total induced charge
  average_total_induced_charge_left /= samples;		// average total induced charge
  average_total_induced_charge_right /= samples;		// average total induced charge

  cout << "Number of samples used to estimate induced charges " << samples << endl;
  cout << "Total induced charge on average " << average_total_induced_charge << endl;
  cout << "Total induced charge leftplane on average " << average_total_induced_charge_left << endl;
  cout << "Total induced charge rightplane on average " << average_total_induced_charge_right << endl;

  return;
}
