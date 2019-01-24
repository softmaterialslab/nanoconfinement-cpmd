// This is thermostat class

#ifndef _THERMOSTAT_H
#define _THERMOSTAT_H

#include "utility.h"

class THERMOSTAT 
{
  public:

  // members
  double Q;			// mass of the thermostat
  double T;			// temperature to be set by thermostat
  int dof;			// degrees of freedom
  double xi;			// thermostat variable
  double eta;			// thermostat variable - useful for computing extended energy
  double pe;			// potential energy of thermostat
  double ke;			// kinetic energy of thermostat
  double hold;			// hold of the bath, which bath holds to a fixed value
  
  // member functions
  
  // make a thermostat
  THERMOSTAT(double mass = 0.0, double temperature = 0.0, int dof = 0, double xi = 0.0, double eta = 0.0, double hold = 0.0) : Q(mass), T(temperature), dof(dof), xi(xi), eta(eta), hold(hold)
  {
  }
  
  // update xi
  void update_xi(double IT, double dt) 
  {
    if (Q == 0)
      return;
    xi = xi + 0.5 * dt * (1.0 / Q) * (IT - dof * kB * T);	// not used. valid only for a solitary thermostat.
    return;
  }

  // update eta
  void update_eta(double dt) 
  {
    if (Q == 0)
      return;
    eta = eta + 0.5 * dt * xi;
    return;
  }

  // calculate potential energy
  void potential_energy() 
  {
    pe = dof * kB * T * eta;				// eta is zero for dummy making pe 0
    return;
  }

  // calculate kinetic energy
  void kinetic_energy() 
  {
    ke = 0.5 * Q * xi * xi;				// Q is zero for dummy making ke 0
    return;
  }
};

#endif
