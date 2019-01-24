// This is a header file containing the energy routines for induced charges and particles

#ifndef _ENERGIES_H
#define _ENERGIES_H

#include "interface.h"
#include "particle.h"
#include "vertex.h"
#include "thermostat.h"
#include "functions.h"

inline long double fake_kinetic_energy(INTERFACE& box)
{
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    box.leftplane[k].kinetic_energy();
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].kinetic_energy();
  long double kinetic_energy = 0.0;
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    kinetic_energy += box.leftplane[k].ke;  
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    kinetic_energy += box.rightplane[k].ke;  
  return kinetic_energy;
}

inline long double particle_kinetic_energy(vector<PARTICLE>& ion)
{
  for (unsigned int i = 0; i < ion.size(); i++)
    ion[i].kinetic_energy();
  long double kinetic_energy = 0.0;
  for (unsigned int i = 0; i < ion.size(); i++)
    kinetic_energy += ion[i].ke;
  return kinetic_energy;
}

inline double bath_kinetic_energy(vector<THERMOSTAT>& bath)
{
  for (unsigned int j = 0; j < bath.size(); j++)
    bath[j].kinetic_energy();
  double kinetic_energy = 0.0;
  for (unsigned int j = 0; j < bath.size(); j++)
    kinetic_energy += bath[j].ke;
  return kinetic_energy;
}

inline double bath_potential_energy(vector<THERMOSTAT>& bath)
{
  for (unsigned int j = 0; j < bath.size(); j++)
    bath[j].potential_energy();
  double potential_energy = 0.0;
  for (unsigned int j = 0; j < bath.size(); j++)
    potential_energy += bath[j].pe;
  return potential_energy;
}

long double energy_functional(vector<PARTICLE>&, INTERFACE&);

#endif

