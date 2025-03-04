// This is particle class

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"

class PARTICLE 
{
  public:

  // members
  int id;		// id of the particle
  double diameter;	// diameter of the particle
  int valency;		// valency of the ion
  double q;		// charge of the particle
  double m; 		// mass of the particle
  double epsilon;	// dielectric constant of the medium
  VECTOR3D posvec;	// position vector of the particle
  VECTOR3D velvec;	// velocity vector of the particle
  VECTOR3D forvec;	// force vector on the particle
  double pe;		// potential energy
  long double ke;	// kinetic energy
  double energy;	// energy
  double lx;			// box size
  double ly;
  double lz;
  
  // member functions
  
  // make a particle
  PARTICLE(int id = 0, double diameter = 0, int valency = 0, double charge = 0, double mass = 0, double diconst = 0, 
	   VECTOR3D position = VECTOR3D(0,0,0), double lx = 0,double ly = 0, double lz = 0 ): 
	   id(id), diameter(diameter), valency(valency), q(charge), m(mass), epsilon(diconst), posvec(position),lx(lx),ly(ly),lz(lz)
  {
  }
  
  // update position of the particle
  void update_position(double dt)			
  {
    posvec = ( posvec + (velvec ^ dt) );
    if(posvec.x>lx/2)posvec.x -=lx;
    if(posvec.x<-lx/2)posvec.x +=lx;
    if(posvec.y>ly/2)posvec.y -=ly;
    if(posvec.y<-ly/2)posvec.y +=ly;

    return;
  }
  
  // update velocity of the particle
  void update_velocity(double dt)	
  {
    velvec = ( velvec + ( forvec ^ ( 0.5 * (dt/m) ) ) );
    return;
  }
  
  void new_update_velocity(double dt, THERMOSTAT main_bath, long double expfac)
  {
    velvec = ( ( velvec ^ (expfac)  ) + ( forvec ^ (0.5 * (dt/m) * sqrt(expfac)) ) );
    return;
  }
  
  // calculate kinetic energy of a particle
  void kinetic_energy()				
  {
    ke = 0.5 * m * velvec.GetMagnitude() * velvec.GetMagnitude();
    return;
  }
};

#endif
