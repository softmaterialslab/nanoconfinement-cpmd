// This is a header file containing functions

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "utility.h"
#include "interface.h"
#include "particle.h"
#include "vertex.h"
#include "control.h"
#include "BIN.h"
#include "thermostat.h"
#include "energies.h"


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// general functions
// -----------------

// overloaded << to print 3d vectors
ostream& operator<<(ostream&, VECTOR3D);

// make bins
void make_bins(vector<BIN>&, INTERFACE&, double);

// bin ions
void bin_ions(vector<PARTICLE>&, INTERFACE&, vector<double>&, vector<BIN>&);

// initialize particle velocities
void initialize_particle_velocities(vector<PARTICLE>&, vector<THERMOSTAT>&);

// initialize fake velocities
void initialize_fake_velocities(vector<VERTEX>&, vector<THERMOSTAT>&);

// compute and write useful data in cpmd
void compute_n_write_useful_data(int, vector<PARTICLE>&, vector<THERMOSTAT>&, vector<THERMOSTAT>&, INTERFACE&);

// verify with F M D
double verify_with_FMD(int,  vector<PARTICLE>&, INTERFACE&, CONTROL&, CONTROL&);

// make movie
void make_movie(int num, vector<PARTICLE>& ion, INTERFACE& box);

// compute density profile
void compute_density_profile(int, double, vector<double>&, vector<double>&, vector<PARTICLE>&, INTERFACE&, vector<BIN>&, CONTROL&);

// post analysis : compute R
double compute_MD_trust_factor_R(int);

// post analysis : auto correlation function
void auto_correlation_function();

double exact_density_integration(double, double, double);

// functions useful in computing forces and energies
// -------------------------------------------------

// Green's function general
inline long double Greens(VECTOR3D& d_vec)
{
   return 1.0 / d_vec.GetMagnitude();
}
// Green's function G for induced charge within one  interface
inline long double Greens_intra(VECTOR3D& d_vec, vector<VERTEX>& s, unsigned int k , unsigned int l)
{
  if (l != k) return 1.0 / d_vec.GetMagnitude();
  else  return (2 * sqrt(pi) * sqrt(s[k].a) / s[k].a);
}
// Green's function G for induced charge between two different interfaces
inline long double Greens_inter(VECTOR3D& d_vec)
{
   return 1.0 / d_vec.GetMagnitude();
}
// computes gradient of green's fn
inline VECTOR3D Grad_Greens(VECTOR3D& d_vec)
{
  long double r = d_vec.GetMagnitude();
  if (r == 0) return VECTOR3D(0,0,0);
  long double r3 = r * r * r;
  return d_vec ^ ( (-1.0) / r3 );
}
// Gradient of dotGreens functions between ions and interfaces, used on getting forces on ions
inline VECTOR3D Grad_ndot_Greens(VECTOR3D& normal, VECTOR3D& d_vec)
{
    long double factor = d_vec * normal;
    long double r1 = d_vec.GetMagnitude();
    long double r3 = r1*r1*r1;
    long double r5 = r3*r1*r1;
    long double temp = -3 * factor/r5;

    return ((d_vec ^ temp) + (normal ^ (1/r3)));
}

// Green's function G for induced charge within one  interface
inline long double G(vector<VERTEX>& s, unsigned int k , unsigned int l)
{
  if (l != k) return 1.0 / (s[k].posvec - s[l].posvec).GetMagnitude();
  else  return (2 * sqrt(pi) * sqrt(s[k].a) / s[k].a);
}
// Green's function G for induced charge between two different interfaces
inline long double G_inter(vector<VERTEX>& s1, vector<VERTEX>& s2, unsigned int k , unsigned int l)
{
   return 1.0 / (s1[k].posvec - s2[l].posvec).GetMagnitude();
}
// computes gradient of green's fn
inline VECTOR3D Grad(VECTOR3D& vec1, VECTOR3D& vec2)
{
  long double r = (vec1 - vec2).GetMagnitude();
  long double r3 = r * r * r;
  return (vec1 - vec2) ^ ( (-1.0) / r3 );
}

// Gradient of Green's function  within one  interface
inline VECTOR3D Grad_G(vector<VERTEX>& s, unsigned int k, unsigned int l)
{
  if (l != k) return  Grad(s[k].posvec, s[l].posvec);
  else return VECTOR3D(0,0,0);
}

// Gradient of Green's function between two different interfaces
inline VECTOR3D Grad_G_inter(vector<VERTEX>& s1, vector<VERTEX>& s2, unsigned int k, unsigned int l)
{
  return  Grad(s1[k].posvec, s2[l].posvec);
}

// Gradient of dotGreens functions between ions and interfaces, used on getting forces on ions
inline VECTOR3D Grad_dot_G(VECTOR3D& s, VECTOR3D& normal, VECTOR3D & r)
{
    long double factor = (r - s) * normal;
    long double r1 = (r -s).GetMagnitude();
    long double r3 = r1*r1*r1;
    long double r5 = r3*r1*r1;
    long double temp = -3 * factor/r5;

    return (((r - s) ^ temp) + (normal ^ (1/r3)));
}

// Mimimium image convention for periodic boundary condition
inline VECTOR3D Mini_Image(VECTOR3D& d, double L)
{
  
  d.x = d.x - round(d.x/ L) * L;
  d.y = d.y - round(d.y/ L) * L;
  
  return d;
}

// functions useful in implementing constraint
// -------------------------------------------

// constraint equation
inline long double constraint( vector<PARTICLE>& ion, INTERFACE& box)
{
  // for one closed interface
  return (box.total_induced_charge() - box.total_charge_inside(ion)*(1/box.e_left - 1/box.ein));
}

inline long double constraint_left( vector<PARTICLE>& ion, INTERFACE& box)
{
  // for one infinite planar interface
  if(ion.size() ==1 && box.e_right == box.ein) return (box.total_induced_charge_left() - (box.total_charge_inside(ion) / box.ein) * (box.ein - box.e_left) / (box.ein + box.e_left));
  else return (box.total_induced_charge_left() );
}

inline long double constraint_right( vector<PARTICLE>& ion, INTERFACE& box)
{
  // for one infinite planar interface
  if(ion.size() == 1 && box.e_left == box.ein)return (box.total_induced_charge_right() - (box.total_charge_inside(ion) / box.ein) * (box.ein - box.e_right) / (box.ein + box.e_right));
  else return (box.total_induced_charge_right());
}

// SHAKE to ensure constraint is true
inline void SHAKE(vector<PARTICLE>& ion, INTERFACE& box, CONTROL& simremote)	// remote of the considered simulation
{
  long double sigma = constraint(ion, box);
  int N = box.leftplane.size() + box.rightplane.size();
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    box.leftplane[k].vw = box.leftplane[k].vw - (1.0/simremote.timestep) * sigma / (box.leftplane[k].a * N);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].vw = box.rightplane[k].vw - (1.0/simremote.timestep) * sigma / (box.rightplane[k].a * N);
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
   box.leftplane[k].w = box.leftplane[k].w - sigma / (box.leftplane[k].a * N);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
   box.rightplane[k].w = box.rightplane[k].w - sigma / (box.rightplane[k].a * N);

  return;
}

inline void SHAKE_left(vector<PARTICLE>& ion, INTERFACE& box, CONTROL& simremote)	// remote of the considered simulation
{
  long double sigma = constraint_left(ion, box);
  int N = box.leftplane.size();
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    box.leftplane[k].vw = box.leftplane[k].vw - (1.0/simremote.timestep) * sigma / (box.leftplane[k].a * N);
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
   box.leftplane[k].w = box.leftplane[k].w - sigma / (box.leftplane[k].a * N);

  return;
}
inline void SHAKE_right(vector<PARTICLE>& ion, INTERFACE& box, CONTROL& simremote)	// remote of the considered simulation
{
  long double sigma = constraint_right(ion, box);
  int N =  box.rightplane.size();
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].vw = box.rightplane[k].vw - (1.0/simremote.timestep) * sigma / (box.rightplane[k].a * N);
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
   box.rightplane[k].w = box.rightplane[k].w - sigma / (box.rightplane[k].a * N);

  return;
}

// dot constraint equation
inline long double dotconstraint(INTERFACE& box)
{
  long double sigmadot = 0;
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    sigmadot += box.leftplane[k].vw * box.leftplane[k].a;
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    sigmadot += box.rightplane[k].vw * box.rightplane[k].a;
  return sigmadot;
}

inline long double dotconstraint_left(INTERFACE& box)
{
  long double sigmadot = 0;
  for (unsigned int k = 0; k < box.leftplane.size(); k++)
    sigmadot += box.leftplane[k].vw * box.leftplane[k].a;
  return sigmadot;
}
inline long double dotconstraint_right(INTERFACE& box)
{
  long double sigmadot = 0;
  for (unsigned int k = 0; k < box.rightplane.size(); k++)
    sigmadot += box.rightplane[k].vw * box.rightplane[k].a;
  return sigmadot;
}
// RATTLE to ensure time derivative of the constraint is true
inline void RATTLE(INTERFACE& box)
{
 long double sigmadot = dotconstraint(box);
   int N = box.leftplane.size() + box.rightplane.size();
   for (unsigned int k = 0; k < box.leftplane.size(); k++)
    box.leftplane[k].vw = box.leftplane[k].vw - sigmadot / (box.leftplane[k].a * N);
   for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].vw = box.rightplane[k].vw - sigmadot / (box.rightplane[k].a * N);

  return;
}
inline void RATTLE_left(INTERFACE& box)
{
 long double sigmadot = dotconstraint_left(box);
   int N = box.leftplane.size();
   for (unsigned int k = 0; k < box.leftplane.size(); k++)
    box.leftplane[k].vw = box.leftplane[k].vw - sigmadot / (box.leftplane[k].a * N);

  return;
}
inline void RATTLE_right(INTERFACE& box)
{
 long double sigmadot = dotconstraint_right(box);
   int N = box.rightplane.size();
   for (unsigned int k = 0; k < box.rightplane.size(); k++)
    box.rightplane[k].vw = box.rightplane[k].vw - sigmadot / (box.rightplane[k].a * N);

  return;
}

// functions useful in Nose-Hoover chain implementation
// -------------------------------------------

// update bath xi value
inline void update_chain_xi(unsigned int j, vector<THERMOSTAT>& bath, double dt, long double ke)
{
  if (bath[j].Q == 0)
    return;
  if (j != 0)
    bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j+1].xi) + 0.5 * dt * (1.0 / bath[j].Q) * (bath[j-1].Q * bath[j-1].xi * bath[j-1].xi - bath[j].dof * kB * bath[j].T) * exp(-0.25 * dt * bath[j+1].xi);
//     bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j+1].xi) + 0.5 * dt * (1.0 / bath[j].Q) * (bath[j-1].Q * bath[j-1].xi * bath[j-1].xi - bath[j].dof * kB * bath[j].T) * (exp(-0.5 * dt * bath[j+1].xi) - 1) / (-0.5 * dt * bath[j+1].xi);
  else
    bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j+1].xi) + 0.5 * dt * (1.0 / bath[j].Q) * (2*ke - bath[j].dof * kB * bath[j].T) * exp(-0.25 * dt * bath[j+1].xi);
//     bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j+1].xi) + 0.5 * dt * (1.0 / bath[j].Q) * (2*ke - bath[j].dof * kB * bath[j].T) * (exp(-0.5 * dt * bath[j+1].xi) - 1) / (-0.5 * dt * bath[j+1].xi);
  return;
}

// functions that write data files
// -------------------------------

inline void write_basic_files(int cpmdstep, vector<PARTICLE>& ion, vector<THERMOSTAT>& real_bath, INTERFACE& box)
{
  ofstream list_position ("outfiles/ion_position.dat", ios::app);
  ofstream list_velocity ("outfiles/ion_velocity.dat", ios::app);
  ofstream list_force ("outfiles/ion_force.dat", ios::app);
  ofstream list_bath ("outfiles/bath_values.dat", ios::app);
  
  list_position << cpmdstep << setw(15) << ion[0].posvec << setw(15) << ion[1].posvec << endl; 
  list_velocity << cpmdstep << setw(15) << ion[0].velvec << setw(15) << ion[1].velvec << endl; 
  list_force << cpmdstep << setw(15) << ion[0].forvec << setw(15) << ion[1].forvec << endl; 
  list_bath << cpmdstep << setw(15) << real_bath[0].xi << setw(15) << real_bath[0].eta << endl;
  return;
}

// bin ions to get density profile
inline void bin_ions(vector<PARTICLE>& ion, INTERFACE& box, vector<double>& density, vector<BIN>& bin)
{
  double r;
  int bin_number;
  for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
    bin[bin_num].n = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    r = 0.5*box.lz + ion[i].posvec.z;	// r should be positive, counting from half lz and that is the adjustment here; bin 0 corresponds to left wall
    bin_number = int(r/bin[0].width);
    bin[bin_number].n = bin[bin_number].n + 1;
  }
  for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
    density.push_back(bin[bin_num].n / bin[bin_num].volume);			// push_back is the culprit, array goes out of bound
									// volume now measured in inverse Molars, such that density is in M
  return;
}

// compute density profile
inline void compute_density_profile(int cpmdstep, double density_profile_samples,
				    vector<double>& mean_positiveion_density, 
				    vector<double>& mean_sq_positiveion_density,
				    vector<double>& mean_negativeion_density,
				    vector<double>& mean_sq_negativeion_density,
				    vector<PARTICLE>& ion, INTERFACE& box, 
				    vector<BIN>& bin, CONTROL& cpmdremote)
{
  vector<double> sample_positiveion_density;
  vector<double> sample_negativeion_density;
  
  vector<PARTICLE> positiveion;
  vector<PARTICLE> negativeion;
  
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    if (ion[i].valency > 0)
      positiveion.push_back(ion.at(i));
    else if (ion[i].valency < 0)
      negativeion.push_back(ion.at(i));
  }
      
  bin_ions(positiveion, box, sample_positiveion_density, bin);
  bin_ions(negativeion, box, sample_negativeion_density, bin);
  
  for (unsigned int b = 0; b < mean_positiveion_density.size(); b++)
    mean_positiveion_density.at(b) = mean_positiveion_density.at(b) + sample_positiveion_density.at(b);
  for (unsigned int b = 0; b < mean_negativeion_density.size(); b++)
    mean_negativeion_density.at(b) = mean_negativeion_density.at(b) + sample_negativeion_density.at(b);
  for (unsigned int b = 0; b < sample_positiveion_density.size(); b++)
    mean_sq_positiveion_density.at(b) = mean_sq_positiveion_density.at(b) + sample_positiveion_density.at(b)*sample_positiveion_density.at(b);
  for (unsigned int b = 0; b < sample_negativeion_density.size(); b++)
    mean_sq_negativeion_density.at(b) = mean_sq_negativeion_density.at(b) + sample_negativeion_density.at(b)*sample_negativeion_density.at(b);
  
  // write files
  if (cpmdstep % cpmdremote.writedensity == 0)
  {
    char datap[200], datan[200];
    sprintf(datap, "datafiles/_z+_den_%.06d.dat", cpmdstep);
    sprintf(datan, "datafiles/_z-_den_%.06d.dat", cpmdstep);
    ofstream outdenp, outdenn;
    outdenp.open(datap);
    outdenn.open(datan);
    for (unsigned int b = 0; b < mean_positiveion_density.size(); b++)
      outdenp << (-box.lz/2+b*bin[b].width) * unitlength << setw(15) << mean_positiveion_density.at(b)/density_profile_samples << endl;
    for (unsigned int b = 0; b < mean_negativeion_density.size(); b++)
      outdenn << (-box.lz/2+b*bin[b].width) * unitlength << setw(15) << mean_negativeion_density.at(b)/density_profile_samples << endl;
    outdenp.close();
    outdenn.close();
  } 
  return;
}

// display progress bar (code from the internet)
void progressBar(double);

#endif
