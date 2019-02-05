// This is a header file for the INTERFACE class.  

#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "utility.h"
#include "vector3d.h"
#include "particle.h"
#include "vertex.h"

class INTERFACE
{
  public:

  VECTOR3D posvec;		// position vector of the inteface, for sphere its center
  double ein; 			// permittivity of inside medium
  double e_left;            // permittivity of left medium
  double e_right;           // permittivity of right medium
  double em_left;			// permittivity at the left interface, mean
  double em_right;			// permittivity at the right interface, mean
  double ed_left;			// permittivity change at the left inteface, difference scaled by 4 pi
  double ed_right;			// permittivity change at the right inteface, difference scaled by 4 pi
  double lx;			// length of the box in x direction
  double ly;			// length of the box in y direction
  double lz; 			// length of the box in z direction
  double lB_in;			// Bjerrum length inside
  double lB_left;           // Bjerrum length left side
  double lB_right;          // Bjerrum length right side
  double inv_kappa_in;		// debye length inside
  double inv_kappa_left;	// debye length left side
  double inv_kappa_right;	// debye length right side
  double mean_sep_in;		// mean separation inside
  double mean_sep_left;		// mean separation left side
  double mean_sep_right;	// mean separation right side
  
  double width;			// discretization width
  vector<VERTEX> leftplane;
  vector<VERTEX> rightplane;

  vector<long double> ndot_grad_greens_innersum_left;
  vector<long double> greens_innersum_left;
  vector<long double> ndot_grad_greens_innersum_right;
  vector<long double> greens_innersum_right;

  void set_up(double, double,  double, int, int, int,double, double, double);
  vector<VERTEX> compute_exact_induced_density_left (vector<PARTICLE>&);
  vector<VERTEX> compute_exact_induced_density_right (vector<PARTICLE>&);
  double exact_total_energy(vector<VERTEX>& , vector<VERTEX>& ,vector<PARTICLE>& );
  vector<VERTEX> compute_exact_induced_density (vector<PARTICLE>&, char);
  void put_counterions(vector<PARTICLE>&, int, double, vector<PARTICLE>&);
  void put_saltions_inside(vector<PARTICLE>&, int, int, double, double, vector<PARTICLE>&);
  void put_one_saltion_in_center(vector<PARTICLE>&, int, int, double, vector<PARTICLE>&);
  void put_saltions_leftside(vector<PARTICLE>&, int, double, double, vector<PARTICLE>&);
  void put_saltions_rightside(vector<PARTICLE>&, int, double, double, vector<PARTICLE>&);
  void discretize(double, double);
  
  INTERFACE(VECTOR3D posvec = VECTOR3D(0,0,0), double ein = 1) : posvec(posvec),  ein(ein)
  {
  }
  
  // total charge inside
  double total_charge_inside(vector<PARTICLE>& ion)
  {
    double charge = 0;
    for (unsigned int i = 0; i < ion.size(); i++)
    {
      double r = 0.5*ion[i].diameter;
      	if (ion[i].posvec.x < lx-r && ion[i].posvec.x > -lx+r) 
	{	
	  if (ion[i].posvec.y < ly-r && ion[i].posvec.y > -ly+r)
	  {	
	    if (ion[i].posvec.z < lz-r && ion[i].posvec.z > -lz+r) 
	      charge += ion[i].q;
	  }
	}	
    }
    return charge;
  }
  double total_induced_charge_left()
  {
    double charge = 0;
    for (unsigned int k = 0; k < leftplane.size(); k++)
      charge += leftplane[k].w * leftplane[k].a;
    return charge;
  }
  double total_induced_charge_right()
  {
    double charge = 0;
    for (unsigned int k = 0; k < rightplane.size(); k++)
      charge += rightplane[k].w * rightplane[k].a;
    return charge;
  }
  // total induced charge
  double total_induced_charge()
  {
    return(total_induced_charge_left() + total_induced_charge_right());
  }

};

#endif

