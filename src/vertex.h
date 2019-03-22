// This is vertex class

#ifndef _VERTEX_H
#define _VERTEX_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"

class VERTEX 
{
  public:

  // members
  VECTOR3D posvec;			// position vector of the vertex
  long double a;			// area of the vertex
  VECTOR3D normalvec;     		// normal vector on the surface pointing interior to exterior
  double r, theta, phi;			// polar coordinates of the vertex
  long double mu;			// 'mass' of the induced charge
  long double w;			// discretized induced charge on the vertex
  long double vw;			// 'velocity' of the induced charge
  long double fw; 			// 'force' on the induced charge
  vector<long double> fwiwi; 			// 'force' on the induced charge
  vector<long double> fwiwj; 			// 'force' on the induced charge
  vector<long double> fwq_1; 			// 'force' on the induced charge
  vector<long double> fwq_2; 			// 'force' on the induced charge
  vector<long double> fqw_1; 			// 'force' on the induced charge
  vector<long double> fqw_2; 			// 'force' on the induced charge
  long double ke;              		// 'kinetic energy' of the induced charge
  double pe;				// potential energy of the fake degrees associated with vertices
  double wmean;				// mean w computed on fmd
  
  // member vectors
//  vector<long double> presumgwEw;		// result of a precalculation-- gwEw
//  vector<long double> presumgEwEq;		// result of a precalculation-- gEwEq
//  vector<long double> presumgEwEw;		// result of a precalculation-- gEwEw
//  vector<long double> presumfwEw;		// result of a precalculation-- fwEw
//  vector<long double> presumfEwEq;		// result of a precalculation-- fwEw
//  vector<long double> presumhEqEw;		// result of a precalculation-- hEqEw
  vector<long double> Greens;		// Greens functions between induced charges on same interface  Gw1w1 Gw2w2
  vector<long double> Greens_inter;		// Greens functions between induced charges on different interface Gw1w2 Gw2w1
  vector<long double> Greens_bar_1;      // precalculate G_bar_1_ w1w1,G_bar_1_ w2w2
  vector<long double> Greens_bar_1_inter;      //G_bar_1_ w1w2,G_bar_1_ w2w1
  vector<long double> Greens_bar_2;      //G_bar_2_ w1w1,G_bar_2_ w2w2
  vector<long double> Greens_bar_2_inter;      //G_bar_2_ w1w2,G_bar_2_ w2w1
  vector<long double> Greens_bbar_11;  //G_barbar_11_w1w1, G_barbar_11_w2w2
  vector<long double> Greens_bbar_11_inter;  //G_barbar_11_w1w2, G_barbar_11_w2w1
  vector<long double> Greens_bbar_22;  //G_barbar_22_w1w1, G_barbar_22_w2w2
  vector<long double> Greens_bbar_22_inter;  //G_barbar_22_w1w2, G_barbar_22_w2w1
  vector<long double> Greens_bbar_12;  //G_barbar_12_w1w1, G_barbar_12_w2w2
  vector<long double> Greens_bbar_12_inter;  //G_barbar_12_w1w2, G_barbar_12_w2w1
  vector<long double> ndotGradGreens;	// precalculate n.gradGwiwi
  vector<long double> ndotGradGreens_inter;	// precalculate n.gradGw1w2
  vector<long double> Charged_sheet_corrections_intra;	// precalculate n.gradGw1w2
  vector<long double> Charged_sheet_corrections_inter;	// precalculate n.gradGw1w2
//  long double ndot_grad_greens_innersum; //insum used in fake & real forces calculations
//  long double greens_innersum; //insum used in fake & real forces calculations
  vector<long double> inner_wiq_i;	//insum used in real force calculations 
  vector<long double> inner_wjq_i;	//insum used in real force calculations 
  
  // member functions
  
  // make a vertex
  VERTEX(VECTOR3D position = VECTOR3D(0,0,0), double area = 0, VECTOR3D normal = VECTOR3D(0,0,0)) : posvec(position), a(area), normalvec(normal)
  {
  }
  
  // update position of fake degree
  void update_position(double dt)		
  {
    w = w + dt * vw;
    return;
  }
  
  // update velocity of fake degree
  void update_velocity(double dt)		
  {
    vw = vw + 0.5 * dt * fw / (mu);
    return;
  }
  
  void new_update_velocity(double dt, THERMOSTAT main_bath, long double expfac)
  {
//     vw = vw * expfac + 0.5 * dt * (fw / (mu)) * (expfac - 1) / (-0.5 * dt * main_bath.xi);
    vw = vw * expfac + 0.5 * dt * (fw / (mu)) * sqrt(expfac);
    return;
  }
  
  // compute kinetic energy of fake degree
  void kinetic_energy()				
  {
    ke = 0.5 * mu * vw * vw;
    return;
  }
  
  // get the polar coordinates for the vertex
  void get_polar() 				
  {
    r = posvec.GetMagnitude();
    theta = posvec.y > 0 ? acos(posvec.z / r) : 2 * pi - acos(posvec.z / r);
    phi = posvec.y > 0 ? acos(posvec.x / r) : - acos(posvec.x / r);
    return;
  }
};

#endif
