// This file contains the routine that computes the force
// on the induced charge at vertex k and the force on the particle i
// for all k and i

#include "forces.h"

// Total Force on all degrees of freedom: fictitious force on induced charge on left & right interface
void compute_force_fake_dof(vector<PARTICLE>& ion,  INTERFACE& box,  char flag)
{

  //gwiwj:			G (s_i, s_j)
  //gwiEwj_k:			Gbar k (s_i, s_j)
  //gEwiEwj_kl:			Gbarbar k,l (s_i, s_j)
  //gqq				G (r_i, r_j)
  //gqEq_k			Gbar k (r_i, r_j)
  //gEqEq_kl			Gbarbar k,l (r_i , r_j )
  //gwiq			G (s_i, q_j)
  //gqEwi_k			Gbar k (r_j, s_i)
  //gwiEq_k			Gbar k (s_i, r_j)
  //gEqEwi_kl			Gbarbar k, l (r_j, s_i)

  // force calculation for fake degrees of freedom of vertex on left plane
  // gw1q : force due to induced charge on left interface (w1) - ion (q) interaction
  // gw1w1 : force due to induced charge on left interface (w1) - induced charge on left interface (w1) interaction
  // gw1w2 : force due to induced charge on left interface (w1) - induced charge on right interface (w2) interaction
  // gw2w1 : force due to induced charge on right interface (w2) - induced charge on left interface (w1) interaction
  // gwEw : force due to induced charge (w) - Electric field due to induced charge (Ew) interaction
  // gEwEw : force due to Electric field due to induced charge (Ew) - Electric field due to induced charge (Ew) interaction
  // gEwq : force due to electric field of induced charge (Ew) - ion (q) interaction
  // gw1Eq : force due to induced charge on left interface (w1) - Electric field due to ion (Eq) interaction
  // gEwEq : force due to Electric field due to induced charge (Ew) - Electric field due to ion (Eq) interaction


  // force calculation for fake degrees of freedom of vertex on right  plane
  // gw2q : force due to induced charge on right interface (w2) - ion (q) interaction
  // gw2w2 : force due to induced charge on right interface (w2) - induced charge on right  interface (w2) interaction
  // gw1w2 : force due to induced charge on left interface (w1) - induced charge on right interface (w2) interaction
  // gw2w1 : force due to induced charge on right interface (w2) - induced charge on left interface (w1) interaction
  // gwEw : force due to induced charge (w) - Electric field due to induced charge (Ew) interaction,
  // gEwEw : force due to Electric field due to induced charge (Ew) - Electric field due to induced charge (Ew) interaction
  // gEwq : force due to electric field of induced charge (Ew) - ion (q) interaction
  // gw2Eq : force due to induced charge on right interface (w2) - Electric field due to ion (Eq) interaction
  // gEwEq : force due to Electric field due to induced charge (Ew) - Electric field due to ion (Eq) interaction  // declarations (necessary beforehand for parallel implementation)
  
  vector<VERTEX> L = box.leftplane;
  vector<VERTEX> R = box.rightplane;

  unsigned int i, k,l;
  double insum1, insum2;
  double dz, r_1, r_2, gcsh_z, gcsh_inf;
  VECTOR3D vec;
  double left_gw1q, left_gw1q_1_gw1w1, left_gw1q_2_gw1w2;
  double right_gw2q, right_gw2q_1_gw2w1, right_gw2q_2_gw2w2;

  vector<double> inner_lw1(L.size(),0.0);             //precaculation og gw1Eq_1
  vector<double> inner_lw2(R.size(),0.0);             //precaculation og gw1Eq_2
  vector<double> inner_lw3(L.size(),0.0);             //precaculation og gEqw1_1
  vector<double> inner_lw4(R.size(),0.0);             //precaculation og gEqw1_2

  #pragma omp parallel default(shared) private(k, i, l, vec, insum1, insum2, dz, r_1, r_2, gcsh_z, gcsh_inf, \
    left_gw1q, left_gw1q_1_gw1w1, left_gw1q_2_gw1w2, right_gw2q, right_gw2q_1_gw2w1, right_gw2q_2_gw2w2)
  {  
    #pragma omp for schedule(dynamic)
    for (k = 0; k < L.size(); k++)
    {
      insum1 = 0;
      insum2 = 0;

      for (i = 0; i < ion.size(); i++)
      { 
	vec = L[k].posvec - ion[i].posvec;
	vec = Mini_Image(vec, box.lx);
      
       	insum1 += (L[k].normalvec * Grad_Greens(vec))*ion[i].q /ion[i].epsilon;
      	insum2 += Greens(vec)*ion[i].q /ion[i].epsilon;
      }
    
      inner_lw1[k] = insum1; 		//gw1Eq_1				inner_lw1=inner_rw1
      inner_lw3[k] = insum2; 		//gqEw1_1				inner_lw1=inner_rw1
    }

    #pragma omp for schedule(dynamic)
    for (k = 0; k < R.size(); k++)
    {
      insum1 = 0;
      insum2 = 0;

      for (i = 0; i < ion.size(); i++)
      { 
	vec = R[k].posvec - ion[i].posvec; 
	vec = Mini_Image(vec, box.lx);

	insum1 += (R[k].normalvec * Grad_Greens(vec))*ion[i].q /ion[i].epsilon;
	insum2 += Greens(vec)*ion[i].q /ion[i].epsilon;
      }

      inner_lw2[k] = insum1; 		//gw1Eq_1				inner_lw1=inner_rw1
      inner_lw4[k] = insum2; 		//gqEw1_1				inner_lw1=inner_rw1
    }

    //calculate force on vertexes on left plane
    
    #pragma omp for schedule(dynamic) nowait
    for (k = 0; k < L.size(); k++)
    {
      //force from free ions 
      left_gw1q=0;					//Gwq
      for (i = 0; i < ion.size(); i++)
      { 
	vec = L[k].posvec - ion[i].posvec; 
	vec = Mini_Image(vec, box.lx);

        //---------charged sheet method for pbc
     	dz = vec.z;
     	r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
     	r_2 = (0.25 + (dz/box.lx)*(dz/box.lx));
     	gcsh_z = 2 * box.lx * log((0.5 + r_1)*(0.5 + r_1)/r_2) - fabs(dz) * (2 * pi - 4 * atan(4*fabs(dz)*r_1/box.lx));
     	gcsh_inf = -2 * pi * fabs(dz);

     	left_gw1q += (ion[i].q / (box.lx * box.lx)) * (1 - box.em_left/ion[i].epsilon)* (gcsh_inf - gcsh_z);
	//--------------
	left_gw1q += (-0.5) * ion[i].q * (1 - box.em_left/ion[i].epsilon)* Greens(vec);
      }

      left_gw1q_1_gw1w1= 0;                          //G bar for R_rho omega1 (left interface)
      for (l = 0; l < L.size(); l++)
        left_gw1q_1_gw1w1 += L[k].fwq_1[l] * inner_lw1[l] + L[k].fqw_1[l] * inner_lw3[l] + L[k].fwiwi[l] * L[l].w;

      left_gw1q_2_gw1w2 = 0;
      for (l = 0; l < R.size(); l++)
        left_gw1q_2_gw1w2 += L[k].fwq_2[l] * inner_lw2[l] + L[k].fqw_2[l] * inner_lw4[l] + L[k].fwiwj[l] * R[l].w;

      L[k].fw = left_gw1q + left_gw1q_1_gw1w1 + left_gw1q_2_gw1w2;
    }
  
    //calculate force on vertices on right plane
    
    #pragma omp for schedule(dynamic) nowait
    for (k = 0; k < R.size(); k++)
    {
      //force from free ions 
      right_gw2q=0;			//G
      for (i = 0; i < ion.size(); i++)
      { 
	vec = R[k].posvec - ion[i].posvec; 
	vec = Mini_Image(vec, box.lx);

        //---------charged sheet method for pbc
     	dz = vec.z;
     	r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
     	r_2 = (0.25 + (dz/box.lx)*(dz/box.lx));
     	gcsh_z = 2 * box.lx * log((0.5 + r_1)*(0.5 + r_1)/r_2) - fabs(dz) * (2 * pi - 4 * atan(4*fabs(dz)*r_1/box.lx));
     	gcsh_inf = -2 * pi * fabs(dz);

     	right_gw2q += (ion[i].q / (box.lx * box.lx)) * (1 - box.em_right/ion[i].epsilon)* (gcsh_inf - gcsh_z);
	//-----------------
	right_gw2q += (-0.5) * ion[i].q * (1 - box.em_right/ion[i].epsilon)*Greens(vec);
      }


      right_gw2q_1_gw2w1 =0;                                 //G bar for R_rho omega2 (right interface)
      for (l = 0; l < L.size(); l++)
        right_gw2q_1_gw2w1 += R[k].fwq_1[l]* inner_lw1[l] + R[k].fqw_1[l] * inner_lw3[l] + R[k].fwiwj[l] * L[l].w;

      right_gw2q_2_gw2w2 =0;
      for (l = 0; l < R.size(); l++)
        right_gw2q_2_gw2w2 += R[k].fwq_2[l]* inner_lw2[l] + R[k].fqw_2[l] * inner_lw4[l] + R[k].fwiwi[l] * R[l].w;

      R[k].fw = right_gw2q + right_gw2q_1_gw2w1 + right_gw2q_2_gw2w2;
    }
  }
  
  // force on the fake degrees of freedom
  for (unsigned int k = 0; k < L.size(); k++)
    box.leftplane[k].fw =  L[k].a * L[k].fw * scalefactor;
  // force on the fake degrees of freedom
  for (unsigned int k = 0; k < R.size(); k++)
    box.rightplane[k].fw = R[k].a * R[k].fw * scalefactor;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  inner_lw1.clear();
  inner_lw2.clear();
  inner_lw3.clear();
  inner_lw4.clear();

  return; 
}
