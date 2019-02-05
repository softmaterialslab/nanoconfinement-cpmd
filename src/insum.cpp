#include "vertex.h"
#include "particle.h"
#include "interface.h"
#include "functions.h"

void insum_forces(vector<PARTICLE>& ion,INTERFACE& box)
{

  vector<VERTEX> L = box.leftplane;
  vector<VERTEX> R = box.rightplane;

  unsigned int i,k,l;
  double temp, insum1, insum2;
  VECTOR3D vec;

  vector<double> inner_l1(L.size(),0.0);            
  vector<double> inner_l2(L.size(),0.0);         
  vector<double> inner_r1(R.size(),0.0);        
  vector<double> inner_r2(R.size(),0.0);       

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

     inner_l1[k] = insum1;             //gw1Eq_1                               inner_lw1=inner_rw1
     inner_l2[k] = insum2;             //gqEw1_1                               inner_lw1=inner_rw1

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

      inner_r1[k] = insum1;            //gw1Eq_2
      inner_r2[k] = insum2;            //gw1Eq_2
    }

    return;
}
