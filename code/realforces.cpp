// This file contains the routine that computes the force
// on the induced charge at vertex k and the force on the particle i
// for all k and i

#include "forces.h"

// Total Force on all degrees of freedom
void compute_force_real_dof(vector<PARTICLE>& ion, INTERFACE& box, char flag)
{

  // force calculation for ions (electrostatic)
  // hqqc : force due to central charge qc
  // hqq : force due to ion (q) - ion (q) interaction
  // hqw1 : force due to ion (q) - induced charge on left interface (w1) interaction
  // hqw2 : force due to ion (q) - induced charge on right interface (w2) interaction
  // part of hqEq : force due to ion (q) - electric field of ion (Eq) interaction
  // hqEq1: force due to ion (q) - electric field of ion  (Eq1) interaction
  // hqEq2 : force due to ion (q) - electric field of ion (Eq2) interaction
  // hqEw : force due to ion (q) - electric field of induced charge (Ew) interaction
  // other part of hqEq : force due to ion (q) - electric field of ion (Eq) interaction
  // hEq1w1 : force due to electric field of ion (Eq) - induced charge (w) interaction
  // hEq2w1: force due to electric field of ion (Eq) - induced charge (w) interaction
  // hEq1w2 : force due to electric field of ion (Eq) - induced charge (w) interaction
  // hEq2w2: force due to electric field of ion (Eq) - induced charge (w) interaction
  // hEq1Eq1 : force due to electric field of ion (Eq) - electric field  of ion (Eq) interaction
  // hEq2Eq2 : force due to electric field of ion (Eq) - electric field  of ion (Eq) interaction
  // hEq1Eq2 : force due to electric field of ion (Eq) - electric field  of ion (Eq) interaction
  // hEqEw : force due to electric field of ion (Eq) - electric field of induced degree (Ew) interaction
   vector<VERTEX> L = box.leftplane;
   vector<VERTEX> R = box.rightplane;

   unsigned int i, j, k, l;
   int factor;
   double hcsh,E_z, r_1, r_2, dz, gcsh_z, gcsh_inf;
   double insum1, insum2, innersum_1, innersum_2, innersum_3, innersum_4;
  // VECTOR3D insum_1, insum_2, insum_3, insum_4,insum_5, insum_6, insum_7, insum_8;
   VECTOR3D vec, grad_G, grad_ndot_G ;
   double left_gw1q, left_gw1q_1_gw1w1, left_gw1q_2_gw1w2;
   double right_gw2q, right_gw2q_1_gw2w1, right_gw2q_2_gw2w2;
   VECTOR3D hqq_ij, hqq_1_ij, hqq_2_ij, hqq_11_ij, hqq_22_ij, hqq_12_ij, hqq_21_ij;
   //VECTOR3D hqq_ji, hqq_1_ji, hqq_2_ji, hqq_11_ji, hqq_22_ji, hqq_12_ji, hqq_21_ji;
   vector<double> inner_qq_11(L.size(), 0.0);         	//precaculation of hEqEq_11
   vector<double> inner_qq_12(L.size(), 0.0);        	//precaculation of hEqEq_12
   vector<double> inner_qq_21(R.size(), 0.0);       	//precaculation of hEqEq_21
   vector<double> inner_qq_22(R.size(), 0.0);        	//precaculation of hEqEq_22

   //VECTOR3D hqw1, hqEw1_1, hqEw1_2, hw1Eq_1, hw1Eq_2,hEqEw1_11, hEqEw1_22, hEqEw1_12, hEw1Eq_12;
   VECTOR3D hqw1, hqEw1_1, hw1Eq_1;
   vector<double> inner_qw1_1(L.size(), 0.0);         	//precaculation of hqEw1_1
   vector<double> inner_w1q_1(L.size(), 0.0);         	//precaculation of hw1Eq_1

   //VECTOR3D hqw2, hqEw2_1, hqEw2_2, hw2Eq_1, hw2Eq_2,hEqEw2_11, hEqEw2_22, hEw2Eq_21, hEqEw2_21;
   VECTOR3D hqw2, hqEw2_2, hw2Eq_2;
   vector<double> inner_qw2_2(R.size(), 0.0);         	//precaculation of hqEw2_2
   vector<double> inner_w2q_2(R.size(), 0.0);          	//precaculation of hw2Eq_2


  vector<double> inner_lw1(L.size(),0.0);             //precaculation og gw1Eq_1
  vector<double> inner_lw2(R.size(),0.0);             //precaculation og gw1Eq_2
  vector<double> inner_lw3(L.size(),0.0);             //precaculation og gEqw1_1
  vector<double> inner_lw4(R.size(),0.0);             //precaculation og gEqw1_2

  vector<VECTOR3D> lj_ion_ion(ion.size(),VECTOR3D(0,0,0));
  vector<VECTOR3D> lj_ion_leftdummy(ion.size(),VECTOR3D(0,0,0));
  vector<VECTOR3D> lj_ion_left_wall(ion.size(),VECTOR3D(0,0,0));
  vector<VECTOR3D> lj_ion_rightdummy(ion.size(),VECTOR3D(0,0,0));
  vector<VECTOR3D> lj_ion_right_wall(ion.size(),VECTOR3D(0,0,0));
  VECTOR3D fljcc, flj;
  PARTICLE dummy, wall_dummy;
  //double r2, d, d2, elj, r6, r12,d6,d12;
  #pragma omp parallel default(shared) private(k, j, i, l, vec, insum1, insum2, dz, r_1, r_2, gcsh_z, gcsh_inf, factor, hcsh, E_z, \
  left_gw1q, left_gw1q_1_gw1w1, left_gw1q_2_gw1w2, right_gw2q, right_gw2q_1_gw2w1, right_gw2q_2_gw2w2,\
  grad_G, grad_ndot_G,innersum_1, innersum_2, innersum_3, innersum_4,\
  hqq_ij, hqq_1_ij, hqq_2_ij, hqq_11_ij, hqq_22_ij, hqq_12_ij, hqq_21_ij,\
  hqw1, hqEw1_1, hw1Eq_1, hqw2, hqEw2_2, hw2Eq_2, fljcc, flj, dummy, wall_dummy) num_threads(THREADSIZE)
  {  
    #pragma omp for schedule(dynamic)nowait
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

    #pragma omp for schedule(dynamic)nowait
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
    #pragma omp barrier

    #pragma omp for schedule(dynamic)nowait
    for (k = 0; k < L.size(); k++)
    {
      innersum_1 = 0;
      for(l = 0; l < L.size(); l ++)
	innersum_1 += L[k].Greens[l] * inner_lw1[l] * L[l].a;
      inner_qq_11[k] = innersum_1;//precaculation of hEqEq_11

      innersum_2 = 0;
      for(l = 0; l < R.size(); l ++)
	innersum_2 += L[k].Greens_inter[l] * inner_lw2[l] *  R[l].a;
      inner_qq_12[k] = innersum_2 ;//precaculation of hEqEq_12


      //pre calculation of ion - leftplane interaction
      innersum_1 = 0;
      innersum_2 = 0;
      for (l = 0; l < L.size(); l++)
	{

	innersum_1 += L[k].ndotGradGreens[l] * L[l].w * L[l].a;
	innersum_2 += L[k].inner_wiq_i[l] * L[l].w;

	}

      //pre calculation of ion - rightplane interaction
      innersum_3 = 0;
      innersum_4 = 0;
      for (l = 0; l < R.size(); l++)
	{

	innersum_3 += L[k].ndotGradGreens_inter[l] * R[l].w * R[l].a;
	innersum_4 += L[k].inner_wjq_i[l] * R[l].w;

	}
      inner_qw1_1[k] = (innersum_1 + innersum_3) * (-0.5 * box.ed_left); //hqEw1_1
      inner_w1q_1[k] = (innersum_2 + innersum_4); //hqEw1_1

      //calculate force on vertexes on left plane

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
      //force from interfaces
      left_gw1q_1_gw1w1= 0;                      
      for (l = 0; l < L.size(); l++)
        left_gw1q_1_gw1w1 += L[k].fwq_1[l] * inner_lw1[l] + L[k].fqw_1[l] * inner_lw3[l] + L[k].fwiwi[l] * L[l].w;

      left_gw1q_2_gw1w2 = 0;
      for (l = 0; l < R.size(); l++)
        left_gw1q_2_gw1w2 += L[k].fwq_2[l] * inner_lw2[l] + L[k].fqw_2[l] * inner_lw4[l] + L[k].fwiwj[l] * R[l].w;

      L[k].fw = left_gw1q + left_gw1q_1_gw1w1 + left_gw1q_2_gw1w2;
    }

    #pragma omp for schedule(dynamic)nowait
    for (k = 0; k < R.size(); k++)
    {
      innersum_1 = 0;
      for(l = 0; l < L.size(); l ++)
	innersum_1 += R[k].Greens_inter[l] *  inner_lw1[l] *  L[l].a;
      inner_qq_21[k] = innersum_1;//precaculation of hEqEq_21

      innersum_2 = 0;
      for(l = 0; l < R.size(); l ++)
	innersum_2 += R[k].Greens[l] * inner_lw2[l] *  R[l].a;
      inner_qq_22[k] = innersum_2;//precaculation of hEqEq_22


      //pre calculation of ion - leftplane interaction
      innersum_1 = 0;
      innersum_2 = 0;
      for (l = 0; l < L.size(); l++)
	{
	innersum_1 += R[k].ndotGradGreens_inter[l] * L[l].w * L[l].a;
	innersum_2 += R[k].inner_wjq_i[l] * L[l].w;
	}

      //pre calculation of ion - rightplane interaction
      innersum_3 = 0;
      innersum_4 = 0;
      for (l = 0; l < R.size(); l++)
	{
	innersum_3 += R[k].ndotGradGreens[l]* R[l].w * R[l].a;
	innersum_4 += R[k].inner_wiq_i[l] * R[l].w;
	}
      inner_qw2_2[k] = (innersum_1 + innersum_3) * ( -0.5 * box.ed_right);//hqEw2_2
      inner_w2q_2[k] = (innersum_2 + innersum_4);//hqEw2_2

      //calculate force on vertices on right plane

      //force from free ions 
      right_gw2q=0;			//Gwq
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

      //force from interfaces
      right_gw2q_1_gw2w1 =0;                  
      for (l = 0; l < L.size(); l++)
        right_gw2q_1_gw2w1 += R[k].fwq_1[l]* inner_lw1[l] + R[k].fqw_1[l] * inner_lw3[l] + R[k].fwiwj[l] * L[l].w;

      right_gw2q_2_gw2w2 =0;
      for (l = 0; l < R.size(); l++)
        right_gw2q_2_gw2w2 += R[k].fwq_2[l]* inner_lw2[l] + R[k].fqw_2[l] * inner_lw4[l] + R[k].fwiwi[l] * R[l].w;

      R[k].fw = right_gw2q + right_gw2q_1_gw2w1 + right_gw2q_2_gw2w2;
    }

  //pre calculation of ion (q) - ion (q) interaction
    #pragma omp barrier

    #pragma omp for schedule(dynamic)
    for(i = 0; i < ion.size(); i++)
    {
      hqq_ij = VECTOR3D(0,0,0);

      for (j = 0; j < ion.size(); j++)
      {
	if (i == j) continue;

	vec = ion[i].posvec-ion[j].posvec;
	vec = Mini_Image(vec, box.lx);
	//--------- charged sheet method for pbc
     	dz = vec.z;
        if (dz >= 0) factor = 1;
        else factor = -1;
        r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
        r_2 = (0.25 + (dz/box.lx)*(dz/box.lx));
        E_z = 4 * atan(4*fabs(dz)*r_1/box.lx);
        hcsh = (4 / box.lx ) * ( 1/(r_1 *(0.5 + r_1)) - 1/(r_2)) * dz + factor * E_z + 16 * fabs(dz) * (box.lx / ( box.lx * box.lx + 4 * dz * dz * r_1 * r_1)) * ( fabs(dz) * dz / (box.lx * box.lx * r_1) + factor * r_1);

      	hqq_ij.z += 2 * ion[i].q * (ion[j].q/ (box.lx * box.lx)) * 0.5 * (1/ion[i].epsilon + 1/ion[j].epsilon) * hcsh;
	//------------
	
	hqq_ij += Grad_Greens(vec) ^(-0.5 * ion[i].q * ion[j].q /ion[i].epsilon); // hqq_ji = hqq_ij;
      }

      
      hqq_1_ij = VECTOR3D(0,0,0);
      hqq_11_ij = VECTOR3D(0,0,0);
      hqq_12_ij = VECTOR3D(0,0,0);

      hqw1 = VECTOR3D(0,0,0);
      hqEw1_1 = VECTOR3D(0,0,0);
      hw1Eq_1 = VECTOR3D(0,0,0);

      for (k = 0; k < L.size(); k ++)
      {
	vec = ion[i].posvec - L[k].posvec;
	vec = Mini_Image(vec, box.lx);

	grad_G = Grad_Greens(vec) ^ ((ion[i].q / ion[i].epsilon) * L[k].a);
	grad_ndot_G = Grad_ndot_Greens(L[k].normalvec, vec) ^ ((ion[i].q / ion[i].epsilon) * L[k].a);
	
	hqq_1_ij += (grad_G ^ inner_lw1[k]) + (grad_ndot_G ^ inner_lw3[k]);
	hqq_11_ij += grad_ndot_G ^ inner_qq_11[k];
	hqq_12_ij += grad_ndot_G ^ inner_qq_12[k];

        //hqq_11_ji = hqq_11_ij;
        //hqq_12_ji = hqq_21_ij;

        //---------charged sheet method for pbc
        dz = vec.z;
        if (dz >= 0) factor = 1;
        else factor = -1;
        r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
        r_2 = (0.25 + (dz/box.lx)*(dz/box.lx));
        E_z = 4 * atan(4*fabs(dz)*r_1/box.lx);
        hcsh = (4 / box.lx ) * ( 1/(r_1 *(0.5 + r_1)) - 1/(r_2)) * dz + factor * E_z + 16 * fabs(dz) * (box.lx / ( box.lx * box.lx + 4 * dz * dz * r_1 * r_1)) * ( fabs(dz) * dz / (box.lx * box.lx * r_1) + factor * r_1);

        hqw1.z += -2 * (ion[i].q / ion[i].epsilon) * (L[k].w * L[k].a/ (box.lx * box.lx)) * hcsh;
	//-----------
	
	hqw1 += grad_G ^ L[k].w;
	hqEw1_1 += grad_G ^ inner_qw1_1[k];
	hw1Eq_1 += grad_ndot_G ^ inner_w1q_1[k];

      }

      hqq_1_ij = hqq_1_ij ^ (-0.5 * box.ed_left);
      hqq_11_ij = hqq_11_ij ^ (-0.5 * box.ed_left * box.ed_left) ;
      hqq_12_ij = hqq_12_ij ^ (-0.5 * box.ed_left * box.ed_right) ;
      hqw1 = hqw1 ^ ( -0.5 * (1 - box.em_left/ion[i].epsilon) * ion[i].epsilon );

      hqq_2_ij = VECTOR3D(0,0,0);
      hqq_21_ij = VECTOR3D(0,0,0);
      hqq_22_ij = VECTOR3D(0,0,0);

      hqw2 = VECTOR3D(0,0,0);
      hqEw2_2 = VECTOR3D(0,0,0);
      hw2Eq_2 = VECTOR3D(0,0,0);

      for (k = 0; k < R.size(); k ++)
      {
	vec = ion[i].posvec - R[k].posvec;
	vec = Mini_Image(vec, box.lx);

	grad_G = Grad_Greens(vec) ^ ((ion[i].q / ion[i].epsilon) * R[k].a);
	grad_ndot_G = Grad_ndot_Greens(R[k].normalvec, vec) ^ ((ion[i].q / ion[i].epsilon) * R[k].a);

	
	hqq_2_ij += (grad_G ^ inner_lw2[k]) + (grad_ndot_G ^ inner_lw4[k]);
	hqq_21_ij += grad_ndot_G ^ inner_qq_21[k];
	hqq_22_ij += grad_ndot_G ^ inner_qq_22[k];

        //hqq_21_ji = hqq_12_ij;
        //hqq_22_ji = hqq_22_ij;

        //---------charged sheet method for pbc
        dz = vec.z;
        if (dz >= 0) factor = 1;
        else factor = -1;
        r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
        r_2 = (0.25 + (dz/box.lx)*(dz/box.lx));
        E_z = 4 * atan(4*fabs(dz)*r_1/box.lx);
        hcsh = (4 / box.lx ) * ( 1/(r_1 *(0.5 + r_1)) - 1/(r_2)) * dz + factor * E_z + 16 * fabs(dz) * (box.lx / ( box.lx * box.lx + 4 * dz * dz * r_1 * r_1)) * ( fabs(dz) * dz / (box.lx * box.lx * r_1) + factor * r_1);

        hqw2.z += -2 * (ion[i].q / ion[i].epsilon) * (R[k].w * R[k].a/ (box.lx * box.lx)) * hcsh;
	//--------------------

	hqw2 += grad_G ^ R[k].w;
	hqEw2_2 += grad_G ^ inner_qw2_2[k];
	hw2Eq_2 += grad_ndot_G ^ inner_w2q_2[k];
      }

      hqq_2_ij = hqq_2_ij ^ (-0.5 * box.ed_right);
      hqq_21_ij = hqq_21_ij ^ (-0.5 * box.ed_left * box.ed_right) ;
      hqq_22_ij = hqq_22_ij ^ (-0.5 * box.ed_right * box.ed_right) ;
      hqw2 = hqw2 ^ (-0.5 * (1 - box.em_right/ ion[i].epsilon) * ion[i].epsilon);

      ion[i].forvec = hqq_ij + hqq_1_ij + hqq_2_ij + hqq_11_ij + hqq_22_ij + hqq_12_ij + hqq_21_ij \
	  +  hqq_ij +  hqq_11_ij+ hqq_22_ij + hqq_12_ij + hqq_21_ij \
	  + hqw1 + hqEw1_1 + hw1Eq_1 + hqw2 + hqEw2_2 + hw2Eq_2;
  
  }


  // excluded volume interactions given by purely repulsive LJ
  // ion-ion
  #pragma omp for schedule(dynamic) nowait
  for ( i = 0; i < ion.size(); i++)
  {
    fljcc = VECTOR3D(0,0,0);
    for (j = 0; j < ion.size(); j++)
    {
      if (j == i) continue;
      vec = ion[i].posvec - ion[j].posvec;
      vec = Mini_Image(vec, box.lx);
      double r2 = vec.GetMagnitudeSquared();
      double d = 0.5 * (ion[i].diameter + ion[j].diameter);
      double d2 = d * d;
      double elj = 1.0;
      if (r2 < dcut2 * d2)
      {
        double r6 = r2 * r2 * r2;
        double r12 = r6 * r6;
        double d6 = d2 * d2 * d2;
        double d12 = d6 * d6;
        fljcc = fljcc + ( vec ^ ( 48 * elj * (  (d12 / r12)  - 0.5 *  (d6 / r6) ) * ( 1 / r2 ) ) );
      }
      else
        fljcc = fljcc + VECTOR3D(0,0,0);
    }
    lj_ion_ion[i] = fljcc;

  // ion-box

  // interaction with the left plane hard wall

  // make a dummy particle with the same diameter as the ion and touching left of the left wall s. t. it is closest to the ion
    flj = VECTOR3D(0,0,0);
    if(ion[i].posvec.z < -0.5*box.lz + ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk. replacing 1 by diameter. -Yufei -Vikram -Vikram
    {
      dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_left,
      VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, -0.5*box.lz - 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
      vec = ion[i].posvec - dummy.posvec;
      double r2 = vec.GetMagnitudeSquared();
      double d = 0.5 * (ion[i].diameter + dummy.diameter);
      double d2 = d * d;
      double elj = 1.0;
      if (r2 < dcut2 * d2)
      {
        double r6 = r2 * r2 * r2;
        double r12 = r6 * r6;
        double d6 = d2 * d2 * d2;
        double d12 = d6 * d6;
        flj = vec ^ ( 48 * elj * ( (d12 / r12) - 0.5 * (d6 / r6) ) * (1 / r2) );
      }
    }
    lj_ion_leftdummy[i]=flj;

  // ion interacting with discretized left wall
    flj = VECTOR3D(0,0,0);
    if(ion[i].posvec.z < -0.5*box.lz + ion[i].diameter)  // avoiding calculating interactions between left wall and ions in bulk. -Yufei - Vikram
    {
      for (k = 0; k < box.leftplane.size(); k++)
      {
        wall_dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_left,
        VECTOR3D(box.leftplane[k].posvec.x, box.leftplane[k].posvec.y, box.leftplane[k].posvec.z - 0.5*ion[i].diameter),
                                      box.lx,box.ly,box.lz);
        vec = ion[i].posvec - wall_dummy.posvec;
        vec = Mini_Image(vec, box.lx);
        double r2 = vec.GetMagnitudeSquared();
        double d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
        double d2 = d * d;
        double elj = 1.0;
        if (r2 < dcut2 * d2)
        {
          double r6 = r2 * r2 * r2;
          double r12 = r6 * r6;
          double d6 = d2 * d2 * d2;
          double d12 = d6 * d6;
          flj = flj + (vec ^ ( 48 * elj * ( (d12 / r12) - 0.5 * (d6 / r6) ) * (1 / r2) ));
        }
      }
    }
    lj_ion_left_wall[i]=flj;

  // interaction with the right plane hard wall

  // make a dummy particle with the same diameter as the ion and touching right of the right wall s. t. it is closest to the ion
    flj = VECTOR3D(0,0,0);
    if(ion[i].posvec.z > 0.5*box.lz - ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
    {
      dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_right,
      VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, 0.5*box.lz + 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
      vec = ion[i].posvec - dummy.posvec;
      double r2 = vec.GetMagnitudeSquared();
      double d = 0.5 * (ion[i].diameter + dummy.diameter);
      double d2 = d * d;
      double elj = 1.0;
      if (r2 < dcut2 * d2)
      {
        double r6 = r2 * r2 * r2;
        double r12 = r6 * r6;
        double d6 = d2 * d2 * d2;
        double d12 = d6 * d6;
        flj = vec ^ ( 48 * elj * ( (d12 / r12) - 0.5 * (d6 / r6) ) * (1 / r2) );
      }
    }
    lj_ion_rightdummy[i]=flj;

  // ion interacting with discretized right wall
    flj = VECTOR3D(0,0,0);
    if(ion[i].posvec.z > 0.5*box.lz - ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
    {
      for (k = 0; k < box.rightplane.size(); k++)
      {
        wall_dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_right,VECTOR3D(box.rightplane[k].posvec.x, box.rightplane[k].posvec.y, box.rightplane[k].posvec.z + 0.5*ion[i].diameter), box.lx,box.ly,box.lz);
        vec = ion[i].posvec - wall_dummy.posvec;
        vec = Mini_Image(vec, box.lx);
        double r2 = vec.GetMagnitudeSquared();
        double d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
        double d2 = d * d;
        double elj = 1.0;
        if (r2 < dcut2 * d2)
        {
          double r6 = r2 * r2 * r2;
          double r12 = r6 * r6;
          double d6 = d2 * d2 * d2;
          double d12 = d6 * d6;
          flj = flj + (vec ^ ( 48 * elj * ( (d12 / r12) - 0.5 * (d6 / r6) ) * (1 / r2) ));
        }

      }
    }
    lj_ion_right_wall[i]=flj;
  }

}

  for ( i = 0; i < ion.size(); i++)
    ion[i].forvec = ( (ion[i].forvec) ^ (scalefactor) );
  // t tal force on the particle = the electrostatic force + the Lennard-Jones force
  for (unsigned int i = 0; i < ion.size(); i++)
    ion[i].forvec = ion[i].forvec + lj_ion_ion[i] + lj_ion_leftdummy[i] + lj_ion_left_wall[i] + lj_ion_rightdummy[i] + lj_ion_right_wall[i];

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

   inner_qq_11.clear();           //precaculation of hEqEq_11
   inner_qq_12.clear();           //precaculation of hEqEq_12
   inner_qq_21.clear();           //precaculation of hEqEq_21
   inner_qq_22.clear();           //precaculation of hEqEq_22

   inner_qw1_1.clear();           //precaculation of hqEw1_1
   inner_w1q_1.clear();           //precaculation of hw1Eq_1

   inner_qw2_2.clear();           //precaculation of hqEw2_2
   inner_w2q_2.clear();           //precaculation of hw2Eq_2


  return;
}
