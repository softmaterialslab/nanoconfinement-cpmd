// This file contains the routine to compute the total potential energy of the system
// Electrostatic and Excluded volume contributions

#include "energies.h"

// Potential energy
long double energy_functional(vector<PARTICLE>& ion, INTERFACE& box)
{
  // Electrostatic interaction

  //fwiwj:                             G (s_i, s_j)
  //fwiEwj_k:                     Gbar k (s_i, s_j)
  //fEwiEwj_kl:                 Gbarbar k,l (s_i, s_j)
  //fqq                                    G (r_i, r_j)
  //fqEq_k                            Gbar k (r_i, r_j)
  //fEqEq_kl                       Gbarbar k,l (r_i , r_j )
  //fwiq                                 G (s_i, q_j)
  //fqEwi_k                         Gbar k (r_j, s_i)
  //fwiEq_k                         Gbar k (s_i, r_j)
  //fEqEwi_kl                     Gbarbar k, l (r_j, s_i)

  /*
fw1w1,fw1w2,fw2w1,fw2w2 : potential energy due to induced charge and induced charge interaction
  // fwEw : potential energy due to induced charge and E of induced charge interaction
  // fEwEw : potential due to E of induced charge and E of induced charge interaction
  // fqqc : potential energy due to ion and central charge
  // fqq : potential energy due to ion and ion interaction
  // fw1q,fw2q : potential energy due to induced charge and ion interaction
  // fqEq1,fqEq2 : potential due to ion and E due to ion interaction integrated on interface 1(left) or 2(right)
  // fEwq : potential energy due to E of induced charge and ion interaction, 0 here
  // fw1Eq1,fw1Eq2,fw2Eq1,fw2Eq2 : potential energy due to induced charge on interface 1(left) or 2(right) and E of ion interaction  integrated on interface 1(left) or 2(right)
  // fEq1Eq1,fEq2Eq2,fEq1Eq2 : potential due to E of ion and E of ion interaction
  // fEwEq : potential due to E of induced charge and E of ion interaction,0  here
  */

  vector<VERTEX>L = box.leftplane;
  vector<VERTEX> R = box.rightplane;

  unsigned int i, j, k, l;
  double  ion_ion, ind_ind, ion_ind;
  double temp, greens, ndot_grad_greens, insum, insum_1, insum_2, insum_3, insum_4, insum_5, insum_6, insum_7, insum_8;
  VECTOR3D vec;
  double fw1w1, fw1Ew1_1, fw1Ew1_2, fEw1Ew1_11, fEw1Ew1_22, fEw1Ew1_12 ;
  double fw2w2, fw2Ew2_2, fw2Ew2_1, fEw2Ew2_22, fEw2Ew2_11, fEw2Ew2_12 ;
  double fw1w2,fw2w1, fw1Ew2_1, fw1Ew2_2, fw2Ew1_1, fw2Ew1_2, fEw1Ew2_11, fEw2Ew1_11, fEw1Ew2_22, fEw2Ew1_22, fEw1Ew2_12, fEw2Ew1_12 ;
  double fw1q, fw2q, fw1Eq_1, fw1Eq_2, fw2Eq_1, fw2Eq_2, fqEw1_1, fqEw1_2, fqEw2_1, fqEw2_2, fEqEw1_11, fEqEw1_22, fEqEw1_12, fEqEw1_21,  fEqEw2_11, fEqEw2_22, fEqEw2_12, fEqEw2_21;
  double fqq,  fqEq_1, fqEq_2, fEqEq_11, fEqEq_22, fEqEq_12;
  double dz, r_1, r_2, fcsh_z, fcsh_inf;
  double uljcc, ulj, ulj_wall, elj, r, r2, r6, d, d2, d6;
  PARTICLE dummy, wall_dummy;

  vector<double> innersum_qll(L.size(),0.0);;             //precaculation of fEqEq_11
  vector<double> inner_qll(L.size(),0.0);;             //precaculation of fEqEq_11
  vector<double> innersum_qrr(R.size(),0.0);;             //precaculation of fEqEq_11
  vector<double> inner_qrr(R.size(),0.0);;             //precaculation of fEqEq_11
  vector<double> innersum_qlr(L.size(),0.0);;             //precaculation of fEqEq_11
  vector<double> inner_qlr(R.size(),0.0);;             //precaculation of fEqEq_11
  vector<double> inner_ql(L.size(),0.0);;             //precaculation of fqEq_1
  vector<double> inner_qr(R.size(),0.0);;             //precaculation of fqEq_2

  vector<double> inner_lw1(L.size(),0.0);             //precaculation of fw1Eq_1
  vector<double> inner_lw2(R.size(),0.0);             //precaculation of fw1Eq_2
  vector<double> inner_lw3(L.size(),0.0);             //precaculation of fEqw1_1
  vector<double> inner_lw4(R.size(),0.0);             //precaculation of fEqw1_2

  vector<double> inner_rw1(L.size(),0.0);             //precaculation of fw2Eq_1
  vector<double> inner_rw2(R.size(),0.0);             //precaculation of fw2Eq_2
  vector<double> inner_rw3(L.size(),0.0);             //precaculation of fEqw2_1
  vector<double> inner_rw4(R.size(),0.0);             //precaculation of fEqw2_2

  vector<double> inner_lw_11(L.size(),0.0);             //precaculation of fEqEw1_11
  vector<double> inner_lw_12(L.size(),0.0);             //precaculation of fEqEw1_12
  vector<double> inner_lw_22(R.size(),0.0);             //precaculation of fEqEw1_22
  vector<double> inner_lw_21(R.size(),0.0);             //precaculation of fEqEw1_21

  vector<double> inner_rw_11(L.size(),0.0);             //precaculation of fEqEw2_11
  vector<double> inner_rw_12(L.size(),0.0);             //precaculation of fEqEw2_12
  vector<double> inner_rw_22(R.size(),0.0);             //precaculation of fEqEw2_22
  vector<double> inner_rw_21(R.size(),0.0);             //precaculation of fEqEw2_21

  vector<double>  ind_energy_left(L.size(),0.0);
  vector<double>  ind_energy_right(R.size(),0.0);
  vector<double>  ion_energy(ion.size(),0.0);
  vector<double>  ion_ind_energy(ion.size(),0.0);

  vector<double> lj_ion_ion(ion.size(),0.0);
  vector<double> lj_ion_leftdummy(ion.size(),0.0);
  vector<double> lj_ion_leftwall(ion.size(),0.0);
  vector<double> lj_ion_rightdummy(ion.size(),0.0);
  vector<double> lj_ion_rightwall(ion.size(),0.0);
  // push_back is not compatible with pragma, so we have initialize vector in this way.

#pragma omp parallel default(shared) private(k, l, j, i, temp, vec, greens, ndot_grad_greens, insum, insum_1, insum_2, insum_3, insum_4, insum_5, insum_6, insum_7, insum_8,\
  fw1w1, fw1Ew1_1, fw1Ew1_2, fEw1Ew1_11, fEw1Ew1_22, fEw1Ew1_12, fw2w2, fw2Ew2_2, fw2Ew2_1, fEw2Ew2_22, fEw2Ew2_11, fEw2Ew2_12, dz, r_1, r_2, fcsh_z, fcsh_inf, \
  fqq , fqEq_1 , fqEq_2 , fEqEq_11 , fEqEq_22 , fEqEq_12, fw1q, fw2q , fw1Eq_1 , fw1Eq_2 , fqEw1_1 , fqEw1_2 , fw2Eq_1 , fw2Eq_2 , fqEw2_1 , fqEw2_2 ,\
  fEqEw1_11 , fEqEw1_22 , fEqEw1_12 , fEqEw1_21 , fEqEw2_11 , fEqEw2_22 , fEqEw2_12 , fEqEw2_21, uljcc, ulj, ulj_wall, dummy, wall_dummy, elj, r, r2, r6, d, d2, d6)
{
  //calculate induced charge-induced charge electrostatic interacion, all renormalized greens functions are calculated in precal.cpp
  #pragma omp for schedule(dynamic) nowait
    for (k = 0; k < L.size(); k++)
      {
       fw1w1 = 0;
       fw1Ew1_1 = 0;
       fw1Ew1_2 = 0;
       fEw1Ew1_11 = 0;
       fEw1Ew1_22 = 0;
       fEw1Ew1_12 =0;
       for ( l = 0; l < L.size(); l++)
         {
	    temp = 0.5 * L[k].w* L[k].a * L[l].w * L[l].a;

            fw1w1 +=  temp * box.em_left * (box.em_left - 1) * L[k].Charged_sheet_corrections_intra[l];
            fw1w1 +=  temp * box.em_left * (box.em_left - 1) * L[k].Greens[l];
            fw1Ew1_1 +=  temp * box.ed_left * (1- 2*box.em_left ) * L[k].Greens_bar_1[l];
            fw1Ew1_2 += temp * box.ed_right * (1- 2*box.em_left ) * L[k].Greens_bar_2[l];
            fEw1Ew1_11 += temp * box.ed_left * box.ed_left * L[k].Greens_bbar_11[l];
            fEw1Ew1_22 += temp * box.ed_right * box.ed_right * L[k].Greens_bbar_22[l];
            fEw1Ew1_12 += temp * 2 * box.ed_left * box.ed_right * L[k].Greens_bbar_12[l];
         }


       fw1w2 = 0;
       fw1Ew2_1 = 0;
       fw1Ew2_2 = 0;
       fEw1Ew2_11 = 0;
       fEw1Ew2_22 = 0;
       fEw1Ew2_12 =0;
       for ( l = 0; l < R.size(); l++)
         {
	    temp = 0.5 * L[k].w* L[k].a * R[l].w * R[l].a;

            fw1w2 += temp * box.em_left * (box.em_right - 1) * L[k].Charged_sheet_corrections_inter[l];
            fw1w2 += temp * box.em_left * (box.em_right - 1) * L[k].Greens_inter[l];
            fw1Ew2_1 += temp * box.ed_left * (1- 2*box.em_left ) * L[k].Greens_bar_1_inter[l];
            fw1Ew2_2 += temp * box.ed_right * (1- 2*box.em_left ) * L[k].Greens_bar_2_inter[l];
            fEw1Ew2_11 += temp * box.ed_left * box.ed_left * L[k].Greens_bbar_11_inter[l];
            fEw1Ew2_22 += temp * box.ed_right * box.ed_right * L[k].Greens_bbar_22_inter[l];
            fEw1Ew2_12 += temp * 2 * box.ed_left * box.ed_right * L[k].Greens_bbar_12_inter[l];
         }

       ind_energy_left[k]=  fw1w1+ fw1Ew1_1+ fw1Ew1_2+  fEw1Ew1_11 + fEw1Ew1_22 + fEw1Ew1_12 + fw1w2+ fw1Ew2_1+ fw1Ew2_2+  fEw1Ew2_11 + fEw1Ew2_22 + fEw1Ew2_12;

      }

  #pragma omp for schedule(dynamic) nowait
    for (k = 0; k < R.size(); k++)
      {
       fw2w2 = 0;
       fw2Ew2_2 = 0;
       fw2Ew2_1 = 0;
       fEw2Ew2_22 = 0;
       fEw2Ew2_11 = 0;
       fEw2Ew2_12 = 0;
       for ( l = 0; l < R.size(); l++)
         {
	    temp = 0.5 * R[k].w* R[k].a * R[l].w * R[l].a;

            fw2w2 += temp * box.em_right * (box.em_right - 1) * R[k].Charged_sheet_corrections_intra[l];
            fw2w2 += temp * box.em_right * (box.em_right - 1) * R[k].Greens[l];
            fw2Ew2_2 += temp * box.ed_right* (1 - 2 * box.em_right ) * R[k].Greens_bar_2[l];
            fw2Ew2_1 += temp * box.ed_left* (1 - 2 * box.em_right ) * R[k].Greens_bar_1[l];
            fEw2Ew2_22 += temp * box.ed_right* box.ed_right * R[k].Greens_bbar_22[l];
            fEw2Ew2_11 += temp * box.ed_left* box.ed_left * R[k].Greens_bbar_11[l];
            fEw2Ew2_12 += temp *2 *  box.ed_left* box.ed_right * R[k].Greens_bbar_12[l];
         }

       fw2w1 = 0;
       fw2Ew1_2 = 0;
       fw2Ew1_1 = 0;
       fEw2Ew1_22 = 0;
       fEw2Ew1_11 = 0;
       fEw2Ew1_12 = 0;
       for ( l = 0; l < L.size(); l++)
         {
	    temp = 0.5 * R[k].w* R[k].a * L[l].w * L[l].a;

            fw2w1 += temp * box.em_right * (box.em_left - 1) * R[k].Charged_sheet_corrections_inter[l];
            fw2w1 += temp * box.em_right * (box.em_left - 1) * R[k].Greens_inter[l];
            fw2Ew1_2 += temp * box.ed_right* (1 - 2 * box.em_right ) * R[k].Greens_bar_2_inter[l];
            fw2Ew1_1 += temp * box.ed_left* (1 - 2 * box.em_right ) * R[k].Greens_bar_1_inter[l];
            fEw2Ew1_22 += temp * box.ed_right* box.ed_right * R[k].Greens_bbar_22_inter[l];
            fEw2Ew1_11 += temp * box.ed_left* box.ed_left * R[k].Greens_bbar_11_inter[l];
            fEw2Ew1_12 += temp *2 *  box.ed_left* box.ed_right * R[k].Greens_bbar_12_inter[l];
         }
       ind_energy_right[k] =  fw2w2+ fw2Ew2_2+fw2Ew2_1 + fEw2Ew2_22 + fEw2Ew2_11 + fEw2Ew2_12 + fw2w1+ fw2Ew1_2+fw2Ew1_1 + fEw2Ew1_22 + fEw2Ew1_11 + fEw2Ew1_12; //right interface
      }

  //precaculation for inner summation, commutative law of addition is applied. G() need one summation; G_bar need two summations; Gbarbar need three
  //this presum are suitable for calculation of forces on free ions. another presum are used for fictitious force on vertexes, see pfmdforces.cpp.

  #pragma omp for schedule(dynamic)
    for (k = 0; k < L.size(); k++)
      {
        insum = 0;
        for ( j = 0; j < ion.size(); j++)
	{ 
	  vec = L[k].posvec - ion[j].posvec;
	  vec = Mini_Image(vec, box.lx);

          insum +=  (L[k].normalvec * Grad_Greens(vec)) * L[k].a * ion[j].q / ion[j].epsilon;
	}
        innersum_qll[k] = insum;  //fEqEq_11
        inner_ql[k] = insum; //fqEq_1
      }

  #pragma omp for schedule(dynamic)
    for (k = 0; k < R.size(); k++)
      {
        insum = 0;
        for ( j = 0; j < ion.size(); j++)
	{ 
	  vec = R[k].posvec - ion[j].posvec;
	  vec = Mini_Image(vec, box.lx);

          insum += (R[k].normalvec * Grad_Greens(vec)) * R[k].a * ion[j].q / ion[j].epsilon;
	}
        innersum_qrr[k] = insum; // fEqEq_22
        innersum_qlr[k]  = insum; // fEqEq_12
        inner_qr[k]  = insum;
      }

  #pragma omp for schedule(dynamic)
    for (k = 0; k < L.size(); k++)
      {
        insum_1 = 0;
        for (l = 0; l < L.size(); l++)
          insum_1 += L[k].Greens[l] *innersum_qll[l] ;
       	inner_qll[k] = insum_1;  //fEqEq_11

       	insum_2 = 0;
       	for ( l = 0; l < R.size(); l++)
          insum_2 +=L[k].Greens_inter[l] *innersum_qlr[l] ;
       	inner_qlr[k] = insum_2; // fEqEq_12
      }

  #pragma omp for schedule(dynamic)
    for (k = 0; k < R.size(); k++)
      {
    	insum = 0;
    	for ( l = 0; l < R.size(); l++)
      	  insum += R[k].Greens[l] *innersum_qrr[l];
    	inner_qrr[k] = insum; // fEqEq_22
      }

  #pragma omp for schedule(dynamic)
    for (k = 0; k < L.size(); k++)
      {
        insum_1 = 0;
        for (l = 0; l < L.size(); l++)
          insum_1 += L[k].Greens[l] * L[l].w * L[l].a;
        inner_lw1[k] = insum_1; //fw1Eq_1

        insum_2 = 0;
        for (l = 0; l < R.size(); l++)
          insum_2 += L[k].Greens_inter[l] * R[l].w * R[l].a;
        inner_rw1[k] = insum_2; //fw2Eq_1

	insum_3 = 0;
	for (l = 0; l < L.size(); l++)
 	  insum_3 += L[k].ndotGradGreens[l] * L[l].w * L[l].a;
	inner_lw3[k] = insum_3; //fqEw1_1

	insum_4 = 0;
	for (l = 0; l < R.size(); l++)
	  insum_4 += L[k].ndotGradGreens_inter[l] * R[l].w * R[l].a;
	inner_rw3[k] = insum_4; //fqEw2_1

	insum_5 = 0;
	for (l = 0; l < L.size(); l++)
	  insum_5 += L[k].Greens_bar_1[l] * L[l].w * L[l].a;
	inner_lw_11[k] = insum_5; //fEqEw1_11: G_barbar 11 (r , s1) = n1 * Grad(r, s1') * G_bar 1(s1', s1)

	insum_6 = 0;
	for (l = 0; l < L.size(); l++)
	  insum_6 += L[k].Greens_bar_2[l] * L[l].w * L[l].a;
	inner_lw_12[k] = insum_6; //fEqEw1_12: G_barbar 12 (r , s1) = n1 * Grad(r, s1') * G_bar 2(s1', s1)

	insum_7 = 0;
	for (l = 0; l < R.size(); l++)
	  insum_7 += L[k].Greens_bar_1_inter[l] * R[l].w * R[l].a;
	inner_rw_11[k] = insum_7; //fEqEw2_11: G_barbar 11 (r , s2) = n1 * Grad(r, s1') * G_bar 1(s1', s2)

	insum_8 = 0;
	for (l = 0; l < R.size(); l++)
	  insum_8 += L[k].Greens_bar_2_inter[l] * R[l].w * R[l].a;
	inner_rw_12[k] = insum_8; //fEqEw2_12: G_barbar 12 (r , s2) = n1 * Grad(r, s1') * G_bar 2(s1', s2)
      }

  #pragma omp for schedule(dynamic) 
    for (k = 0; k < R.size(); k++)
      {
        insum_1 = 0;
        for (l = 0; l < L.size(); l++)
          insum_1 += R[k].Greens_inter[l] * L[l].w * L[l].a;
        inner_lw2[k] = insum_1; //fw1Eq_2

        insum_2 = 0;
        for (l = 0; l < R.size(); l++)
          insum_2 += R[k].Greens[l] * R[l].w * R[l].a;
        inner_rw2[k] = insum_2; //fw2Eq_2

	insum_3 = 0;
	for (l = 0; l < L.size(); l++)
	  insum_3 += R[k].ndotGradGreens_inter[l] * L[l].w * L[l].a;
	inner_lw4[k] = insum_3; //fqEw1_2

	insum_4 = 0;
	for (l = 0; l < R.size(); l++)
	  insum_4 += R[k].ndotGradGreens[l] * R[l].w * R[l].a;
	inner_rw4[k] = insum_4; //fqEw2_2

        insum_5 = 0;
        for (l = 0; l < L.size(); l++)
          insum_5 += R[k].Greens_bar_2_inter[l] * L[l].w * L[l].a;
        inner_lw_22[k] = insum_5; //fEqEw1_22: G_barbar 22 (r , s1) = n2 * Grad(r, s2') * G_bar 2(s2', s1)

        insum_6 = 0;
        for (l = 0; l < L.size(); l++)
          insum_6 += R[k].Greens_bar_1_inter[l] * L[l].w * L[l].a;
        inner_lw_21[k] = insum_6; //fEqEw1_21: G_barbar 21 (r , s1) = n2 * Grad(r, s2') * G_bar 1(s2', s1)

        insum_7 = 0;
        for (l = 0; l < R.size(); l++)
          insum_7 += R[k].Greens_bar_2[l] * R[l].w * R[l].a;
        inner_rw_22[k] = insum_7; //fEqEw2_22: G_barbar 22 (r , s2) = n2 * Grad(r, s2') * G_bar 2(s2', s2)

        insum_8 = 0;
        for (l = 0; l < R.size(); l++)
          insum_8 += R[k].Greens_bar_1[l] * R[l].w * R[l].a;
        inner_rw_21[k] = insum_8; //fEqEw2_21: G_barbar 22 (r , s2) = n2 * Grad(r, s2') * G_bar 1(s2', s2)
      }

  //calculate free charge - free charge interaction
  #pragma omp for schedule(dynamic) nowait
    for (i = 0; i < ion.size(); i++)
    {
      fqq = 0.0;
      for (j = 0; j < ion.size(); j++)
      {
        vec = ion[i].posvec-ion[j].posvec;
	vec = Mini_Image(vec, box.lx);

        //------charged sheet method for pbc
     	dz = vec.z;
     	r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
     	r_2 = sqrt(0.25 + (dz/box.lx)*(dz/box.lx));
     	fcsh_z = 4 * box.lx * log((0.5 + r_1)/r_2) - fabs(dz) * (2 * pi - 4 * atan(4*fabs(dz)*r_1/box.lx));
     	fcsh_inf = -2 * pi * fabs(dz);
	//-------
     	fqq += ion[i].q * (ion[j].q/ (box.lx * box.lx)) * 0.5 * (1/ion[i].epsilon + 1/ion[j].epsilon) * (fcsh_inf - fcsh_z);

        if (i == j) continue;


	fqq += 0.5 * ion[i].q * ion[j].q * Greens(vec)/ ion[i].epsilon;

      }

      fqEq_1=0;
      fEqEq_11 = 0;
      fEqEq_12 = 0;
      for (k = 0; k < L.size(); k++)
      {
        vec = L[k].posvec-ion[i].posvec;
	vec = Mini_Image(vec, box.lx);

        fqEq_1 += Greens(vec) * 0.5 * box.ed_left * inner_ql[k];
        fEqEq_11 += (L[k].normalvec * Grad_Greens(vec)) * 0.5 * box.ed_left * box.ed_left * inner_qll[k]  * L[k].a;
        fEqEq_12 += (L[k].normalvec * Grad_Greens(vec)) * 0.5 * 2 * box.ed_left * box.ed_right * inner_qlr[k] * L[k].a;
      }
      fqEq_1 = fqEq_1 * ion[i].q / ion[i].epsilon;
      fEqEq_11 = fEqEq_11 * ion[i].q / ion[i].epsilon;
      fEqEq_12 = fEqEq_12 * ion[i].q / ion[i].epsilon;


      fqEq_2 = 0;
      fEqEq_22 =0;
      for (k= 0; k < R.size(); k++)
      {
        vec = R[k].posvec-ion[i].posvec;
	vec = Mini_Image(vec, box.lx);

        fqEq_2 += Greens(vec) * 0.5 * box.ed_right *  inner_qr[k] ;
        fEqEq_22 += (R[k].normalvec * Grad_Greens(vec)) * 0.5 * box.ed_right  * box.ed_right *   inner_qrr[k]*R[k].a;
      }
      fqEq_2 = fqEq_2 * ion[i].q / ion[i].epsilon;
      fEqEq_22 = fEqEq_22 * ion[i].q / ion[i].epsilon;


      ion_energy[i] = fqq + fqEq_1 + fqEq_2+  fEqEq_11 + fEqEq_22 +  fEqEq_12 ;
      }

  //calculate free charge - free charge/induced charge interaction
  #pragma omp for schedule(dynamic) nowait
    for (i = 0; i < ion.size(); i++)
       {
      fw1q=0;					//G
      fw1Eq_1= 0;				//G bar for R_rho omega1 (left interface)
      fqEw1_1 = 0;
      fw2Eq_1 =0;				//G bar for R_rho omega2 (right interface)
      fqEw2_1 = 0;
      fEqEw1_11= 0;				//G barbar for R_rho omega1 (left interface)
      fEqEw1_12= 0;
      fEqEw2_11= 0;
      fEqEw2_12= 0;
      for (k = 0 ; k < L.size() ; k++)
      {
 	vec = L[k].posvec - ion[i].posvec;
	vec = Mini_Image(vec, box.lx);

	greens = Greens(vec) * L[k].a * ion[i].q / ion[i].epsilon;
	ndot_grad_greens = (L[k].normalvec * Grad_Greens(vec)) * L[k].a * ion[i].q / ion[i].epsilon;
	
	//------charged sheet method for pbc
        dz = vec.z;
        r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
        r_2 = sqrt(0.25 + (dz/box.lx)*(dz/box.lx));
        fcsh_z = 4 * box.lx * log((0.5 + r_1)/r_2) - fabs(dz) * (2 * pi - 4 * atan(4*fabs(dz)*r_1/box.lx));
        fcsh_inf = -2 * pi * fabs(dz);

        fw1q += ion[i].q * (L[k].w * L[k].a/ (box.lx * box.lx)) * (1 - box.em_left/ion[i].epsilon) * (fcsh_inf - fcsh_z);
	//----------
        fw1q += 0.5 * L[k].w * (1 - box.em_left/ion[i].epsilon) * greens * ion[i].epsilon;
        fw1Eq_1 += ndot_grad_greens * (- 0.5) *  box.ed_left * (2 * box.em_left -1) * inner_lw1[k];
        fqEw1_1 += greens * 0.5 * box.ed_left * inner_lw3[k] * L[k].a;
        fw2Eq_1 += ndot_grad_greens * (- 0.5) *  box.ed_left * (2 * box.em_right -1) * inner_rw1[k];
        fqEw2_1 += greens * 0.5 * box.ed_left * inner_rw3[k];
        fEqEw1_11 += ndot_grad_greens * 0.5 * 2 *  box.ed_left * box.ed_left * inner_lw_11[k];
        fEqEw1_12 += ndot_grad_greens * 0.5 * 2 *  box.ed_left * box.ed_right * inner_lw_12[k];
        fEqEw2_11 += ndot_grad_greens * 0.5 * 2 *  box.ed_left * box.ed_left * inner_rw_11[k];
        fEqEw2_12 += ndot_grad_greens * 0.5 * 2 *  box.ed_left * box.ed_right * inner_rw_12[k];
      }

      fw2q=0;
      fw1Eq_2 = 0;
      fqEw1_2 = 0;
      fw2Eq_2 =0;
      fqEw2_2 = 0;
      fEqEw1_22= 0;
      fEqEw1_21= 0;
      fEqEw2_22= 0;				//G barbar for R_rho omega2 (right interface)
      fEqEw2_21= 0;
      for (k = 0 ; k < R.size() ; k++)
      {
        vec = R[k].posvec - ion[i].posvec;
	vec = Mini_Image(vec, box.lx);
 
        greens = Greens(vec) * R[k].a * ion[i].q / ion[i].epsilon;
        ndot_grad_greens = (R[k].normalvec * Grad_Greens(vec)) * R[k].a * ion[i].q / ion[i].epsilon;

        //------charged sheet method for pbc
        dz = vec.z;
        r_1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
        r_2 = sqrt(0.25 + (dz/box.lx)*(dz/box.lx));
        fcsh_z = 4 * box.lx * log((0.5 + r_1)/r_2) - fabs(dz) * (2 * pi - 4 * atan(4*fabs(dz)*r_1/box.lx));
        fcsh_inf = -2 * pi * fabs(dz);

        fw1q += ion[i].q * (R[k].w * R[k].a/ (box.lx * box.lx)) * (1 - box.em_right/ion[i].epsilon) * (fcsh_inf - fcsh_z);
	//---------
        fw2q += 0.5 * R[k].w * (1 - box.em_right/ion[i].epsilon) * greens * ion[i].epsilon;
        fw1Eq_2 += ndot_grad_greens * (- 0.5) *  box.ed_right * (2 * box.em_left -1) * inner_lw2[k];
        fqEw1_2 += greens * 0.5 * box.ed_right * inner_lw4[k];
        fw2Eq_2 += ndot_grad_greens * (- 0.5) *  box.ed_right * (2 * box.em_right -1) * inner_rw2[k];
        fqEw2_2 += greens * 0.5 * box.ed_right * inner_rw4[k];
        fEqEw1_22 += ndot_grad_greens * 0.5 * 2 *  box.ed_right * box.ed_right * inner_lw_22[k];
        fEqEw1_21 += ndot_grad_greens * 0.5 * 2 *  box.ed_left * box.ed_right * inner_lw_21[k];
        fEqEw2_22 += ndot_grad_greens * 0.5 * 2 *  box.ed_right * box.ed_right * inner_rw_22[k];
        fEqEw2_21 += ndot_grad_greens * 0.5 * 2 *  box.ed_left * box.ed_right * inner_rw_21[k];
      } 

      ion_ind_energy[i] = fw1q + fw2q + fw1Eq_1 + fw1Eq_2 + fqEw1_1 + fqEw1_2 + fw2Eq_1 + fw2Eq_2 + fqEw2_1 + fqEw2_2 + fEqEw1_11 + fEqEw1_22 + fEqEw1_12 + fEqEw1_21 + fEqEw2_11 + fEqEw2_22 + fEqEw2_12 + fEqEw2_21;
      }



  // Excluded volume interaction energy given by purely repulsive LJ

  // ion-ion
  #pragma omp for schedule(dynamic) nowait
  for (i = 0; i < ion.size(); i++)
  {
    uljcc = 0.0;
    for (j = 0; j < ion.size(); j++)
    {
      if (j == i) continue;
      vec = ion[i].posvec - ion[j].posvec;
      vec = Mini_Image(vec, box.lx);

      r = vec.GetMagnitude();
      d = 0.5 * (ion[i].diameter + ion[j].diameter);
      elj = 1.0;
      if (r < dcut * d)
      {
        r2 = r * r;
        r6 = r2 * r2 * r2;
        d2 = d * d;
        d6 = d2 * d2 * d2;
        uljcc = uljcc +  4 * elj * (d6 / r6) * ( ( d6 / r6 ) - 1 ) + elj;
      }
      else
        uljcc = uljcc + 0.0;
    }
    lj_ion_ion[i] = uljcc;
  }

  // ion-box

  // left wall

  // ion interacting with left wall directly (self, closest)
  #pragma omp for schedule(dynamic) nowait
  for (i = 0; i < ion.size(); i++)
  {
    ulj = 0;
    if(ion[i].posvec.z < -0.5*box.lz + ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk. replacing 1 by ion[i].diameter -Yufei -Vikram -Vikram
    {
      dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_left,VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, -0.5*box.lz - 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
      vec = ion[i].posvec - dummy.posvec;
      r2 = vec.GetMagnitudeSquared();
      d = 0.5 * (ion[i].diameter + dummy.diameter);
      d2 = d * d;
      elj = 1.0;
      if (r2 < dcut2 * d2)
      {
        r6 = r2 * r2 * r2;
        d6 = d2 * d2 * d2;
        ulj = 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
      }
    }
    lj_ion_leftdummy[i] = ulj;

  // ion interacting with discretized left wall

    ulj_wall = 0.0;
    if(ion[i].posvec.z < -0.5*box.lz + ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk. -Yufei -Vikram
    {
      for (k = 0; k < L.size(); k++)
      {
        wall_dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_left, VECTOR3D(L[k].posvec.x, L[k].posvec.y, L[k].posvec.z - 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
        vec = ion[i].posvec - wall_dummy.posvec;
        vec = Mini_Image(vec,box.lx);

	r2 = vec.GetMagnitudeSquared();
	d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
	d2 = d * d;
	elj = 1.0;
	if (r2 < dcut2 * d2)
	{
	  r6 = r2 * r2 * r2;
	  d6 = d2 * d2 * d2;
	  ulj_wall += 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
	}
      }
    }
    lj_ion_leftwall[i] = ulj_wall;
  }

  // right wall

  // ion interacting with right wall directly (self, closest)
  #pragma omp for schedule(dynamic) nowait
  for (i = 0; i < ion.size(); i++)
  {
    ulj = 0;
    if(ion[i].posvec.z > 0.5*box.lz - ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
    {
      dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_right,VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, 0.5*box.lz + 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
      vec = ion[i].posvec - dummy.posvec;
      r2 = vec.GetMagnitudeSquared();
      d = 0.5 * (ion[i].diameter + dummy.diameter);
      d2 = d * d;
      elj = 1.0;
      if (r2 < dcut2 * d2)
      {
        r6 = r2 * r2 * r2;
        d6 = d2 * d2 * d2;
        ulj = 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
      }
    }
    lj_ion_rightdummy[i] = ulj;

  // ion interacting with discretized right wall

    ulj_wall = 0;
    if(ion[i].posvec.z > 0.5*box.lz - ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
    {
      for (k = 0; k < R.size(); k++)
      {
        wall_dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.e_right, VECTOR3D(R[k].posvec.x, R[k].posvec.y, R[k].posvec.z + 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
        vec = ion[i].posvec - wall_dummy.posvec;
        vec = Mini_Image(vec, box.lx);

	r2 = vec.GetMagnitudeSquared();
	d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
	d2 = d * d;
	elj = 1.0;
	if (r2 < dcut2 * d2)
	{
	  r6 = r2 * r2 * r2;
	  d6 = d2 * d2 * d2;
	  ulj_wall += 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
	}
      }
    }
    lj_ion_rightwall[i] = ulj_wall;
  }
}
  ind_ind = 0;
  for (k = 0 ; k < L.size(); k ++)
    ind_ind = ind_ind +ind_energy_left[k];
  for (k = 0 ; k < R.size(); k ++)
    ind_ind = ind_ind +  ind_energy_right[k];

  ion_ion = 0;
  for (unsigned int i=0; i< ion.size(); i++)
    ion_ion = ion_ion + ion_energy[i]; //ion ion electrostatic energy R_rhorho

  ion_ind = 0;
  for (unsigned int i =0; i < ion.size(); i++)
    ion_ind = ion_ind + ion_ind_energy[i]; //ion - induced charge electrostatic interactions R_rho omega

  // electrostatic potential energy
  double coulomb = (ind_ind + ion_ion + ion_ind) * scalefactor;

  double total_lj_ion_ion = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_ion += lj_ion_ion[i];
  total_lj_ion_ion = 0.5 * total_lj_ion_ion;
  // factor of half for double counting, same reasoning as electrostatic energy

  double total_lj_ion_leftdummy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_leftdummy += lj_ion_leftdummy[i];

  double total_lj_ion_leftwall = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_leftwall += lj_ion_leftwall[i];

  double total_lj_ion_rightdummy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_rightdummy += lj_ion_rightdummy[i];

  double total_lj_ion_rightwall = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_rightwall += lj_ion_rightwall[i];

  double potential = coulomb + total_lj_ion_ion + total_lj_ion_leftdummy + total_lj_ion_leftwall
                            + total_lj_ion_rightdummy + total_lj_ion_rightwall;



  inner_lw1.clear();
  inner_lw2.clear();
  inner_lw3.clear();
  inner_lw4.clear();
  inner_rw1.clear();
  inner_rw2.clear();
  inner_rw3.clear();
  inner_rw4.clear();
  inner_lw_11.clear();
  inner_lw_12.clear();
  inner_lw_22.clear();
  inner_lw_21.clear();

  inner_rw_11.clear();
  inner_rw_12.clear();
  inner_rw_22.clear();
  inner_rw_21.clear();

  return potential;


}
