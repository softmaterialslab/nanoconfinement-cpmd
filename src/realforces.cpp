// This file contains the routine that computes the force
// on the induced charge at vertex k and the force on the particle i
// for all k and i

#include "forces.h"

// Total Force on all degrees of freedom
void compute_force_real_dof(vector<PARTICLE> &ion, INTERFACE &box, char flag) {

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

//    cout <<"L.size()" << L.size() << endl;
//    cout <<"R.size()" << R.size() << endl;
//    cout <<"ion.size()" << ion.size() << endl;
//    cout <<"sizFVecMesh" << sizFVecMesh << endl;
//    cout <<"sizFVecIons" << sizFVecIons << endl;
//    cout <<"lowerBoundMesh" << lowerBoundMesh << endl;
//    cout <<"upperBoundMesh" << upperBoundMesh << endl;
//    cout <<"lowerBoundIons" << lowerBoundIons << endl;
//    cout <<"upperBoundIons" << upperBoundIons << endl;

    unsigned int i, j, k, l;
    int factor;
    double hcsh, E_z, r_1, r_2, dz, gcsh_z, gcsh_inf;
    double insum1, insum2, innersum_1, innersum_2, innersum_3, innersum_4;
    // VECTOR3D insum_1, insum_2, insum_3, insum_4,insum_5, insum_6, insum_7, insum_8;
    VECTOR3D vec, grad_G, grad_ndot_G;
    double left_gw1q, left_gw1q_1_gw1w1, left_gw1q_2_gw1w2;
    double right_gw2q, right_gw2q_1_gw2w1, right_gw2q_2_gw2w2;
    VECTOR3D hqq_ij, hqq_1_ij, hqq_2_ij, hqq_11_ij, hqq_22_ij, hqq_12_ij, hqq_21_ij;
    //VECTOR3D hqq_ji, hqq_1_ji, hqq_2_ji, hqq_11_ji, hqq_22_ji, hqq_12_ji, hqq_21_ji;
    vector<double> inner_qq_11(sizFVecMesh, 0.0);            //precaculation of hEqEq_11
    vector<double> inner_qq_12(sizFVecMesh, 0.0);            //precaculation of hEqEq_12
    vector<double> inner_qq_21(sizFVecMesh, 0.0);        //precaculation of hEqEq_21
    vector<double> inner_qq_22(sizFVecMesh, 0.0);            //precaculation of hEqEq_22

    //VECTOR3D hqw1, hqEw1_1, hqEw1_2, hw1Eq_1, hw1Eq_2,hEqEw1_11, hEqEw1_22, hEqEw1_12, hEw1Eq_12;
    VECTOR3D hqw1, hqEw1_1, hw1Eq_1;
    vector<double> inner_qw1_1(sizFVecMesh, 0.0);            //precaculation of hqEw1_1
    vector<double> inner_w1q_1(sizFVecMesh, 0.0);            //precaculation of hw1Eq_1

    //VECTOR3D hqw2, hqEw2_1, hqEw2_2, hw2Eq_1, hw2Eq_2,hEqEw2_11, hEqEw2_22, hEw2Eq_21, hEqEw2_21;
    VECTOR3D hqw2, hqEw2_2, hw2Eq_2;
    vector<double> inner_qw2_2(sizFVecMesh, 0.0);            //precaculation of hqEw2_2
    vector<double> inner_w2q_2(sizFVecMesh, 0.0);            //precaculation of hw2Eq_2


    vector<double> inner_lw1(sizFVecMesh, 0.0);          //precaculation og gw1Eq_1
    vector<double> inner_lw2(sizFVecMesh, 0.0);             //precaculation og gw1Eq_2
    vector<double> inner_lw3(sizFVecMesh, 0.0);             //precaculation og gEqw1_1
    vector<double> inner_lw4(sizFVecMesh, 0.0);             //precaculation og gEqw1_2

    vector<VECTOR3D> lj_ion_ion(ion.size(), VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj_ion_leftdummy(ion.size(), VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj_ion_left_wall(ion.size(), VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj_ion_rightdummy(ion.size(), VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj_ion_right_wall(ion.size(), VECTOR3D(0, 0, 0));
    VECTOR3D fljcc, flj;
    PARTICLE dummy, wall_dummy;
    //double r2, d, d2, elj, r6, r12,d6,d12;


    ///////////// MPI Message objects
    vector<double> inner_lw1Gather(L.size(), 0.0);
    vector<double> inner_lw2Gather(R.size(), 0.0);
    vector<double> inner_lw3Gather(L.size(), 0.0);
    vector<double> inner_lw4Gather(R.size(), 0.0);

    vector<double> inner_qq_11Gather(L.size(), 0.0);            //precaculation of hEqEq_11
    vector<double> inner_qq_12Gather(L.size(), 0.0);            //precaculation of hEqEq_12
    vector<double> inner_qq_21Gather(R.size(), 0.0);        //precaculation of hEqEq_21
    vector<double> inner_qq_22Gather(R.size(), 0.0);            //precaculation of hEqEq_22


    vector<double> inner_qw1_1Gather(L.size(), 0.0);            //precaculation of hqEw1_1
    vector<double> inner_w1q_1Gather(L.size(), 0.0);            //precaculation of hw1Eq_1
    vector<double> fwLDist(sizFVecMesh, 0.0);
    vector<double> fwLDistGather(L.size(), 0.0);


    vector<double> inner_qw2_2Gather(R.size(), 0.0);            //precaculation of hqEw2_2
    vector<double> inner_w2q_2Gather(R.size(), 0.0);            //precaculation of hw2Eq_2
    vector<double> fwRDist(sizFVecMesh, 0.0);
    vector<double> fwRDistGather(R.size(), 0.0);

    vector<VECTOR3D> forceVec(sizFVecIons, 0.0);
    vector<VECTOR3D> forceVecGather(ion.size(), 0.0);

#pragma omp for schedule(dynamic) private(k, i, insum1, insum2, vec)
    for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
        insum1 = 0;
        insum2 = 0;

        for (i = 0; i < ion.size(); i++) {
            vec = L[k].posvec - ion[i].posvec;
            vec = Mini_Image(vec, box.lx);

            insum1 += (L[k].normalvec * Grad_Greens(vec)) * ion[i].q / ion[i].epsilon;
            insum2 += Greens(vec) * ion[i].q / ion[i].epsilon;
        }

        inner_lw1[k - lowerBoundMesh] = insum1;        //gw1Eq_1				inner_lw1=inner_rw1
        inner_lw3[k - lowerBoundMesh] = insum2;        //gqEw1_1				inner_lw1=inner_rw1
    }

#pragma omp for schedule(dynamic) private(k, i, insum1, insum2, vec)
    for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
        insum1 = 0;
        insum2 = 0;

        for (i = 0; i < ion.size(); i++) {
            vec = R[k].posvec - ion[i].posvec;
            vec = Mini_Image(vec, box.lx);

            insum1 += (R[k].normalvec * Grad_Greens(vec)) * ion[i].q / ion[i].epsilon;
            insum2 += Greens(vec) * ion[i].q / ion[i].epsilon;
        }

        inner_lw2[k - lowerBoundMesh] = insum1;        //gw1Eq_1				inner_lw1=inner_rw1
        inner_lw4[k - lowerBoundMesh] = insum2;        //gqEw1_1				inner_lw1=inner_rw1
    }

//inner_lw1,inner_lw3,inner_lw2, and inner_lw4 broadcasting using all gather = gather + broadcast
    if (world.size() > 1) {

        //cout <<"This is proc : " << world.rank() << ", LowerBound"<<lowerBoundMesh<< ", UpperBound"<<upperBoundMesh<<endl;
        //cout  << world.rank() <<" : Size of the data trying to send : " << saveinsum.size() <<"saveinsumGather size : " << saveinsumGather.size() <<endl;
        all_gather(world, &inner_lw1[0], inner_lw1.size(), inner_lw1Gather);
        all_gather(world, &inner_lw3[0], inner_lw3.size(), inner_lw3Gather);
        all_gather(world, &inner_lw2[0], inner_lw2.size(), inner_lw2Gather);
        all_gather(world, &inner_lw4[0], inner_lw4.size(), inner_lw4Gather);

        //cout <<"Done proc : " << world.rank() <<endl;


    } else {
        for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
            inner_lw1Gather[k] = inner_lw1[k - lowerBoundMesh];
            inner_lw3Gather[k] = inner_lw3[k - lowerBoundMesh];
            inner_lw2Gather[k] = inner_lw2[k - lowerBoundMesh];
            inner_lw4Gather[k] = inner_lw4[k - lowerBoundMesh];
        }
    }

#pragma omp for schedule(dynamic) private(k, i, l, vec, dz, r_1, r_2, gcsh_z, gcsh_inf, \
     left_gw1q, left_gw1q_1_gw1w1, left_gw1q_2_gw1w2, \
  innersum_1, innersum_2, innersum_3, innersum_4)
    for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
        innersum_1 = 0;
        for (l = 0; l < L.size(); l++)
            innersum_1 += L[k].Greens[l] * inner_lw1Gather[l] * L[l].a;
        inner_qq_11[k - lowerBoundMesh] = innersum_1;//precaculation of hEqEq_11

        innersum_2 = 0;
        for (l = 0; l < R.size(); l++)
            innersum_2 += L[k].Greens_inter[l] * inner_lw2Gather[l] * R[l].a;
        inner_qq_12[k - lowerBoundMesh] = innersum_2;//precaculation of hEqEq_12


        //pre calculation of ion - leftplane interaction
        innersum_1 = 0;
        innersum_2 = 0;
        for (l = 0; l < L.size(); l++) {

            innersum_1 += L[k].ndotGradGreens[l] * L[l].w * L[l].a;
            innersum_2 += L[k].inner_wiq_i[l] * L[l].w;

        }

        //pre calculation of ion - rightplane interaction
        innersum_3 = 0;
        innersum_4 = 0;
        for (l = 0; l < R.size(); l++) {

            innersum_3 += L[k].ndotGradGreens_inter[l] * R[l].w * R[l].a;
            innersum_4 += L[k].inner_wjq_i[l] * R[l].w;

        }
        inner_qw1_1[k - lowerBoundMesh] = (innersum_1 + innersum_3) * (-0.5 * box.ed_left); //hqEw1_1
        inner_w1q_1[k - lowerBoundMesh] = (innersum_2 + innersum_4); //hqEw1_1

        //calculate force on vertexes on left plane

        //force from free ions
        left_gw1q = 0;                    //Gwq
        for (i = 0; i < ion.size(); i++) {
            vec = L[k].posvec - ion[i].posvec;
            vec = Mini_Image(vec, box.lx);

            //---------charged sheet method for pbc
            dz = vec.z;
            r_1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
            r_2 = (0.25 + (dz / box.lx) * (dz / box.lx));
            gcsh_z = 2 * box.lx * log((0.5 + r_1) * (0.5 + r_1) / r_2) -
                     fabs(dz) * (2 * pi - 4 * atan(4 * fabs(dz) * r_1 / box.lx));
            gcsh_inf = -2 * pi * fabs(dz);

            left_gw1q += (ion[i].q / (box.lx * box.lx)) * (1 - box.em_left / ion[i].epsilon) * (gcsh_inf - gcsh_z);
            //--------------
            left_gw1q += (-0.5) * ion[i].q * (1 - box.em_left / ion[i].epsilon) * Greens(vec);
        }
        //force from interfaces
        left_gw1q_1_gw1w1 = 0;
        for (l = 0; l < L.size(); l++)
            left_gw1q_1_gw1w1 +=
                    L[k].fwq_1[l] * inner_lw1Gather[l] + L[k].fqw_1[l] * inner_lw3Gather[l] + L[k].fwiwi[l] * L[l].w;

        left_gw1q_2_gw1w2 = 0;
        for (l = 0; l < R.size(); l++)
            left_gw1q_2_gw1w2 +=
                    L[k].fwq_2[l] * inner_lw2Gather[l] + L[k].fqw_2[l] * inner_lw4Gather[l] + L[k].fwiwj[l] * R[l].w;

        //L[k].fw = left_gw1q + left_gw1q_1_gw1w1 + left_gw1q_2_gw1w2;
        fwLDist[k - lowerBoundMesh] = left_gw1q + left_gw1q_1_gw1w1 + left_gw1q_2_gw1w2;
    }


#pragma omp for schedule(dynamic) private(k, i, l, vec, dz, r_1, r_2, gcsh_z, gcsh_inf, \
  right_gw2q, right_gw2q_1_gw2w1, right_gw2q_2_gw2w2, \
  innersum_1, innersum_2, innersum_3, innersum_4)
    for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
        innersum_1 = 0;
        for (l = 0; l < L.size(); l++)
            innersum_1 += R[k].Greens_inter[l] * inner_lw1Gather[l] * L[l].a;
        inner_qq_21[k - lowerBoundMesh] = innersum_1;//precaculation of hEqEq_21

        innersum_2 = 0;
        for (l = 0; l < R.size(); l++)
            innersum_2 += R[k].Greens[l] * inner_lw2Gather[l] * R[l].a;
        inner_qq_22[k - lowerBoundMesh] = innersum_2;//precaculation of hEqEq_22


        //pre calculation of ion - leftplane interaction
        innersum_1 = 0;
        innersum_2 = 0;
        for (l = 0; l < L.size(); l++) {
            innersum_1 += R[k].ndotGradGreens_inter[l] * L[l].w * L[l].a;
            innersum_2 += R[k].inner_wjq_i[l] * L[l].w;
        }

        //pre calculation of ion - rightplane interaction
        innersum_3 = 0;
        innersum_4 = 0;
        for (l = 0; l < R.size(); l++) {
            innersum_3 += R[k].ndotGradGreens[l] * R[l].w * R[l].a;
            innersum_4 += R[k].inner_wiq_i[l] * R[l].w;
        }
        inner_qw2_2[k - lowerBoundMesh] = (innersum_1 + innersum_3) * (-0.5 * box.ed_right);//hqEw2_2
        inner_w2q_2[k - lowerBoundMesh] = (innersum_2 + innersum_4);//hqEw2_2

        //calculate force on vertices on right plane

        //force from free ions
        right_gw2q = 0;            //Gwq
        for (i = 0; i < ion.size(); i++) {
            vec = R[k].posvec - ion[i].posvec;
            vec = Mini_Image(vec, box.lx);

            //---------charged sheet method for pbc
            dz = vec.z;
            r_1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
            r_2 = (0.25 + (dz / box.lx) * (dz / box.lx));
            gcsh_z = 2 * box.lx * log((0.5 + r_1) * (0.5 + r_1) / r_2) -
                     fabs(dz) * (2 * pi - 4 * atan(4 * fabs(dz) * r_1 / box.lx));
            gcsh_inf = -2 * pi * fabs(dz);

            right_gw2q +=
                    (ion[i].q / (box.lx * box.lx)) * (1 - box.em_right / ion[i].epsilon) * (gcsh_inf - gcsh_z);
            //-----------------
            right_gw2q += (-0.5) * ion[i].q * (1 - box.em_right / ion[i].epsilon) * Greens(vec);
        }

        //force from interfaces
        right_gw2q_1_gw2w1 = 0;
        for (l = 0; l < L.size(); l++)
            right_gw2q_1_gw2w1 +=
                    R[k].fwq_1[l] * inner_lw1Gather[l] + R[k].fqw_1[l] * inner_lw3Gather[l] + R[k].fwiwj[l] * L[l].w;

        right_gw2q_2_gw2w2 = 0;
        for (l = 0; l < R.size(); l++)
            right_gw2q_2_gw2w2 +=
                    R[k].fwq_2[l] * inner_lw2Gather[l] + R[k].fqw_2[l] * inner_lw4Gather[l] + R[k].fwiwi[l] * R[l].w;

        //R[k].fw = right_gw2q + right_gw2q_1_gw2w1 + right_gw2q_2_gw2w2;
        fwRDist[k - lowerBoundMesh] = right_gw2q + right_gw2q_1_gw2w1 + right_gw2q_2_gw2w2;
    }


    //inner_qq_11, inner_qq_12, inner_qw1_1, inner_w1q_1, fwLDist, inner_qq_21, inner_qq_22, inner_qw2_2, inner_w2q_2 and fwRDist
    // broadcasting using all gather = gather + broadcast

    if (world.size() > 1) {

        all_gather(world, &inner_qq_11[0], inner_qq_11.size(), inner_qq_11Gather);
        all_gather(world, &inner_qq_12[0], inner_qq_12.size(), inner_qq_12Gather);
        all_gather(world, &inner_qw1_1[0], inner_qw1_1.size(), inner_qw1_1Gather);
        all_gather(world, &inner_w1q_1[0], inner_w1q_1.size(), inner_w1q_1Gather);
        all_gather(world, &fwLDist[0], fwLDist.size(), fwLDistGather);
        all_gather(world, &inner_qq_21[0], inner_qq_21.size(), inner_qq_21Gather);
        all_gather(world, &inner_qq_22[0], inner_qq_22.size(), inner_qq_22Gather);
        all_gather(world, &inner_qw2_2[0], inner_qw2_2.size(), inner_qw2_2Gather);
        all_gather(world, &inner_w2q_2[0], inner_w2q_2.size(), inner_w2q_2Gather);
        all_gather(world, &fwRDist[0], fwRDist.size(), fwRDistGather);

    } else {
        for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
            inner_qq_11Gather[k] = inner_qq_11[k - lowerBoundMesh];
            inner_qq_12Gather[k] = inner_qq_12[k - lowerBoundMesh];
            inner_qw1_1Gather[k] = inner_qw1_1[k - lowerBoundMesh];
            inner_w1q_1Gather[k] = inner_w1q_1[k - lowerBoundMesh];
            fwLDistGather[k] = fwLDist[k - lowerBoundMesh];
            inner_qq_21Gather[k] = inner_qq_21[k - lowerBoundMesh];
            inner_qq_22Gather[k] = inner_qq_22[k - lowerBoundMesh];
            inner_qw2_2Gather[k] = inner_qw2_2[k - lowerBoundMesh];
            inner_w2q_2Gather[k] = inner_w2q_2[k - lowerBoundMesh];
            fwRDistGather[k] = fwRDist[k - lowerBoundMesh];
        }
    }

    //Moving updated plane forces --- MPI related update.
    for (k = 0; k < fwRDistGather.size(); k++) {

        L[k].fw = fwLDistGather[k];
        R[k].fw = fwRDistGather[k];

    }
    //// --- MPI related clear.
    fwLDist.clear();
    fwLDistGather.clear();
    fwRDist.clear();
    fwRDistGather.clear();

#pragma omp for schedule(dynamic) private(k, j, i, vec, dz, r_1, r_2, factor, hcsh, E_z, \
  grad_G, grad_ndot_G, \
  hqq_ij, hqq_1_ij, hqq_2_ij, hqq_11_ij, hqq_22_ij, hqq_12_ij, \
  hqw1, hqEw1_1, hw1Eq_1, hqw2, hqEw2_2, hw2Eq_2, dummy)
    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        hqq_ij = VECTOR3D(0, 0, 0);

        for (j = 0; j < ion.size(); j++) {
            if (i == j) continue;

            vec = ion[i].posvec - ion[j].posvec;
            vec = Mini_Image(vec, box.lx);
            //--------- charged sheet method for pbc
            dz = vec.z;
            if (dz >= 0) factor = 1;
            else factor = -1;
            r_1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
            r_2 = (0.25 + (dz / box.lx) * (dz / box.lx));
            E_z = 4 * atan(4 * fabs(dz) * r_1 / box.lx);
            hcsh = (4 / box.lx) * (1 / (r_1 * (0.5 + r_1)) - 1 / (r_2)) * dz + factor * E_z +
                   16 * fabs(dz) * (box.lx / (box.lx * box.lx + 4 * dz * dz * r_1 * r_1)) *
                   (fabs(dz) * dz / (box.lx * box.lx * r_1) + factor * r_1);

            hqq_ij.z += 2 * ion[i].q * (ion[j].q / (box.lx * box.lx)) * 0.5 *
                        (1 / ion[i].epsilon + 1 / ion[j].epsilon) * hcsh;
            //------------

            hqq_ij += Grad_Greens(vec) ^ (-0.5 * ion[i].q * ion[j].q / ion[i].epsilon); // hqq_ji = hqq_ij;
        }


        hqq_1_ij = VECTOR3D(0, 0, 0);
        hqq_11_ij = VECTOR3D(0, 0, 0);
        hqq_12_ij = VECTOR3D(0, 0, 0);

        hqw1 = VECTOR3D(0, 0, 0);
        hqEw1_1 = VECTOR3D(0, 0, 0);
        hw1Eq_1 = VECTOR3D(0, 0, 0);

        for (k = 0; k < L.size(); k++) {
            vec = ion[i].posvec - L[k].posvec;
            vec = Mini_Image(vec, box.lx);

            grad_G = Grad_Greens(vec) ^ ((ion[i].q / ion[i].epsilon) * L[k].a);
            grad_ndot_G = Grad_ndot_Greens(L[k].normalvec, vec) ^ ((ion[i].q / ion[i].epsilon) * L[k].a);

            hqq_1_ij += (grad_G ^ inner_lw1Gather[k]) + (grad_ndot_G ^ inner_lw3Gather[k]);
            hqq_11_ij += grad_ndot_G ^ inner_qq_11Gather[k];
            hqq_12_ij += grad_ndot_G ^ inner_qq_12Gather[k];

            //hqq_11_ji = hqq_11_ij;
            //hqq_12_ji = hqq_21_ij;

            //---------charged sheet method for pbc
            dz = vec.z;
            if (dz >= 0) factor = 1;
            else factor = -1;
            r_1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
            r_2 = (0.25 + (dz / box.lx) * (dz / box.lx));
            E_z = 4 * atan(4 * fabs(dz) * r_1 / box.lx);
            hcsh = (4 / box.lx) * (1 / (r_1 * (0.5 + r_1)) - 1 / (r_2)) * dz + factor * E_z +
                   16 * fabs(dz) * (box.lx / (box.lx * box.lx + 4 * dz * dz * r_1 * r_1)) *
                   (fabs(dz) * dz / (box.lx * box.lx * r_1) + factor * r_1);

            hqw1.z += -2 * (ion[i].q / ion[i].epsilon) * (L[k].w * L[k].a / (box.lx * box.lx)) * hcsh;
            //-----------

            hqw1 += grad_G ^ L[k].w;
            hqEw1_1 += grad_G ^ inner_qw1_1Gather[k];
            hw1Eq_1 += grad_ndot_G ^ inner_w1q_1Gather[k];

        }

        hqq_1_ij = hqq_1_ij ^ (-0.5 * box.ed_left);
        hqq_11_ij = hqq_11_ij ^ (-0.5 * box.ed_left * box.ed_left);
        hqq_12_ij = hqq_12_ij ^ (-0.5 * box.ed_left * box.ed_right);
        hqw1 = hqw1 ^ (-0.5 * (1 - box.em_left / ion[i].epsilon) * ion[i].epsilon);

        hqq_2_ij = VECTOR3D(0, 0, 0);
        hqq_21_ij = VECTOR3D(0, 0, 0);
        hqq_22_ij = VECTOR3D(0, 0, 0);

        hqw2 = VECTOR3D(0, 0, 0);
        hqEw2_2 = VECTOR3D(0, 0, 0);
        hw2Eq_2 = VECTOR3D(0, 0, 0);

        for (k = 0; k < R.size(); k++) {
            vec = ion[i].posvec - R[k].posvec;
            vec = Mini_Image(vec, box.lx);

            grad_G = Grad_Greens(vec) ^ ((ion[i].q / ion[i].epsilon) * R[k].a);
            grad_ndot_G = Grad_ndot_Greens(R[k].normalvec, vec) ^ ((ion[i].q / ion[i].epsilon) * R[k].a);


            hqq_2_ij += (grad_G ^ inner_lw2Gather[k]) + (grad_ndot_G ^ inner_lw4Gather[k]);
            hqq_21_ij += grad_ndot_G ^ inner_qq_21Gather[k];
            hqq_22_ij += grad_ndot_G ^ inner_qq_22Gather[k];

            //hqq_21_ji = hqq_12_ij;
            //hqq_22_ji = hqq_22_ij;

            //---------charged sheet method for pbc
            dz = vec.z;
            if (dz >= 0) factor = 1;
            else factor = -1;
            r_1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
            r_2 = (0.25 + (dz / box.lx) * (dz / box.lx));
            E_z = 4 * atan(4 * fabs(dz) * r_1 / box.lx);
            hcsh = (4 / box.lx) * (1 / (r_1 * (0.5 + r_1)) - 1 / (r_2)) * dz + factor * E_z +
                   16 * fabs(dz) * (box.lx / (box.lx * box.lx + 4 * dz * dz * r_1 * r_1)) *
                   (fabs(dz) * dz / (box.lx * box.lx * r_1) + factor * r_1);

            hqw2.z += -2 * (ion[i].q / ion[i].epsilon) * (R[k].w * R[k].a / (box.lx * box.lx)) * hcsh;
            //--------------------

            hqw2 += grad_G ^ R[k].w;
            hqEw2_2 += grad_G ^ inner_qw2_2Gather[k];
            hw2Eq_2 += grad_ndot_G ^ inner_w2q_2Gather[k];
        }

        hqq_2_ij = hqq_2_ij ^ (-0.5 * box.ed_right);
        hqq_21_ij = hqq_21_ij ^ (-0.5 * box.ed_left * box.ed_right);
        hqq_22_ij = hqq_22_ij ^ (-0.5 * box.ed_right * box.ed_right);
        hqw2 = hqw2 ^ (-0.5 * (1 - box.em_right / ion[i].epsilon) * ion[i].epsilon);

        forceVec[i - lowerBoundIons] = ((hqq_ij + hqq_1_ij + hqq_2_ij + hqq_11_ij + hqq_22_ij + hqq_12_ij + hqq_21_ij \
 + hqq_ij + hqq_11_ij + hqq_22_ij + hqq_12_ij + hqq_21_ij \
 + hqw1 + hqEw1_1 + hw1Eq_1 + hqw2 + hqEw2_2 + hw2Eq_2) ^ (scalefactor));

    }

    // excluded volume interactions given by purely repulsive LJ
    // ion-ion
#pragma omp for schedule(dynamic) private(k, j, i, vec, r_1, \
  fljcc, flj, dummy, wall_dummy)
    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        fljcc = VECTOR3D(0, 0, 0);
        for (j = 0; j < ion.size(); j++) {
            if (j == i) continue;
            vec = ion[i].posvec - ion[j].posvec;
            vec = Mini_Image(vec, box.lx);
            double r2 = vec.GetMagnitudeSquared();
            double d = 0.5 * (ion[i].diameter + ion[j].diameter);
            double d2 = d * d;
            double elj = 1.0;
            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                fljcc = fljcc + (vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            } else
                fljcc = fljcc + VECTOR3D(0, 0, 0);
        }
        forceVec[i - lowerBoundIons] = forceVec[i - lowerBoundIons] + fljcc;

        // ion-box

        // interaction with the left plane hard wall

        // make a dummy particle with the same diameter as the ion and touching left of the left wall s. t. it is closest to the ion
        flj = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.z < -0.5 * box.lz +
                              ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk. replacing 1 by diameter. -Yufei -Vikram -Vikram
        {
            dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.e_left,
                             VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, -0.5 * box.lz - 0.5 * ion[i].diameter),
                             box.lx, box.ly, box.lz);
            vec = ion[i].posvec - dummy.posvec;
            double r2 = vec.GetMagnitudeSquared();
            double d = 0.5 * (ion[i].diameter + dummy.diameter);
            double d2 = d * d;
            double elj = 1.0;
            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                flj = vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
            }
        }
        forceVec[i - lowerBoundIons] = forceVec[i - lowerBoundIons] + flj;

        // ion interacting with discretized left wall
        flj = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.z < -0.5 * box.lz +
                              ion[i].diameter)  // avoiding calculating interactions between left wall and ions in bulk. -Yufei - Vikram
        {
            for (k = 0; k < box.leftplane.size(); k++) {
                wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.e_left,
                                      VECTOR3D(box.leftplane[k].posvec.x, box.leftplane[k].posvec.y,
                                               box.leftplane[k].posvec.z - 0.5 * ion[i].diameter),
                                      box.lx, box.ly, box.lz);
                vec = ion[i].posvec - wall_dummy.posvec;
                vec = Mini_Image(vec, box.lx);
                double r2 = vec.GetMagnitudeSquared();
                double d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
                double d2 = d * d;
                double elj = 1.0;
                if (r2 < dcut2 * d2) {
                    double r6 = r2 * r2 * r2;
                    double r12 = r6 * r6;
                    double d6 = d2 * d2 * d2;
                    double d12 = d6 * d6;
                    flj = flj + (vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
                }
            }
        }
        forceVec[i - lowerBoundIons] = forceVec[i - lowerBoundIons] + flj;

        // interaction with the right plane hard wall

        // make a dummy particle with the same diameter as the ion and touching right of the right wall s. t. it is closest to the ion
        flj = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.z > 0.5 * box.lz -
                              ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
        {
            dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.e_right,
                             VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, 0.5 * box.lz + 0.5 * ion[i].diameter),
                             box.lx, box.ly, box.lz);
            vec = ion[i].posvec - dummy.posvec;
            double r2 = vec.GetMagnitudeSquared();
            double d = 0.5 * (ion[i].diameter + dummy.diameter);
            double d2 = d * d;
            double elj = 1.0;
            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                flj = vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
            }
        }
        forceVec[i - lowerBoundIons] = forceVec[i - lowerBoundIons] + flj;

        // ion interacting with discretized right wall
        flj = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.z > 0.5 * box.lz -
                              ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
        {
            for (k = 0; k < box.rightplane.size(); k++) {
                wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.e_right,
                                      VECTOR3D(box.rightplane[k].posvec.x, box.rightplane[k].posvec.y,
                                               box.rightplane[k].posvec.z + 0.5 * ion[i].diameter), box.lx, box.ly,
                                      box.lz);
                vec = ion[i].posvec - wall_dummy.posvec;
                vec = Mini_Image(vec, box.lx);
                double r2 = vec.GetMagnitudeSquared();
                double d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
                double d2 = d * d;
                double elj = 1.0;
                if (r2 < dcut2 * d2) {
                    double r6 = r2 * r2 * r2;
                    double r12 = r6 * r6;
                    double d6 = d2 * d2 * d2;
                    double d12 = d6 * d6;
                    flj = flj + (vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
                }

            }
        }
        forceVec[i - lowerBoundIons] = forceVec[i - lowerBoundIons] + flj;
    }

    //inner_qq_11, inner_qq_12, inner_qw1_1, inner_w1q_1, fwLDist, inner_qq_21, inner_qq_22, inner_qw2_2, inner_w2q_2 and fwRDist
    // broadcasting using all gather = gather + broadcast

    if (world.size() > 1) {

        all_gather(world, &forceVec[0], forceVec.size(), forceVecGather);

    } else {
        for (k = lowerBoundIons; k <= upperBoundIons; k++) {
            forceVecGather[k] = forceVec[k - lowerBoundMesh];
        }
    }

    // reassign force
    for (unsigned int k = 0; k < ion.size(); k++)
        ion[k].forvec = forceVecGather[k];

    forceVec.clear();
    forceVecGather.clear();


    // force on the fake degrees of freedom
    for (unsigned int k = 0; k < L.size(); k++)
        box.leftplane[k].fw = L[k].a * L[k].fw * scalefactor;
    // force on the fake degrees of freedom
    for (unsigned int k = 0; k < R.size(); k++)
        box.rightplane[k].fw = R[k].a * R[k].fw * scalefactor;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inner_lw1.clear();
    inner_lw2.clear();
    inner_lw3.clear();
    inner_lw4.clear();
    inner_lw1Gather.clear();
    inner_lw2Gather.clear();
    inner_lw3Gather.clear();
    inner_lw4Gather.clear();

    inner_qq_11.clear();           //precaculation of hEqEq_11
    inner_qq_11Gather.clear();
    inner_qq_12.clear();           //precaculation of hEqEq_12
    inner_qq_12Gather.clear();
    inner_qq_21.clear();           //precaculation of hEqEq_21
    inner_qq_21Gather.clear();
    inner_qq_22.clear();           //precaculation of hEqEq_22
    inner_qq_22Gather.clear();

    inner_qw1_1.clear();           //precaculation of hqEw1_1
    inner_qw1_1Gather.clear();
    inner_w1q_1.clear();           //precaculation of hw1Eq_1
    inner_w1q_1Gather.clear();

    inner_qw2_2.clear();           //precaculation of hqEw2_2
    inner_qw2_2Gather.clear();
    inner_w2q_2.clear();           //precaculation of hw2Eq_2
    inner_w2q_2Gather.clear();


    return;
}
