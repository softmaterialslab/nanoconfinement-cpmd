// This routine sets up and performs precalculations

#include "precalculations.h"

void precalculate(INTERFACE &box) {

    if (world.rank() == 0)
        cout << "precalculate Greens function on interface" << endl;

    vector<VERTEX> L = box.leftplane;
    vector<VERTEX> R = box.rightplane;
    char data_left[200], data_right[200];
    sprintf(data_left, "datafile/discre_left_%.06d.dat", int(box.leftplane.size()));
    sprintf(data_right, "datafile/discre_right_%.06d.dat", int(box.rightplane.size()));
    ifstream greens_left(data_left);
    ifstream greens_right(data_right);
    long double Greens;
    long double left_gw1w1, left_gw1Ew1_1, left_gw1Ew1_2, left_gEw1Ew1_11, left_gEw1Ew1_22, left_gEw1Ew1_12;
    long double left_gw1w2, left_gw1Ew2_1, left_gw1Ew2_2, left_gEw1Ew2_11, left_gEw1Ew2_22, left_gEw1Ew2_12;
    long double left_gw2w1, left_gw2Ew1_2, left_gw2Ew1_1, left_gEw2Ew1_22, left_gEw2Ew1_11, left_gEw2Ew1_12;
    long double right_gw2w2, right_gw2Ew2_1, right_gw2Ew2_2, right_gEw2Ew2_11, right_gEw2Ew2_22, right_gEw2Ew2_12;
    long double right_gw1w2, right_gw1Ew2_1, right_gw1Ew2_2, right_gEw1Ew2_11, right_gEw1Ew2_22, right_gEw1Ew2_12;
    long double right_gw2w1, right_gw2Ew1_2, right_gw2Ew1_1, right_gEw2Ew1_22, right_gEw2Ew1_11, right_gEw2Ew1_12;

    for (unsigned int k = 0; k < L.size(); k++) {
        for (unsigned int l = 0; l < L.size(); l++)  //Interaction within leftplane
        {

            greens_left >> Greens;
            box.leftplane[k].Greens.push_back(Greens); //Gw1w1

            greens_left >> Greens;
            box.leftplane[k].Greens_bar_1.push_back(Greens); //Gbar_1_w1w1

            greens_left >> Greens;
            box.leftplane[k].Greens_bar_2.push_back(Greens);//Gbar_2_w1w1

            greens_left >> Greens;
            box.leftplane[k].Greens_bbar_11.push_back(Greens);//Gbarbar_11_w1w1

            greens_left >> Greens;
            box.leftplane[k].Greens_bbar_22.push_back(Greens);//Gbarbar_22_w1w1

            greens_left >> Greens;
            box.leftplane[k].Greens_bbar_12.push_back(Greens);//Gbarbar_12_w1w1

            greens_left >> Greens;
            box.leftplane[k].ndotGradGreens.push_back(Greens);//Gbarbar_12_w1w1

            greens_left >> Greens;
            box.leftplane[k].Charged_sheet_corrections_intra.push_back(Greens);
        }

        for (unsigned int l = 0; l < R.size(); l++)  // Interactions between leftplane and rightplane
        {
            greens_left >> Greens;
            box.leftplane[k].Greens_inter.push_back(Greens);

            greens_left >> Greens;
            box.leftplane[k].Greens_bar_1_inter.push_back(Greens); //Gbar_1_w1w2

            greens_left >> Greens;
            box.leftplane[k].Greens_bar_2_inter.push_back(Greens);//Gbar_2_w1w2

            greens_left >> Greens;
            box.leftplane[k].Greens_bbar_11_inter.push_back(Greens);//Gbarbar_11_w1w2

            greens_left >> Greens;
            box.leftplane[k].Greens_bbar_22_inter.push_back(Greens);//Gbarbar_22_w1w2

            greens_left >> Greens;
            box.leftplane[k].Greens_bbar_12_inter.push_back(Greens);//Gbarbar_12_w1w1

            greens_left >> Greens;
            box.leftplane[k].ndotGradGreens_inter.push_back(Greens);//Gbarbar_12_w1w1

            greens_left >> Greens;
            box.leftplane[k].Charged_sheet_corrections_inter.push_back(Greens);
        }

    }


    for (unsigned int k = 0; k < R.size(); k++) {
        for (unsigned int l = 0; l < R.size(); l++) {
            greens_right >> Greens;
            box.rightplane[k].Greens.push_back(Greens);

            greens_right >> Greens;
            box.rightplane[k].Greens_bar_1.push_back(Greens); //Gbar_1_w2w2

            greens_right >> Greens;
            box.rightplane[k].Greens_bar_2.push_back(Greens);//Gbar_2_w2w2

            greens_right >> Greens;
            box.rightplane[k].Greens_bbar_11.push_back(Greens);//Gbarbar_11_w2w2

            greens_right >> Greens;
            box.rightplane[k].Greens_bbar_22.push_back(Greens);//Gbarbar_22_w2w2

            greens_right >> Greens;
            box.rightplane[k].Greens_bbar_12.push_back(Greens);//Gbarbar_12_w1w2

            greens_right >> Greens;
            box.rightplane[k].ndotGradGreens.push_back(Greens);//Gbarbar_12_w1w2

            greens_right >> Greens;
            box.rightplane[k].Charged_sheet_corrections_intra.push_back(Greens);
        }

        for (unsigned int l = 0; l < L.size(); l++) {
            greens_right >> Greens;
            box.rightplane[k].Greens_inter.push_back(Greens);

            greens_right >> Greens;
            box.rightplane[k].Greens_bar_1_inter.push_back(Greens); //Gbar_1_w2w1

            greens_right >> Greens;
            box.rightplane[k].Greens_bar_2_inter.push_back(Greens);//Gbar_2_w2w1

            greens_right >> Greens;
            box.rightplane[k].Greens_bbar_11_inter.push_back(Greens);//Gbarbar_11_w2w1

            greens_right >> Greens;
            box.rightplane[k].Greens_bbar_22_inter.push_back(Greens);//Gbarbar_22_w2w2

            greens_right >> Greens;
            box.rightplane[k].Greens_bbar_12_inter.push_back(Greens);//Gbarbar_12_w2w1

            greens_right >> Greens;
            box.rightplane[k].ndotGradGreens_inter.push_back(Greens);//Gbarbar_12_w2w1

            greens_right >> Greens;
            box.rightplane[k].Charged_sheet_corrections_inter.push_back(Greens);
        }
    }

    L = box.leftplane;
    R = box.rightplane;

    // precalculations for inner sums for fake forces calculations
    for (unsigned int k = 0; k < L.size(); k++) {
        for (unsigned int l = 0; l < L.size(); l++) {
            left_gw1w1 = (-1.0) * box.em_left * (box.em_left - 1) * L[k].Greens[l] * L[l].a;        //gw1w1
            left_gw1Ew1_1 = (-0.5) * box.ed_left * (1 - 2 * box.em_left) * L[k].Greens_bar_1[l] * L[l].a;    //gw1Ew1_1
            left_gw1Ew1_1 += (-0.5) * box.ed_left * (1 - 2 * box.em_left) * L[l].Greens_bar_1[k] * L[l].a;    //gw1Ew1_1
            left_gw1Ew1_2 = (-0.5) * box.ed_right * (1 - 2 * box.em_left) * L[k].Greens_bar_2[l] * L[l].a;    //gw1Ew1_2
            left_gw1Ew1_2 +=
                    (-0.5) * box.ed_right * (1 - 2 * box.em_left) * L[l].Greens_bar_2[k] * L[l].a;    //gw1Ew1_2
            left_gEw1Ew1_11 = (-1.0) * box.ed_left * box.ed_left * L[k].Greens_bbar_11[l] * L[l].a;    //gEw1Ew1_11
            left_gEw1Ew1_22 = (-1.0) * box.ed_right * box.ed_right * L[k].Greens_bbar_22[l] * L[l].a;    //gEw1Ew1_22
            left_gEw1Ew1_12 = (-1.0) * box.ed_left * box.ed_right * L[k].Greens_bbar_12[l] * L[l].a;    //gEw1Ew1_12;
            left_gEw1Ew1_12 += (-1.0) * box.ed_left * box.ed_right * L[l].Greens_bbar_12[k] * L[l].a;    //gEw1Ew1_12;

            //left_gw1w1 += left_gw1Ew1_1 + left_gw1Ew1_2 + left_gEw1Ew1_11 + left_gEw1Ew1_22 + left_gEw1Ew1_12;
            box.leftplane[k].fwiwi.push_back(
                    left_gw1w1 + left_gw1Ew1_1 + left_gw1Ew1_2 + left_gEw1Ew1_11 + left_gEw1Ew1_22 + left_gEw1Ew1_12);
        }

        for (unsigned int l = 0; l < R.size(); l++) {
            left_gw1w2 = (-0.5) * box.em_left * (box.em_right - 1) * L[k].Greens_inter[l] * R[l].a;            //gw1w2
            left_gw1Ew2_1 = (-0.5) * box.ed_left * (1 - 2 * box.em_left) * L[k].Greens_bar_1_inter[l] *
                            R[l].a;        //gw1Ew2_1
            left_gw1Ew2_2 = (-0.5) * box.ed_right * (1 - 2 * box.em_left) * L[k].Greens_bar_2_inter[l] *
                            R[l].a;        //gw1Ew2_2
            left_gEw1Ew2_11 =
                    (-0.5) * box.ed_left * box.ed_left * L[k].Greens_bbar_11_inter[l] * R[l].a;            //gEw1Ew2_11
            left_gEw1Ew2_22 =
                    (-0.5) * box.ed_right * box.ed_right * L[k].Greens_bbar_22_inter[l] * R[l].a;        //gEw1Ew2_22
            left_gEw1Ew2_12 =
                    (-0.5) * 2 * box.ed_left * box.ed_right * L[k].Greens_bbar_12_inter[l] * R[l].a;        //gEw1Ew2_12

            left_gw2w1 = (-0.5) * R[l].a * box.em_right * (box.em_left - 1) * R[l].Greens_inter[k];        //gw2w1
            left_gw2Ew1_2 =
                    (-0.5) * R[l].a * box.ed_right * (1 - 2 * box.em_right) * R[l].Greens_bar_2_inter[k];    //gw2Ew1_2
            left_gw2Ew1_1 =
                    (-0.5) * R[l].a * box.ed_left * (1 - 2 * box.em_right) * R[l].Greens_bar_1_inter[k];    //gw2Ew1_1
            left_gEw2Ew1_22 =
                    (-0.5) * R[l].a * box.ed_right * box.ed_right * R[l].Greens_bbar_22_inter[k];        //gEw2Ew1_22
            left_gEw2Ew1_11 =
                    (-0.5) * R[l].a * box.ed_left * box.ed_left * R[l].Greens_bbar_11_inter[k];        //gEw2Ew1_11
            left_gEw2Ew1_12 =
                    (-0.5) * R[l].a * 2 * box.ed_left * box.ed_right * R[l].Greens_bbar_12_inter[k];    //gEw2Ew1_12

            //left_gw1w2 += left_gw1Ew2_1 + left_gw1Ew2_2 + left_gEw1Ew2_11 + left_gEw1Ew2_22 + left_gEw1Ew2_12 +  left_gw2w1 + left_gw2Ew1_2 + left_gw2Ew1_1 + left_gEw2Ew1_22 + left_gEw2Ew1_11 + left_gEw2Ew1_12;
            box.leftplane[k].fwiwj.push_back(
                    left_gw1w2 + left_gw1Ew2_1 + left_gw1Ew2_2 + left_gEw1Ew2_11 + left_gEw1Ew2_22 + left_gEw1Ew2_12 +
                    left_gw2w1 + left_gw2Ew1_2 + left_gw2Ew1_1 + left_gEw2Ew1_22 + left_gEw2Ew1_11 + left_gEw2Ew1_12);
        }

    }

    for (unsigned int k = 0; k < R.size(); k++) {
        for (unsigned int l = 0; l < R.size(); l++) {
            right_gw2w2 = (-1.0) * box.em_right * (box.em_right - 1) * R[k].Greens[l] * R[l].a;            //gw2w2
            right_gw2Ew2_2 =
                    (-0.5) * box.ed_right * (1 - 2 * box.em_right) * R[k].Greens_bar_2[l] * R[l].a;    //gw2Ew2_2
            right_gw2Ew2_2 +=
                    (-0.5) * box.ed_right * (1 - 2 * box.em_right) * R[l].Greens_bar_2[k] * R[l].a;    //gw2Ew2_2
            right_gw2Ew2_1 =
                    (-0.5) * box.ed_left * (1 - 2 * box.em_right) * R[k].Greens_bar_1[l] * R[l].a;    //gw2Ew2_1
            right_gw2Ew2_1 +=
                    (-0.5) * box.ed_left * (1 - 2 * box.em_right) * R[l].Greens_bar_1[k] * R[l].a;    //gw2Ew2_1
            right_gEw2Ew2_22 =
                    (-1.0) * box.ed_right * box.ed_right * R[k].Greens_bbar_22[l] * R[l].a;        //gEw2Ew2_22
            right_gEw2Ew2_11 = (-1.0) * box.ed_left * box.ed_left * R[k].Greens_bbar_11[l] * R[l].a;        //gEw2Ew2_11
            right_gEw2Ew2_12 =
                    (-1.0) * box.ed_left * box.ed_right * R[k].Greens_bbar_12[l] * R[l].a;        //gEw2Ew2_12
            right_gEw2Ew2_12 +=
                    (-1.0) * box.ed_left * box.ed_right * R[l].Greens_bbar_12[k] * R[l].a;        //gEw2Ew2_12

            box.rightplane[k].fwiwi.push_back(
                    right_gw2w2 + right_gw2Ew2_2 + right_gw2Ew2_1 + right_gEw2Ew2_22 + right_gEw2Ew2_11 +
                    right_gEw2Ew2_12);
        }

        for (unsigned int l = 0; l < L.size(); l++) {
            right_gw2w1 = (-0.5) * box.em_right * (box.em_left - 1) * R[k].Greens_inter[l] * L[l].a;            //gw2w1
            right_gw2Ew1_2 =
                    (-0.5) * box.ed_right * (1 - 2 * box.em_right) * R[k].Greens_bar_2_inter[l] * L[l].a;    //gw2Ew1_2
            right_gw2Ew1_1 =
                    (-0.5) * box.ed_left * (1 - 2 * box.em_right) * R[k].Greens_bar_1_inter[l] * L[l].a;    //gw2Ew1_1
            right_gEw2Ew1_22 =
                    (-0.5) * box.ed_right * box.ed_right * R[k].Greens_bbar_22_inter[l] * L[l].a;        //gEw2Ew1_22
            right_gEw2Ew1_11 =
                    (-0.5) * box.ed_left * box.ed_left * R[k].Greens_bbar_11_inter[l] * L[l].a;        //gEw2Ew1_11
            right_gEw2Ew1_12 =
                    (-0.5) * 2 * box.ed_left * box.ed_right * R[k].Greens_bbar_12_inter[l] * L[l].a;        //gEw2Ew1_12

            right_gw1w2 = (-0.5) * L[l].a * box.em_left * (box.em_right - 1) * L[l].Greens_inter[k];        //gw1w2
            right_gw1Ew2_1 =
                    (-0.5) * L[l].a * box.ed_left * (1 - 2 * box.em_left) * L[l].Greens_bar_1_inter[k];    //gw1Ew2_1
            right_gw1Ew2_2 =
                    (-0.5) * L[l].a * box.ed_right * (1 - 2 * box.em_left) * L[l].Greens_bar_2_inter[k];    //gw1Ew2_2
            right_gEw1Ew2_11 =
                    (-0.5) * L[l].a * box.ed_left * box.ed_left * L[l].Greens_bbar_11_inter[k];        //gEw1Ew2_11
            right_gEw1Ew2_22 =
                    (-0.5) * L[l].a * box.ed_right * box.ed_right * L[l].Greens_bbar_22_inter[k];    //gEw1Ew2_22
            right_gEw1Ew2_12 =
                    (-0.5) * L[l].a * 2 * box.ed_left * box.ed_right * L[l].Greens_bbar_12_inter[k];    //gEw1Ew2_12

            box.rightplane[k].fwiwj.push_back(
                    right_gw2w1 + right_gw2Ew1_2 + right_gw2Ew1_1 + right_gEw2Ew1_22 + right_gEw2Ew1_11 +
                    right_gEw2Ew1_12 + right_gw1w2 + right_gw1Ew2_1 + right_gw1Ew2_2 + right_gEw1Ew2_11 +
                    right_gEw1Ew2_22 + right_gEw1Ew2_12);
        }
    }

    for (unsigned int k = 0; k < L.size(); k++) {
        for (unsigned int l = 0; l < L.size(); l++)
            box.leftplane[k].fwq_1.push_back(((-1.0) * L[k].Greens[l] * (-0.5) * box.ed_left * (2 * box.em_left - 1) \
 + (-1.0) * L[l].Greens_bar_1[k] * 0.5 * 2 * box.ed_left * box.ed_left \
 + (-1.0) * L[l].Greens_bar_2[k] * 0.5 * 2 * box.ed_left * box.ed_right) * L[l].a);
        //left_gw1Eq_1 += (-1.0) * L[k].Greens[l] * (- 0.5) *  box.ed_left * (2 * box.em_left -1) * inner_lw1[l] * L[l].a;
        //left_gEqEw1_11 += (-1.0) * L[l].Greens_bar_1[k] * 0.5 * 2 *  box.ed_left * box.ed_left * inner_lw_11[l] * L[l].a;
        //left_gEqEw1_12 += (-1.0) * L[l].Greens_bar_2[k] * 0.5 * 2 *  box.ed_left * box.ed_right * inner_lw_12[l] * L[l].a;

        for (unsigned int l = 0; l < R.size(); l++)
            box.leftplane[k].fwq_2.push_back(
                    ((-1.0) * L[k].Greens_inter[l] * (-0.5) * box.ed_right * (2 * box.em_left - 1) \
 + (-1.0) * R[l].Greens_bar_2_inter[k] * 0.5 * 2 * box.ed_right * box.ed_right\
 + (-1.0) * R[l].Greens_bar_1_inter[k] * 0.5 * 2 * box.ed_left * box.ed_right) * R[l].a);
        //left_gw1Eq_2 += (-1.0) * L[k].Greens_inter[l] * (- 0.5) *  box.ed_right * (2 * box.em_left -1) * inner_lw2[l] * R[l].a;
        //left_gEqEw1_22 += (-1.0) * R[l].Greens_bar_2_inter[k] * 0.5 * 2 *  box.ed_right * box.ed_right * inner_lw_22[l] * R[l].a;
        //left_gEqEw1_21 += (-1.0) * R[l].Greens_bar_1_inter[k] * 0.5 * 2 *  box.ed_left * box.ed_right * inner_lw_21[l] * R[l].a;

        for (unsigned int l = 0; l < L.size(); l++)
            box.leftplane[k].fqw_1.push_back((-1.0) * L[l].ndotGradGreens[k] * 0.5 * box.ed_left * L[l].a);
        //left_gqEw1_1 += (-1.0) * L[l].ndotGradGreens[k] * 0.5 * box.ed_left * inner_lw3[l] * L[l].a;


        for (unsigned int l = 0; l < R.size(); l++)
            box.leftplane[k].fqw_2.push_back((-1.0) * R[l].ndotGradGreens_inter[k] * 0.5 * box.ed_right * R[l].a);
        //left_gqEw1_2 += (-1.0) * R[l].ndotGradGreens_inter[k] * 0.5 * box.ed_right * inner_lw4[l] * R[l].a;
    }

    for (unsigned int k = 0; k < R.size(); k++) {
        for (unsigned int l = 0; l < L.size(); l++)
            box.rightplane[k].fwq_1.push_back(
                    ((-1.0) * R[k].Greens_inter[l] * (-0.5) * box.ed_left * (2 * box.em_right - 1) \
 + (-1.0) * L[l].Greens_bar_1_inter[k] * 0.5 * 2 * box.ed_left * box.ed_left \
 + (-1.0) * L[l].Greens_bar_2_inter[k] * 0.5 * 2 * box.ed_left * box.ed_right) * L[l].a);
        //right_gw2Eq_1 += ((-1.0) * R[k].Greens_inter[l] * (- 0.5) *  box.ed_left * (2 * box.em_right -1) * inner_rw1[l]  * L[l].a;
        //right_gEqEw2_11 += (-1.0) * L[l].Greens_bar_1_inter[k] * 0.5 * 2 *  box.ed_left * box.ed_left * inner_rw_11[l] * L[l].a;
        //right_gEqEw2_12 += (-1.0) * L[l].Greens_bar_2_inter[k] * 0.5 * 2 *  box.ed_left * box.ed_right * inner_rw_12[l] * L[l].a;

        for (unsigned int l = 0; l < R.size(); l++)
            box.rightplane[k].fwq_2.push_back(((-1.0) * R[k].Greens[l] * (-0.5) * box.ed_right * (2 * box.em_right - 1) \
 + (-1.0) * R[l].Greens_bar_2[k] * 0.5 * 2 * box.ed_right * box.ed_right \
 + (-1.0) * R[l].Greens_bar_1[k] * 0.5 * 2 * box.ed_left * box.ed_right) * R[l].a);
        //right_gw2Eq_2 += ((-1.0) * R[k].Greens[l] * (- 0.5) *  box.ed_right * (2 * box.em_right -1) * inner_rw2[l] * R[l].a;
        //right_gEqEw2_22 += (-1.0) * R[l].Greens_bar_2[k]  * 0.5 * 2 *  box.ed_right * box.ed_right * inner_rw_22[l] * R[l].a;
        //right_gEqEw2_21 += (-1.0) * R[l].Greens_bar_1[k] * 0.5 * 2 *  box.ed_left * box.ed_right * inner_rw_21[l] * R[l].a;

        for (unsigned int l = 0; l < L.size(); l++)
            box.rightplane[k].fqw_1.push_back((-1.0) * (L[l].ndotGradGreens_inter[k]) * 0.5 * box.ed_left * L[l].a);
        //right_gqEw2_1 +=  (-1.0) * (L[l].ndotGradGreens_inter[k]) * 0.5 * box.ed_left * inner_rw3[l] * L[l].a;

        for (unsigned int l = 0; l < R.size(); l++)
            box.rightplane[k].fqw_2.push_back((-1.0) * (R[l].ndotGradGreens[k]) * 0.5 * box.ed_right * R[l].a);
        //right_gqEw2_2 += (-1.0) * (R[l].ndotGradGreens[k]) * 0.5 * box.ed_right * inner_rw4[l] * R[l].a;

    }

    // precalculations for inner sums for real forces calculations
    for (unsigned int k = 0; k < L.size(); k++) {
        //pre calculation of ion - leftplane interaction
        for (unsigned int l = 0; l < L.size(); l++) {
//	box.leftplane[k].inner_qwi_i.push_back(( -0.5 * box.ed_left) * L[k].ndotGradGreens[l] * L[l].a);
            box.leftplane[k].inner_wiq_i.push_back(
                    L[k].Greens[l] * L[l].a * (-0.5 * box.ed_left * (1 - 2 * box.em_left)) +
                    L[k].Greens_bar_1[l] * L[l].a * (-0.5 * 2 * box.ed_left * box.ed_left) +
                    L[k].Greens_bar_2[l] * L[l].a * (-0.5 * 2 * box.ed_left * box.ed_right));
        }

        //pre calculation of ion - rightplane interaction
        for (unsigned int l = 0; l < R.size(); l++) {
//	box.leftplane[k].inner_qwj_i.push_back((-0.5 * box.ed_left) * L[k].ndotGradGreens_inter[l] * R[l].a);
            box.leftplane[k].inner_wjq_i.push_back(
                    L[k].Greens_inter[l] * R[l].a * (-0.5 * box.ed_left * (1 - 2 * box.em_right)) +
                    L[k].Greens_bar_1_inter[l] * R[l].a * (-0.5 * 2 * box.ed_left * box.ed_left) +
                    L[k].Greens_bar_2_inter[l] * R[l].a * (-0.5 * 2 * box.ed_left * box.ed_right));
        }
    }

    for (unsigned int k = 0; k < R.size(); k++) {

        //pre calculation of ion - leftplane interaction
        for (unsigned int l = 0; l < L.size(); l++) {
            //box.rightplane[k].inner_qwj_i.push_back(( -0.5 * box.ed_right) * R[k].ndotGradGreens_inter[l] * L[l].a);
            box.rightplane[k].inner_wjq_i.push_back(
                    R[k].Greens_inter[l] * L[l].a * (-0.5 * box.ed_right * (1 - 2 * box.em_left)) +
                    R[k].Greens_bar_2_inter[l] * L[l].a * (-0.5 * 2 * box.ed_right * box.ed_right) +
                    R[k].Greens_bar_1_inter[l] * L[l].a * (-0.5 * 2 * box.ed_left * box.ed_right));
        }

        //pre calculation of ion - rightplane interaction

        for (unsigned int l = 0; l < R.size(); l++) {

            //box.rightplane[k].inner_qwi_i.push_back((-0.5 * box.ed_right) * R[k].Greens[l] * R[l].a);
            box.rightplane[k].inner_wiq_i.push_back(
                    R[k].Greens[l] * R[l].a * (-0.5 * box.ed_right * (1 - 2 * box.em_right)) +
                    R[k].Greens_bar_2[l] * R[l].a * (-0.5 * 2 * box.ed_right * box.ed_right) +
                    R[k].Greens_bar_1[l] * R[l].a * (-0.5 * 2 * box.ed_left * box.ed_right));
        }
    }

    return;
}


