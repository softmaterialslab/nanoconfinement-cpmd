// This file contains the routines 

#include "functions.h"
#include "fmd.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

// overload out
ostream &operator<<(ostream &os, VECTOR3D vec) {
    os << vec.x << setw(15) << vec.y << setw(15) << vec.z;
    return os;
}

// make bins
void make_bins(vector<BIN> &bin, INTERFACE &box, double bin_width) {
    unsigned int number_of_bins = int(box.lz / bin_width);
    bin.resize(number_of_bins);
    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
        bin[bin_num].set_up(bin_num, bin_width, box.lx, box.ly, box.lz);
    if (world.rank() == 0) {
        ofstream listbin("outfiles/listbin.dat");
        for (unsigned int num = 0; num < bin.size(); num++)
            listbin << bin[num].n << setw(15) << bin[num].width << setw(15) << bin[num].volume << setw(15)
                    << bin[num].lower
                    << setw(15) << bin[num].higher << endl;
        listbin.close();
    }
    return;
}

// initialize velocities of particles to start simulation
void initialize_particle_velocities(vector<PARTICLE> &ion, vector<THERMOSTAT> &bath) {
    if (bath.size() == 1) {
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].velvec = VECTOR3D(0, 0, 0);                    // initialized velocities
        if (world.rank() == 0)
            cout << "Velocities initialized to 0" << endl;
        return;
    }
    double p_sigma = sqrt(kB * bath[0].T / (2.0 * ion[0].m));        // Maxwell distribution width

    // same random numbers used to generate the gaussian distribution every time. seed is fixed.
    // let me know if you need to change the rnd nums every run.
    UTILITY ugsl;

    for (unsigned int i = 0; i < ion.size(); i++)
        ion[i].velvec = VECTOR3D(gsl_ran_gaussian(ugsl.r, p_sigma), gsl_ran_gaussian(ugsl.r, p_sigma),
                                 gsl_ran_gaussian(ugsl.r, p_sigma));    // initialized velocities
    VECTOR3D average_velocity_vector = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < ion.size(); i++)
        average_velocity_vector = average_velocity_vector + ion[i].velvec;
    average_velocity_vector = average_velocity_vector ^ (1.0 / ion.size());
    for (unsigned int i = 0; i < ion.size(); i++)
        ion[i].velvec = ion[i].velvec - average_velocity_vector;
    return;
}

// initialize velocities of fake degrees to start simulation
void initialize_fake_velocities(vector<VERTEX> &s, vector<THERMOSTAT> &fake_bath) {
    if (fake_bath.size() == 1)                            // only the dummy is there
    {
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].vw = 0.0;
        if (world.rank() == 0)
            cout << "No thermostat for fake system" << endl;
        return;
    }
    UTILITY ugsl;
    for (unsigned int k = 0; k < s.size(); k++) {
        double w_sigma = sqrt(kB * fake_bath[0].T / (2.0 * s[0].mu));        // Maxwell distribution width
        s[k].vw = gsl_ran_gaussian(ugsl.r, w_sigma);                            // initialized velocities
    }
    double average_fake_momentum = 0.0;
    for (unsigned int k = 0; k < s.size(); k++)
        average_fake_momentum = average_fake_momentum + s[k].mu * s[k].vw;
    average_fake_momentum = average_fake_momentum * (1.0 / s.size());
    for (unsigned int k = 0; k < s.size(); k++)
        s[k].vw = s[k].vw - average_fake_momentum / (s[k].mu);
    return;
}

// make movie
void make_movie(int num, vector<PARTICLE> &ion, INTERFACE &box) {
    if (world.rank() == 0) {
        ofstream outdump("outfiles/p.lammpstrj", ios::app);
        outdump << "ITEM: TIMESTEP" << endl;
        outdump << num - 1 << endl;
        outdump << "ITEM: NUMBER OF ATOMS" << endl;
        outdump << ion.size() << endl;
        outdump << "ITEM: BOX BOUNDS" << endl;
        outdump << -0.5 * box.lx << "\t" << 0.5 * box.lx << endl;
        outdump << -0.5 * box.ly << "\t" << 0.5 * box.ly << endl;
        outdump << -0.5 * box.lz << "\t" << 0.5 * box.lz << endl;
//  outdump << -dsphere.box_radius << "\t" << dsphere.box_radius << endl;
        outdump << "ITEM: ATOMS index type x y z" << endl;
        string type;
        for (unsigned int i = 0; i < ion.size(); i++) {
            if (ion[i].valency > 0)
                type = "1";
            else
                type = "-1";
            outdump << setw(6) << i << "\t" << type << "\t" << setw(8) << ion[i].posvec.x << "\t" << setw(8)
                    << ion[i].posvec.y << "\t" << setw(8) << ion[i].posvec.z << endl;
        }
        outdump.close();
    }
    return;
}

// compute additional quantities
void compute_n_write_useful_data(int cpmdstep, vector<PARTICLE> &ion, vector<THERMOSTAT> &real_bath,
                                 vector<THERMOSTAT> &fake_bath, INTERFACE &box) {

    if (world.rank() == 0) {
        ofstream list_tic("outfiles/total_induced_charge.dat", ios::app);
        ofstream list_temperature("outfiles/temperature.dat", ios::app);
        ofstream list_energy("outfiles/energy.dat", ios::app);
        list_temperature << cpmdstep << setw(15) << 2 * particle_kinetic_energy(ion) / (real_bath[0].dof * kB)
                         << setw(15)
                         << real_bath[0].T << setw(15) << 2 * fake_kinetic_energy(box) / (fake_bath[0].dof * kB)
                         << setw(15)
                         << fake_bath[0].T << endl;
        list_tic << cpmdstep << setw(15) << box.total_induced_charge() << endl;
        double fake_ke = fake_kinetic_energy(box);
        double particle_ke = particle_kinetic_energy(ion);
        double potential_energy = energy_functional(ion, box);
        double real_bath_ke = bath_kinetic_energy(real_bath);
        double real_bath_pe = bath_potential_energy(real_bath);
        double fake_bath_ke = bath_kinetic_energy(fake_bath);
        double fake_bath_pe = bath_potential_energy(fake_bath);
        double extenergy =
                fake_ke + particle_ke + potential_energy + real_bath_ke + real_bath_pe + fake_bath_ke + fake_bath_pe;
        list_energy << cpmdstep << setw(15) << extenergy << setw(15) << particle_ke << setw(15) << potential_energy
                    << setw(15) << particle_ke + potential_energy + real_bath_ke + real_bath_pe << setw(15) << fake_ke
                    << setw(15) << fake_ke + fake_bath_ke + fake_bath_pe << setw(15) << real_bath_ke << setw(15)
                    << real_bath_pe << setw(15) << fake_bath_ke << setw(15) << fake_bath_pe << endl;
/*
  ofstream list_temperature ("outfiles_windows/temperature.dat", ios::app);
  ofstream list_energy ("outfiles_windows/energy.dat", ios::app);
  list_temperature << cpmdstep << setw(15) << 2*particle_kinetic_energy(ion)/(real_bath[0].dof*kB) << setw(15) << real_bath[0].T << setw(15) << endl;
  double particle_ke = particle_kinetic_energy(ion);
  double potential_energy = energy_functional(ion, box);
  double real_bath_ke = bath_kinetic_energy(real_bath);
  double real_bath_pe = bath_potential_energy(real_bath);
  double extenergy = particle_ke + potential_energy + real_bath_ke + real_bath_pe;
  list_energy << cpmdstep << setw(15) << extenergy << setw(15) << particle_ke << setw(15) << potential_energy << setw(15) << particle_ke + potential_energy + real_bath_ke + real_bath_pe << setw(15) << real_bath_ke << setw(15) << real_bath_pe << endl;*/
        list_tic.close();
        list_temperature.close();
        list_energy.close();
    }

}


void progressBar(double fraction_completed) {

    if (world.rank() == 0) {
        int val = (int) (fraction_completed * 100);
        int lpad = (int) (fraction_completed * PBWIDTH);
        int rpad = PBWIDTH - lpad;
        printf("\r%3d%% |%.*s%*s|", val, lpad, PBSTR, rpad, "");
        fflush(stdout);
    }
}

//integration used to calculate the exact induce density for e_1, e_2, e_1
/*
struct exact_density_params {double a; double b; double c;};

double exact_density_func (double x, void * p)
{
  struct exact_density_params *params = (struct exact_density_params *)p;

  double R = (params -> a);
  double L = (params -> b);
  double e_factor = (params -> c);

  double J0 = gsl_sf_bessel_J0(x*R);
  return (x*J0/(exp(x*L) + e_factor * exp(-x*L)));
}

double exact_density_integration(double R, double L, double e_factor)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(10000);

  double result, error;

  gsl_function func;
  struct exact_density_params params = { R, L, e_factor };

  func.function = &exact_density_func;
  func.params = &params;

  gsl_integration_qagiu(&func, 0, 0,1e-6,10000,w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}
*/

// // compute density profile
// void compute_density_profile(int cpmdstep, double density_profile_samples, vector<double>& mean_density, vector<double>& mean_sq_density, vector<PARTICLE>& ion, INTERFACE& box, vector<BIN>& bin, CONTROL& cpmdremote)
// {
//   vector<double> sample_density;
// //   ofstream file_for_auto_corr ("outfiles/for_auto_corr.dat", ios::app);
//       
//   bin_ions(ion, box, sample_density, bin);
//   for (unsigned int b = 0; b < mean_density.size(); b++)
//     mean_density.at(b) = mean_density.at(b) + sample_density.at(b);
//   for (unsigned int b = 0; b < sample_density.size(); b++)
//     mean_sq_density.at(b) = mean_sq_density.at(b) + sample_density.at(b)*sample_density.at(b);
//   
//   // write a file for post analysis to get auto correlation time		// NOTE this is assuming ions do not cross the interface
// //  if (ion[0].posvec.GetMagnitude() > dsphere.radius)
//  //   file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t" << sample_density[int((dsphere.radius+2)/bin[0].width)] << endl;
//  // else
//  //   file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t" << sample_density[int((dsphere.radius-2)/bin[0].width)] << endl;
// 
//   // write files
//   if (cpmdstep % cpmdremote.writedensity == 0)
//   {
//     char data[200];
//     sprintf(data, "datafiles/_den_%.06d.dat", cpmdstep);
//     ofstream outden;
//     outden.open(data);
//     for (unsigned int b = 0; b < mean_density.size(); b++)
//       outden << (-box.lz/2+b*bin[b].width) << setw(15) << mean_density.at(b)/density_profile_samples << endl;
//     outden.close();
//   } 
//   return;
// }


// verify on the fly properties with exact
double verify_with_FMD(int cpmdstep, vector<PARTICLE> &ion, INTERFACE &box, CONTROL &fmdremote, CONTROL &cpmdremote) {
    INTERFACE box_exact = box;
    fmdremote.verify = cpmdstep;
    for (unsigned int i = 0; i < box_exact.leftplane.size(); i++)
        box_exact.leftplane[i].w = 0.0;
    for (unsigned int i = 0; i < box_exact.rightplane.size(); i++)
        box_exact.rightplane[i].w = 0.0;
    fmd(ion, box_exact, fmdremote, cpmdremote);
    for (unsigned int k = 0; k < box_exact.leftplane.size(); k++)
        box_exact.leftplane[k].w = box_exact.leftplane[k].wmean;
    for (unsigned int k = 0; k < box_exact.rightplane.size(); k++)
        box_exact.rightplane[k].w = box_exact.rightplane[k].wmean;

    double exact_functional = energy_functional(ion, box_exact);
    double on_the_fly_functional = energy_functional(ion, box);
    double functional_deviation = 0;
    functional_deviation = (on_the_fly_functional - exact_functional) / exact_functional;

    if (world.rank() == 0) {
        ofstream track_density_left("outfiles/track_density_left.dat", ios::app);
        ofstream track_density_right("outfiles/track_density_right.dat", ios::app);
        ofstream track_functional("outfiles/track_functional.dat", ios::app);
        ofstream track_functional_deviation("outfiles/track_deviation.dat", ios::app);
        track_density_left << cpmdstep << setw(15) << box.leftplane[0].w << setw(15) << box_exact.leftplane[0].w
                           << endl;
        track_density_right << cpmdstep << setw(15) << box.rightplane[0].w << setw(15) << box_exact.rightplane[0].w
                            << endl;
        track_functional << cpmdstep << setw(15) << on_the_fly_functional << setw(15) << exact_functional << endl;
        track_functional_deviation << cpmdstep << setw(15) << functional_deviation << endl;


        // write exact induced density
        char data_left[200], data_right[200];
        sprintf(data_left, "outfiles/verifiles/_ind_%.06d_left.dat", cpmdstep);
        ofstream out_correct_ind_left;
        out_correct_ind_left.open(data_left);
        for (unsigned int k = 0; k < box_exact.leftplane.size(); k++)
            out_correct_ind_left << k + 1 << setw(15) << box_exact.leftplane[k].posvec << setw(15)
                                 << sqrt(box_exact.leftplane[k].posvec.x * box_exact.leftplane[k].posvec.x +
                                         box_exact.leftplane[k].posvec.y * box_exact.leftplane[k].posvec.y) << setw(15)
                                 << box_exact.leftplane[k].wmean << endl;
        sprintf(data_right, "outfiles/verifiles/_ind_%.06d_right.dat", cpmdstep);
        ofstream out_correct_ind_right;
        out_correct_ind_right.open(data_right);
        for (unsigned int k = 0; k < box_exact.rightplane.size(); k++)
            out_correct_ind_right << k + 1 << setw(15) << box_exact.rightplane[k].posvec << setw(15)
                                  << sqrt(box_exact.rightplane[k].posvec.x * box_exact.rightplane[k].posvec.x +
                                          box_exact.rightplane[k].posvec.y * box_exact.rightplane[k].posvec.y)
                                  << setw(15)
                                  << box_exact.rightplane[k].wmean << endl;

        // write cpmd computed induced density
        sprintf(data_left, "outfiles/verifiles/_cpmdind_%.06d_left.dat", cpmdstep);
        ofstream out_cpmd_ind_left;
        out_cpmd_ind_left.open(data_left);
        for (unsigned int k = 0; k < box.leftplane.size(); k++)
            out_cpmd_ind_left << k + 1 << setw(15) << box.leftplane[k].posvec << setw(15)
                              << sqrt(box.leftplane[k].posvec.x * box.leftplane[k].posvec.x +
                                      box.leftplane[k].posvec.y * box.leftplane[k].posvec.y) << setw(15)
                              << box.leftplane[k].w << endl;
        sprintf(data_right, "outfiles/verifiles/_cpmdind_%.06d_right.dat", cpmdstep);
        ofstream out_cpmd_ind_right;
        out_cpmd_ind_right.open(data_right);
        for (unsigned int k = 0; k < box.rightplane.size(); k++)
            out_cpmd_ind_right << k + 1 << setw(15) << box.rightplane[k].posvec << setw(15)
                               << sqrt(box.rightplane[k].posvec.x * box.rightplane[k].posvec.x +
                                       box.rightplane[k].posvec.y * box.rightplane[k].posvec.y) << setw(15)
                               << box.rightplane[k].w << endl;


        track_density_left.close();
        track_density_right.close();
        track_functional.close();
        track_functional_deviation.close();
        out_correct_ind_left.close();
        out_correct_ind_right.close();
        out_cpmd_ind_left.close();
        out_cpmd_ind_right.close();
    }


    return functional_deviation;
}

// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm) {
    char filename[200];
    sprintf(filename, "outfiles/energy.dat");
    ifstream in(filename, ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        return 0;
    }

    int col1;
    double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
    vector<double> ext, ke, pe, fake, real;
    while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11) {
//     if (col1 < hiteqm) continue;
        ext.push_back(col2);
        ke.push_back(col3);
        pe.push_back(col4);
        fake.push_back(col6);
        real.push_back(col5);
    }

//   cout << "Sizes of ext and ke arrays" << setw(10) << ext.size() << setw(10) << ke.size() << endl;

    double ext_mean = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_mean += ext[i];
    ext_mean = ext_mean / ext.size();
    double ke_mean = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_mean += ke[i];
    ke_mean = ke_mean / ke.size();

//   cout << "Mean ext and ke" << setw(10) << ext_mean << setw(10) << ke_mean << endl;

    double ext_sd = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
    ext_sd = ext_sd / ext.size();
    ext_sd = sqrt(ext_sd);

    double ke_sd = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
    ke_sd = ke_sd / ke.size();
    ke_sd = sqrt(ke_sd);

//   cout << "Standard deviations in ext and ke" << setw(10) << ext_sd << setw(10) << ke_sd << endl;

    double R = ext_sd / ke_sd;
//   cout << "R" << setw(15) <<  R << endl;
    if (world.rank() == 0) {
        ofstream out("outfiles/R.dat");
        out << "Sample size " << ext.size() << endl;
        out << "Sd: ext, kinetic energy and R" << endl;
        out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;

    }

    return R;
}
/*
// auto correlation function
void auto_correlation_function()
{
  char filename[200];
  sprintf(filename, "outfiles_windows/for_auto_corr.dat");
  ifstream in(filename, ios::in);
  if (!in) 
  {
    cout << "File could not be opened" << endl; 
    return;
  }

  double col1, col2;
  vector<double> n, autocorr;
  while (in >> col1 >> col2)
    n.push_back(col2);
//   cout << "Number of samples after eqm" << setw(10) << n.size() << endl;

  double avg = 0;
  for (unsigned int j = 0; j< n.size(); j++)
    avg = avg + n[j];
  avg = avg / n.size();

  int ntau = 5000;			// time to which the auto correlation function is computed 
  
  for (int i = 0; i < ntau; i++)
  {
    double A = 0;
    for (unsigned int j = 0; j< n.size(); j++)
      A = A + n[j+i]*n[j];
    A = A / n.size();
    autocorr.push_back(A - avg*avg);
  }

  ofstream out ("outfiles_windows/auto_correlation.dat");
  for (int i = 0; i < ntau; i++)
    out << i << setw(15) << autocorr[i]/autocorr[0] << endl;
  
  cout << "Auto correlation function generated" << endl;
  return;
}
*/
