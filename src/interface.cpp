// This file contains member functions for interface class

#include "interface.h"
#include "functions.h"

void INTERFACE::set_up(double salt_conc_in, double salt_conc_left, double salt_conc_right, int salt_valency_in,
                       int salt_valency_left, int salt_valency_right, double bx, double by, double bz) {
    // useful combinations of different dielectric constants (inside and outside)
    em_left = 0.5 * (ein + e_left);
    ed_left = (ein - e_left) / (4 * pi);
    em_right = 0.5 * (ein + e_right);
    ed_right = (ein - e_right) / (4 * pi);

    // useful length scales signifying competition between electrostatics and entropy
    lB_in = (lB_water * epsilon_water / ein) / unitlength;
    lB_left = (lB_water * epsilon_water / e_left) / unitlength;
    lB_right = (lB_water * epsilon_water / e_right) / unitlength;

    if (salt_conc_in != 0) {
        inv_kappa_in = (0.257 / (salt_valency_in * sqrt(lB_in * unitlength * salt_conc_in))) / unitlength;
        mean_sep_in = pow(1.2 * salt_conc_in, -1.0 / 3.0) / unitlength;
    } else {
        inv_kappa_in = 0;
        mean_sep_in = 0;
    }
    if (salt_conc_left != 0) {
        inv_kappa_left = (0.257 / (salt_valency_left * sqrt(lB_left * unitlength * salt_conc_left))) / unitlength;
        mean_sep_left = pow(1.2 * salt_conc_left, -1.0 / 3.0) / unitlength;
    } else {
        inv_kappa_left = 0;
        mean_sep_left = 0;
    }
    if (salt_conc_right != 0) {
        inv_kappa_right = (0.257 / (salt_valency_right * sqrt(lB_right * unitlength * salt_conc_right))) / unitlength;
        mean_sep_right = pow(1.2 * salt_conc_right, -1.0 / 3.0) / unitlength;
    } else {
        inv_kappa_right = 0;
        mean_sep_right = 0;
    }
    // simulation box size
    lx = bx;
    ly = by;
    lz = bz;

    return;
}

void INTERFACE::put_one_saltion_in_center(vector<PARTICLE> &saltion_in, int pz, int nz, double diameter,
                                          vector<PARTICLE> &ion) {

    // express diameter in consistent units
    diameter = diameter / unitlength;
    //VECTOR3D posvec_positive = VECTOR3D(27.7265,-16.6373,-0.0962194);
    //VECTOR3D posvec_negative = VECTOR3D(28.3636,15.1469,0.256414);
    VECTOR3D posvec_positive = VECTOR3D(0, 0, 2);
    VECTOR3D posvec_negative = VECTOR3D(0, 0, -3);
    PARTICLE freshion;


    freshion = PARTICLE(int(ion.size()) + 1, diameter, pz, pz * 1.0, 1.0, ein, posvec_positive, lx, ly, lz);
    saltion_in.push_back(freshion);        // create a salt ion
    ion.push_back(freshion);            // copy the salt ion to the stack of all ions

    freshion = PARTICLE(int(ion.size()) + 1, diameter, nz, nz * 1.0, 1.0, ein, posvec_negative, lx, ly, lz);
    saltion_in.push_back(freshion);        // create a salt ion
    ion.push_back(freshion);            // copy the salt ion to the stack of all ions

    if (world.rank() == 0) {
        ofstream list_salt_ions_inside("outfiles/salt_ions_inside.xyz", ios::out);
        list_salt_ions_inside << saltion_in.size() << endl;
        list_salt_ions_inside << "salt ions inside" << endl;

        for (unsigned int i = 0; i < saltion_in.size(); i++)
            list_salt_ions_inside << "Si" << setw(15) << saltion_in[i].posvec << endl;
        list_salt_ions_inside.close();
    }


    return;
}

void INTERFACE::put_saltions_inside(vector<PARTICLE> &saltion_in, int pz, int nz, double concentration, double diameter,
                                    vector<PARTICLE> &ion) {
    // establish the number of inside salt ions first
    // Note: salt concentration is the concentration of one kind of ions, so for total ions a factor of 2 needs to be multiplied.
    // also the factor of 0.6 is there in order to be consistent with units.

    double volume_box = lx * ly * lz;
    unsigned int total_nions_inside = int(
            (concentration * 0.6022) * (volume_box * unitlength * unitlength * unitlength));
    if (total_nions_inside % pz != 0)
        total_nions_inside = total_nions_inside - (total_nions_inside % pz) + pz;

    unsigned int total_pions_inside = abs(nz) * total_nions_inside / pz;
    unsigned int total_saltions_inside = total_nions_inside + total_pions_inside;

    // express diameter in consistent units
    diameter = diameter / unitlength;

    // distance of closest approach between the ion and the interface
    double r0_x = 0.5 * lx - 0.5 * diameter;
    double r0_y = 0.5 * ly - 0.5 * diameter;
    double r0_z = 0.5 * lz - 0.5 * diameter;

//   UTILITY ugsl;													// utility used for making initial configuration 

    const gsl_rng_type *rnT;
    gsl_rng *rnr;

    rnT = gsl_rng_default;
    rnr = gsl_rng_alloc(rnT);

//   unsigned long int s = 23897897;  // gsl_rng_uniform will eventually
    // want a non-negative "long" integer
    gsl_rng_env_setup();

    // if you want to start with the same initial configuration comment three lines below
    // else uncomment them. uncommenting will generate a new randomized distribution of salt ions every run
    // srand((time(0)));                // srand & time are built-in
    // unsigned long int s = random();  // gsl_rng_uniform will eventually
    // gsl_rng_set(rnr,s); // seed the random number generator;

    // generate salt ions inside
    while (saltion_in.size() != total_saltions_inside) {
        double x = gsl_rng_uniform(rnr);
        x = (1 - x) * (-r0_x) + x * (r0_x);
        double y = gsl_rng_uniform(rnr);
        y = (1 - y) * (-r0_y) + y * (r0_y);
        double z = gsl_rng_uniform(rnr);
        z = (1 - z) * (-r0_z) + z * (r0_z);
        VECTOR3D posvec = VECTOR3D(x, y, z);
        if (x > r0_x - diameter || y > r0_y - diameter ||
            z > r0_z - diameter)                        // putting an extra ion diameter length away from interface
            continue;
        bool continuewhile = false;
        for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
            if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5 * diameter + 0.5 * ion[i].diameter))
                continuewhile = true;
        if (continuewhile == true)
            continue;
        PARTICLE freshion;
        if (saltion_in.size() < total_pions_inside)
            freshion = PARTICLE(int(ion.size()) + 1, diameter, pz, pz * 1.0, 1.0, ein, posvec, lx, ly, lz);
        else
            freshion = PARTICLE(int(ion.size()) + 1, diameter, nz, nz * 1.0, 1.0, ein, posvec, lx, ly, lz);
        saltion_in.push_back(freshion);        // create a salt ion
        ion.push_back(freshion);            // copy the salt ion to the stack of all ions
    }

    if (world.rank() == 0) {
        ofstream list_salt_ions_inside("outfiles/salt_ions_inside.xyz", ios::out);
        list_salt_ions_inside << saltion_in.size() << endl;
        list_salt_ions_inside << "salt ions inside" << endl;
        for (unsigned int i = 0; i < saltion_in.size(); i++)
            list_salt_ions_inside << "Si" << setw(15) << saltion_in[i].posvec << endl;
        list_salt_ions_inside.close();
    }

    gsl_rng_free(rnr);

    return;
}

// discretize interface
void INTERFACE::discretize(double ion_diameter, double f) {
    // width of the discretization (f is typically 1 or 1/2 or 1/4 or 1/8 ...)
    width = f * ion_diameter;

    unsigned int nx = int(lx / width);
    unsigned int ny = int(ly / width);

    // creating a discretized hard wall interface at z = - l/2
    for (unsigned int j = 0; j < ny; j++) {
        for (unsigned int i = 0; i < nx; i++) {
            VECTOR3D position = VECTOR3D(-0.5 * lx + 0.5 * width + i * width, -0.5 * ly + 0.5 * width + j * width,
                                         -0.5 * lz);
            double area = width * width;
            VECTOR3D normal = VECTOR3D(0, 0, 1);
            leftplane.push_back(VERTEX(position, area, normal));
        }
    }

    // creating a discretized hard wall interface at z = l/2
    for (unsigned int j = 0; j < ny; j++) {
        for (unsigned int i = 0; i < nx; i++) {
            VECTOR3D position = VECTOR3D(-0.5 * lx + 0.5 * width + i * width, -0.5 * ly + 0.5 * width + j * width,
                                         0.5 * lz);
            double area = width * width;
            VECTOR3D normal = VECTOR3D(0, 0, -1);
            rightplane.push_back(VERTEX(position, area, normal));
        }
    }

    ndot_grad_greens_innersum_left.resize(leftplane.size());
    greens_innersum_left.resize(leftplane.size());
    ndot_grad_greens_innersum_right.resize(rightplane.size());
    greens_innersum_right.resize(rightplane.size());
    if (world.rank() == 0) {
        ofstream listleftplane("outfiles/leftplane.xyz", ios::out);
        listleftplane << leftplane.size() << endl;
        listleftplane << "leftinterface" << endl;
        for (unsigned int k = 0; k < leftplane.size(); k++)
            listleftplane << "I" << setw(15) << leftplane[k].posvec << endl;
        listleftplane.close();

        ofstream listrightplane("outfiles/rightplane.xyz", ios::out);
        listrightplane << rightplane.size() << endl;
        listrightplane << "rightinterface" << endl;
        for (unsigned int k = 0; k < rightplane.size(); k++)
            listrightplane << "I" << setw(15) << rightplane[k].posvec << endl;
        listrightplane.close();
    }
    return;
}

void INTERFACE::put_counterions(vector<PARTICLE> &counterion, int valency, double diameter, vector<PARTICLE> &ion) {
/*  // establish the number of counterions first
  unsigned int total_counterions = int(abs(colloid.valency/valency));
  // express diameter in consistent units
  diameter = diameter / unitlength;
  // distance of closest approach between counterion and the interface (colloid)
  double r0 = 0.5 * colloid.diameter + 0.5 * diameter;					// this would the same as radius + 0.5 * diameter
  // distance of closest approach between counterion and the box
  double r0_x = box_x - 0.5 * diameter;
  double r0_y = box_y - 0.5 * diameter;
  double r0_z = box_z - 0.5 * diameter;
  
  UTILITY ugsl;
  
  // generate counterions in the box
  while (counterion.size() != total_counterions)
  {
    double x = gsl_rng_uniform(ugsl.r);
    x = (1 - x) * (-r0_x) + x * (r0_x);
    double y = gsl_rng_uniform(ugsl.r);
    y = (1 - y) * (-r0_y) + y * (r0_y);
    double z = gsl_rng_uniform(ugsl.r);
    z = (1 - z) * (-r0_z) + z * (r0_z);
    VECTOR3D posvec = VECTOR3D(x,y,z);
    if (posvec.GetMagnitude() < r0 + diameter)						// giving an extra ion diameter length to be safe with the initial generation 
      continue;
    if (posvec.GetMagnitude() >= box_radius - diameter)
      continue;
    bool continuewhile = false;
    for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
      if ((posvec - ion[i].posvec).GetMagnitude() <= diameter) continuewhile = true;		// avoid overlap of initial ion positions
    if (continuewhile == true)
      continue;
    counterion.push_back(PARTICLE(int(ion.size())+1,diameter,valency,valency*1.0,1.0,eout,posvec));		// create a counterion
    ion.push_back(PARTICLE(int(ion.size())+1,diameter,valency,valency*1.0,1.0,eout,posvec));			// copy the counterion as a general ion
  }
  ofstream listcounterions("outfiles/counterions.xyz", ios::out);
  listcounterions << counterion.size() << endl;
  listcounterions << "counterions" << endl;
  for (unsigned int i = 0; i < counterion.size(); i++)
    listcounterions << "C" << setw(15) << counterion[i].posvec << setw(15) << counterion[i].posvec.GetMagnitude() << endl;
  listcounterions.close();
*/
    return;
}
/*
// computing exact induced charge density for a set of charges near a sphere
vector<VERTEX> INTERFACE::compute_exact_induced_density(vector<PARTICLE>& ion, char plane)
{
  vector<VERTEX> exact;
  double prefactor = 0;
  
  if (plane == 'l')
  {
    exact = leftplane;
    prefactor = - (e_left - ein  ) /(2*pi * ein *( e_left + ein) );
  }
  else if (plane == 'r')
  {
    exact = rightplane;
    prefactor = - (e_right - ein  ) /(2*pi * ein *( e_right + ein) );
  }
  
  if (e_left == e_right && e_left != ein)
  {
    double  sum, R;
    for (unsigned int k = 0; k < exact.size(); k++)
    {
      exact[k].w = 0;					// initialize the exact result to zero
      R = sqrt(exact[k].posvec.x*exact[k].posvec.x + exact[k].posvec.y*exact[k].posvec.y);
      for (unsigned int i = 0; i < ion.size(); i ++)	// each ion will contribute to the induced charge at the kth grid point
      {
	sum = exact_density_integration(R, 0.5*lz, ((e_left - ein)/(ein + e_left)));
	exact[k].w += prefactor * sum;
      }
    }
  }
  
  else
  {
    double  sum;
    for (unsigned int k = 0; k < exact.size(); k++)
    {
      exact[k].w = 0;					// initialize the exact result to zero
      for (unsigned int i = 0; i < ion.size(); i ++)	// each ion will contribute to the induced charge at the kth grid point
      {
	sum = ion[i].q *fabs (ion[i].posvec.z - exact[k].posvec.z) /pow((ion[i].posvec - exact[k].posvec).GetMagnitude(),3);
	exact[k].w += prefactor * sum;
      }
    }
  }
  return exact;
}

double INTERFACE::exact_total_energy(vector<VERTEX>& exact_left, vector<VERTEX>& exact_right,vector<PARTICLE>& ion)
{
  double total_energy=0, energy_left=0;
  for (unsigned int i = 0; i< ion.size(); i ++ )
    {
      for (unsigned int k = 0 ; k < exact_left.size(); k ++)
        total_energy += exact_left[k].w * exact_left[k].a/(exact_left[k].posvec - ion[i].posvec).GetMagnitude();
      for (unsigned int k = 0 ; k < exact_right.size(); k ++)
        total_energy += exact_right[k].w * exact_right[k].a/(exact_right[k].posvec - ion[i].posvec).GetMagnitude();

      total_energy = total_energy * 0.5 * ion[i].q;
    }
  if (e_right == ein)
	for (unsigned int i = 0; i < ion.size(); i ++)
	  energy_left += - ion[i].q * ion[i].q *(e_left - e_right)/(4 * ein *0.5 * lz *(e_left + ein));
  energy_left = energy_left * scalefactor;
  cout << " energy for left interface turn on only : " << energy_left << endl;
  return (total_energy * scalefactor);
}
*/

// vector<VERTEX> INTERFACE::compute_exact_induced_density_right( vector<PARTICLE>& ion)
// {
//   vector<VERTEX> exact;
//   exact = rightplane;						// assign exact to s to start with, so that it has all the information held by s
//   double  sum;
//   double prefactor = - (e_right - ein  ) /(2*pi * ein *( e_right + ein) );
// 
//   for (unsigned int k = 0; k <rightplane.size(); k++)
//   {
//     exact[k].w = 0;					// initialize the exact result to zero
//     for (unsigned int i = 0; i < ion.size(); i ++)	// each ion will contribute to the induced charge at the kth grid point
//     {
//        sum = ion[i].q * fabs (ion[i].posvec.z - rightplane[k].posvec.z) /pow((ion[i].posvec - rightplane[k].posvec).GetMagnitude(),3);
//       exact[k].w += prefactor * sum;
//     }
//   }
//   return exact;
// }

void
INTERFACE::put_saltions_leftside(vector<PARTICLE> &saltion_left, int valency, double concentration, double diameter,
                                 vector<PARTICLE> &ion) {
    /*  // establish the number of leftside salt ions first
      // Note: salt concentration is the concentration of one kind of ions, so for total ions a factor of 2 needs to be multiplied. also some factors appear to be consistent with units.
      double volume_box = (4.0/3.0) * pi * (box_radius * box_radius * box_radius - radius * radius * radius);
      unsigned int total_saltions_leftside = int(2 * (concentration * 0.6) * (volume_box * unitlength * unitlength * unitlength));		// NOTE there is a change in the factor, i think it was wrong before
      if (total_saltions_leftside % 2 !=0)
        total_saltions_leftside = total_saltions_leftside + 1;
      // express diameter in consistent units
      diameter = diameter / unitlength;
      // distance of closest approach between the ion and the interface
      double r0 = radius + 0.5 * diameter;
      // distance of closest approach between the ion and the box
      double r0_box = box_radius - 0.5 * diameter;

      UTILITY ugsl;														// utility used for making initial configuration

      // generate salt ions leftside
      while (saltion_left.size() != total_saltions_leftside)
      {
        double x = gsl_rng_uniform(ugsl.r);
        x = (1 - x) * (-r0_box) + x * (r0_box);
        double y = gsl_rng_uniform(ugsl.r);
        y = (1 - y) * (-r0_box) + y * (r0_box);
        double z = gsl_rng_uniform(ugsl.r);
        z = (1 - z) * (-r0_box) + z * (r0_box);
        VECTOR3D posvec = VECTOR3D(x,y,z);
        if (posvec.GetMagnitude() < r0 + diameter)										// giving an extra ion diameter length
          continue;
        if (posvec.GetMagnitude() >= r0_box)
          continue;
        bool continuewhile = false;
        for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
          if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5*diameter+0.5*ion[i].diameter)) continuewhile = true;		// avoid overlap of initial ion positions
        if (continuewhile == true)
          continue;
        saltion_left.push_back(PARTICLE(int(ion.size())+1,diameter,valency,valency*1.0,1.0,eout,posvec));			// create a salt ion
        ion.push_back(PARTICLE(int(ion.size())+1,diameter,valency,valency*1.0,1.0,eout,posvec));				// copy the salt ion to the stack of all ions
        valency = (-1) * valency;												// switch between creating positive and negative ion
      }
      ofstream list_salt_ions_leftside("outfiles/salt_ions_leftside.xyz", ios::out);
      list_salt_ions_leftside << saltion_left.size() << endl;
      list_salt_ions_leftside << "salt ions leftside" << endl;
      for (unsigned int i = 0; i < saltion_left.size(); i++)
        list_salt_ions_leftside << "So" << setw(15) << saltion_left[i].posvec << endl;
      list_salt_ions_leftside.close();
      */
    return;
}

void
INTERFACE::put_saltions_rightside(vector<PARTICLE> &saltion_right, int valency, double concentration, double diameter,
                                  vector<PARTICLE> &ion) {
    /*  // establish the number of rightside salt ions first
      // Note: salt concentration is the concentration of one kind of ions, so for total ions a factor of 2 needs to be multiplied. also some factors appear to be consistent with units.
      double volume_box = (4.0/3.0) * pi * (box_radius * box_radius * box_radius - radius * radius * radius);
      unsigned int total_saltions_rightside = int(2 * (concentration * 0.6) * (volume_box * unitlength * unitlength * unitlength));		// NOTE there is a change in the factor, i think it was wrong before
      if (total_saltions_rightside % 2 !=0)
        total_saltions_rightside = total_saltions_rightside + 1;
      // express diameter in consistent units
      diameter = diameter / unitlength;
      // distance of closest approach between the ion and the interface
      double r0 = radius + 0.5 * diameter;
      // distance of closest approach between the ion and the box
      double r0_box = box_radius - 0.5 * diameter;

      UTILITY ugsl;														// utility used for making initial configuration

      // generate salt ions rightside
      while (saltion_right.size() != total_saltions_rightside)
      {
        double x = gsl_rng_uniform(ugsl.r);
        x = (1 - x) * (-r0_box) + x * (r0_box);
        double y = gsl_rng_uniform(ugsl.r);
        y = (1 - y) * (-r0_box) + y * (r0_box);
        double z = gsl_rng_uniform(ugsl.r);
        z = (1 - z) * (-r0_box) + z * (r0_box);
        VECTOR3D posvec = VECTOR3D(x,y,z);
        if (posvec.GetMagnitude() < r0 + diameter)										// giving an extra ion diameter length
          continue;
        if (posvec.GetMagnitude() >= r0_box)
          continue;
        bool continuewhile = false;
        for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
          if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5*diameter+0.5*ion[i].diameter)) continuewhile = true;		// avoid overlap of initial ion positions
        if (continuewhile == true)
          continue;
        saltion_right.push_back(PARTICLE(int(ion.size())+1,diameter,valency,valency*1.0,1.0,eout,posvec));			// create a salt ion
        ion.push_back(PARTICLE(int(ion.size())+1,diameter,valency,valency*1.0,1.0,eout,posvec));				// copy the salt ion to the stack of all ions
        valency = (-1) * valency;												// switch between creating positive and negative ion
      }
      ofstream list_salt_ions_rightside("outfiles/salt_ions_rightside.xyz", ios::out);
      list_salt_ions_rightside << saltion_right.size() << endl;
      list_salt_ions_rightside << "salt ions rightside" << endl;
      for (unsigned int i = 0; i < saltion_right.size(); i++)
        list_salt_ions_rightside << "So" << setw(15) << saltion_right[i].posvec << endl;
      list_salt_ions_rightside.close();
      */
    return;
}
