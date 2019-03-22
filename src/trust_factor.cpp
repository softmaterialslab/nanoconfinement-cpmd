#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>

using namespace std;

int main(int argc, char* argv[])
{
  char filename[200];
  sprintf(filename, "outfiles/energy.dat");
  ifstream in(filename, ios::in);
  if (!in)
  {
    cout << "File could not be opened" << endl;
    return 0;
  }

  int col1;
  double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
  vector<double> ext, ke, pe;
  while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11)
  {
//     if (col1 < hiteqm) continue;
    ext.push_back(col2);
    ke.push_back(col3);
    pe.push_back(col4);
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
  double pe_mean = 0;
  for (unsigned int i = 0; i < pe.size(); i++)
    pe_mean += pe[i];
  pe_mean = pe_mean / pe.size();

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

  double pe_sd = 0;
  for (unsigned int i = 0; i < pe.size(); i++)
    pe_sd += (pe[i] - pe_mean) * (pe[i] - pe_mean);
  pe_sd = pe_sd / pe.size();
  pe_sd = sqrt(pe_sd);
//   cout << "Standard deviations in ext and ke" << setw(10) << ext_sd << setw(10) << ke_sd << endl;

  double R_ke = ext_sd / ke_sd;
  double R_pe = ext_sd / pe_sd;
//   cout << "R" << setw(15) <<  R << endl;

  cout << "Sample size " << ext.size() << endl;
  cout << "Sd: ext, kinetic energy , potentail energy and R" << endl;
  cout << ext_sd << setw(15) << ke_sd << "	R_ke:	" << R_ke << endl;
  cout << ext_sd << setw(15) << pe_sd << "	R_pe:	" << R_pe << endl;

  return 0;
}

