// This is bin class

#ifndef _BIN_H
#define _BIN_H

#include "utility.h"

class BIN
{
  public:
    
    // members
    int n;		// number of ions in the bin
    double width;	// width of the bin
    double volume;	// volume of the bin
    double lower;	// lower value of bin
    double higher;	// higher value of bin
    
    // member functions
    
    // make a bin
    BIN(int n = 0, double width = 0, double volume = 0, double lower = 0, double higher = 0 ) : n(n), width(width), volume(volume), lower(lower), higher(higher)
    {
    }
    
    // set up bin
    void set_up(int bin_num, double bin_width, double lx, double ly, double lz)
    {
      n = 0;				// initialize number of ions in bin to be zero
      width = bin_width;
      volume = (bin_width * lx * ly) * (unitlength * unitlength * unitlength) * 0.6022;
      lower = 0.5*(-lz) + bin_num * bin_width;
      higher = 0.5*(-lz) + (bin_num + 1) * bin_width;
    }
};

#endif

