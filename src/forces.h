// This is a header file containing the forces routines on induced charges 
// and particles respectively

#ifndef _FORCES_H
#define _FORCES_H

#include "vertex.h"
#include "particle.h"
#include "interface.h"
#include "functions.h"

void compute_force_fake_dof(vector<PARTICLE>&, INTERFACE&, char);
void compute_force_real_dof(vector<PARTICLE>&, INTERFACE&, char);

#endif 
