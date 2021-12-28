#ifndef GUARD_changerule_h
#define GUARD_changerule_h

#include "rand_gen.h"
#include "info.h"
#include "tree.h"
#include "logging.h"

bool changerule(tree& x, xinfo& xi, dinfo& di, double*phi, pinfo& pi, RNG& gen, Logger logger, std::vector<size_t>& ivcnt);

#endif /* changerule_h */

