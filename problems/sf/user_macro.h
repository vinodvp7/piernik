
#include "piernik.h"
/* #include "macros.h" */

#if defined MASS_COMPENS || defined OVLP_ON_BNDS
#define INIT_STATE_REFERENCE
#endif /* MASS_COMPENS || OVLP_ON_BNDS */

#if defined MASS_COMPENS
#define MC_OR_OVLP
#endif /* MASS_COMPENS */

#if defined(SNE_DISTR) && !defined(SN_DISTRIBUTION) && !defined(RANDOMIZE)
#define RANDOMIZE
#endif /* SNE_DISTR ** !SN_DISTRIBUTION && !RANDOMIZE */
