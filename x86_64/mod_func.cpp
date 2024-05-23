#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _kv_reg(void);
extern void _na_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mechanisms/kv.mod\"");
    fprintf(stderr, " \"mechanisms/na.mod\"");
    fprintf(stderr, "\n");
  }
  _kv_reg();
  _na_reg();
}

#if defined(__cplusplus)
}
#endif
