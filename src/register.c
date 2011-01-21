#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "h.h"
#include "hinverse.h"


R_CallMethodDef callRoutines[] = {
    { "hNormalCopula", (DL_FUNC) &hNormalCopula, 3 },
    { "hinverseNormalCopula", (DL_FUNC) &hinverseNormalCopula, 3 },
    { NULL, NULL, 0 }
};

void R_init_vines(DllInfo *info)
{
	R_registerRoutines(info, NULL, callRoutines, NULL, NULL);
}
