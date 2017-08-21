/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2017
  ##
  ##   This file is part of the R package npcp.
  ##
  ##   The R package npcp is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 3 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package npcp is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with the R package npcp. If not, see <http://www.gnu.org/licenses/>.
  ##
  #################################################################################*/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "npcp.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


// ./cdftest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestF_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, REALSXP
};

// ./ectest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestC_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
};
static R_NativePrimitiveArgType cpTestAutocop_t[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, INTSXP
};

// ./set_utils.c ///////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType k_power_set_t[] = {
    INTSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType natural2binary_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP
};

// ./rhotest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestRho_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP
};

// ./meantest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestMean_t[] = {
    REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP
};

// ./Utest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestU_t[] = {
    REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, REALSXP
};

// ./bmtest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestBM_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType fitGEV_t[] = {
    REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP
};


static const R_CMethodDef CEntries[]  = {
    CDEF(cpTestF),

    CDEF(cpTestC),
    CDEF(cpTestAutocop),

    CDEF(k_power_set),
    CDEF(natural2binary),

    CDEF(cpTestRho),

    CDEF(cpTestMean),

    CDEF(cpTestU),

    CDEF(cpTestBM),
    CDEF(fitGEV),

    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {

    //CALLDEF(rF01Frank_vec_c, 5),

    {NULL, NULL, 0}
};

/**
 * register routines
 * @param dll pointer
 * @return none
 * @author Martin Maechler
 */
void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_copula(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
