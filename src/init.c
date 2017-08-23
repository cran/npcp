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

/* 
   Automatic generation from Dev/
   tools::package_native_routine_registration_skeleton('npcp',,,FALSE)
*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "npcp.h"

// ./cdftest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestF_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, REALSXP
};

// ./ectest.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cpTestC_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
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

static R_NativePrimitiveArgType influRho_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
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

// ./utilities.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType pdf_sum_unif_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP
};

static const R_CMethodDef CEntries[] = {

    {"cpTestAutocop",  (DL_FUNC) &cpTestAutocop,  12, cpTestAutocop_t},

    {"cpTestBM",       (DL_FUNC) &cpTestBM,       14, cpTestBM_t},

    {"cpTestC",        (DL_FUNC) &cpTestC,        10, cpTestC_t},

    {"cpTestF",        (DL_FUNC) &cpTestF,        12, cpTestF_t},

    {"cpTestMean",     (DL_FUNC) &cpTestMean,     10, cpTestMean_t},

    {"cpTestRho",      (DL_FUNC) &cpTestRho,      14, cpTestRho_t},

    {"cpTestU",        (DL_FUNC) &cpTestU,        11, cpTestU_t},

    {"fitGEV",         (DL_FUNC) &fitGEV,          9, fitGEV_t},

    {"influRho",       (DL_FUNC) &influRho,        5, influRho_t},

    {"k_power_set",    (DL_FUNC) &k_power_set,     3, k_power_set_t},

    {"natural2binary", (DL_FUNC) &natural2binary,  4, natural2binary_t},

    {"pdf_sum_unif",   (DL_FUNC) &pdf_sum_unif,    4, pdf_sum_unif_t},

    {NULL, NULL, 0}
};

void R_init_npcp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

