/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2014
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
b  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with the R package npcp. If not, see <http://www.gnu.org/licenses/>.
  ##
  #################################################################################*/


#ifndef SET_UTILS_H
#define SET_UTILS_H

/**
 * @file   set_utils.h
 * @author Michel Grabisch and Ivan Kojadinovic
 * @date   May 2007
 *
 * @brief  Set function utilities adapted from the R package kappalab
 *
 *
 */

int card(int n);
double sum_binom(int n, int k);
void k_power_set(int *n, int *k, int *power_set);
void k_power_set_char(int *n, int *k, int *k_power_set, char **subset);
void natural2binary(int *n, double *sf, int *power_set, double *sf_out);

#endif
