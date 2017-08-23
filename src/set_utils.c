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
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with the R package npcp. If not, see <http://www.gnu.org/licenses/>.
  ##
  #################################################################################*/

/**
 * @file   set_utils.c
 * @date   May 2007
 *
 * @brief  Set function utilities written while working on package kappalab
 *
 *
 */

#include <R.h>
#include <Rmath.h>

#include "set_utils.h"

/**
 * Set cardinality of n
 *
 * @param n integer representing a set; uses binary notation
 * @return the number of bits in n equal to 1
 */
int card(int n) {
  int i;
  for(i=0; n; n >>= 1)
    i += n & 1;
  return(i);
}

/**
 * Computes sum_{i=0}^k choose(n,i)
 *
 * @param n
 * @param k
 * @return sum_{i=0}^k choose(n,i)
 */
double sum_binom(int n, int k) {
  double s = 1.;
  for (int i=1; i <= k; i++)
    s += choose(n,i);
  return s;
}

/**
 * Generation of the first k + 1 levels of the power set of X
 * in the "natural" order. Recursive function.
 *
 * @param n cardinality of X
 * @param k
 * @param last
 * @param power_set
 * @param b
 */
static
void k_power_set_rec(int n, int k, int last, int *power_set, int *b) {
  /* look for the leftmost 1 in b and start to fill blank cases
     with 1 left from this position */
  int istart = n;

  if(*b > 0)
    while(!(*b & 1<<(istart-1)))
      istart--;
  else
    istart = 0;

  for(int i=istart; i<n; i++) {

    last++;
    *(power_set+last) = *b + (1<<i);
  }

  if(last != (int)sum_binom(n,k) - 1)
    k_power_set_rec(n, k, last, power_set, b+1);
}

/**
 * Generation of the first k + 1 levels of the power set of X
 * in the "natural" order. Wrapps the previous function.
 *
 * @param n cardinality of X
 * @param k
 * @param power_set contains the k + 1 levels of the power set of X
 */
void k_power_set(int *n, int *k, int *power_set) {
  power_set[0] = 0;
  k_power_set_rec(*n, *k, 0, power_set, power_set);
}


/**
 * Converts an integer representing a set into an array ("list of elements")
 *
 * @param n cardinality of X
 * @param b binary code to be converted
 * @param x result in an array (ordered liste of elements)
 */
static void binary2subset(int n, int b, int *x) {
  for(int i=0; i<n; i++)
    if(b & 1 << i) {
      *x = i;
      x++;
    }
}

/**
 * Used to specify the maximum number of characters for printing a subset
 * Clearly suboptimal
 */
#define SET_MAX 4

/**
 * Converts the k power set of X in the "natural" order to char**
 * Function is suboptimal, partly because of SET_MAX
 *
 * @param n cardinality of X
 * @param sb the length of k_power_set
 * @param k_power_set array representing the k power set of X
 * @param subset converted k power set
 */
void k_power_set_char(int *n, int *sb, int *k_power_set, char **subset) {

  subset[0] = (char *) R_alloc(3, sizeof(char));
  sprintf(subset[0],"{}");

  for(int i=1; i<*sb; i++) {
    int j, x[32];

    for(j=0; j<*n; j++)
      x[j]=0;

    binary2subset(*n, k_power_set[i], x);

    subset[i] = (char *) R_alloc(SET_MAX * (*n), sizeof(char));
    sprintf(subset[i],"{%d",x[0]+1);

    for(j=1; j < card(k_power_set[i]); j++) {
      char string[255];
      sprintf(string,",%d", x[j]+1);
      strcat(subset[i],string);
    }

    strcat(subset[i],"}");
  }
}

/**
 * Writes a set function given in "natural" order in the binary order
 *
 * @param n cardinality of X
 * @param sf array representing the set function in "natural" order
 * @param power_set contains the power_set of X in "natural" order
 * @param sf_out array representing the set function in binary order
 */
void natural2binary(int *n, double *sf, int *power_set, double *sf_out) {
  for(int i=0; i<(1<<*n); i++)
    sf_out[power_set[i]] = sf[i];
}

