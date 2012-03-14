/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at
  your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
  Normal and Gamma random number generators. They are stripped from GSL
  library and further modified by me.

  - Attractive Chaos
 */

#ifndef LH3_GENRAN_H_
#define LH3_GENRAN_H_

#include <stdlib.h>
#include <time.h>

#ifndef _WIN32 /* POSIX: rand48 family */

#include <sys/types.h>
#include <unistd.h>

#define ran_seed() srand48(time(0) * (long)getpid())
#define ran_uniform() drand48()

#else /* Windows: this will be pretty BAD. */

#define ran_seed() srand(time(0))
#define ran_uniform() ((double)rand() / RAND_MAX)

#endif

#ifdef __cpluplus
extern "C" {
#endif

	double ran_normal(); /* Ziggurat method */
	double ran_normal2(); /* slower implementation */
	/* g(x;a,b) = x^{a-1} b^a e^{-bx} / \Gamma(a) */
	double ran_gamma(double alpha, double beta);

#ifdef __cplusplus
}
#endif

#endif
