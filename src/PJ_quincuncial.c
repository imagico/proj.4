/*
 * Copyright (c) 2016 Christoph Hormann
 * 
 * based on tcl implementation (mapproj.tcl):
 * https://core.tcl.tk/tcllib/finfo?name=modules/mapproj/mapproj.tcl
 * Copyright (c) 2007 by Kevin B. Kenny
 * and licensed under the following terms:
 * 
 * This software is copyrighted by Ajuba Solutions and other parties.
 * The following terms apply to all files associated with the software unless
 * explicitly disclaimed in individual files.
 * 
 * The authors hereby grant permission to use, copy, modify, distribute,
 * and license this software and its documentation for any purpose, provided
 * that existing copyright notices are retained in all copies and that this
 * notice is included verbatim in any distributions. No written agreement,
 * license, or royalty fee is required for any of the authorized uses.
 * Modifications to this software may be copyrighted by their authors
 * and need not follow the licensing terms described here, provided that
 * the new terms are clearly indicated on the first page of each file where
 * they apply.
 * 
 * IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
 * FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
 * ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
 * DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
 * IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 * 
 * GOVERNMENT USE: If you are acquiring this software on behalf of the
 * U.S. government, the Government shall have only "Restricted Rights"
 * in the software and related documentation as defined in the Federal 
 * Acquisition Regulations (FARs) in Clause 52.227.19 (c) (2).  If you
 * are acquiring the software on behalf of the Department of Defense, the
 * software shall be classified as "Commercial Computer Software" and the
 * Government shall have only "Restricted Rights" as defined in Clause
 * 252.227-7013 (c) (1) of DFARs.  Notwithstanding the foregoing, the
 * authors grant the U.S. Government and others acting in its behalf
 * permission to use and distribute the software in accordance with the
 * terms specified in this license.
 */
 
/*
 * ...
 */
#define PROJ_PARMS__ \
	double	sinrot, cosrot;
#define PJ_LIB__
# include   <projects.h>
PROJ_HEAD(pq, "Peirce quincuncial") "\n\tMisc Sph";
#define EPS 1.0e-11
#define TOL 1.e-7
#define SQRT2 sqrt(2.0)
#define PeirceQuincuncialScale 3.7081493546027438
#define PeirceQuincuncialLimit 1.8540746773013719
#define NX 128

/*
 Nextk --
     Auxiliary function for computing next value of k

 Arguments:
     k           Parameter
 Return value:
     Next value to be used
*/

static double Nextk(double k) {
    double ksq = sqrt(1.0-k*k);
    return (1.0-ksq)/(1.0+ksq);
}

/*
 IterateUK --
     Auxiliary function to compute the raw value (phi)

 Arguments:
     u           Independent variable
     k           Parameter
 Return value:
     phi
*/

static double IterateUK(double u, double k) {
    double kvalues[NX+1];
    int i = 0;
    while (k > EPS) {
        k = Nextk(k);
        kvalues[i] = k;
        u = u*2.0/(1.0+k);
        i++;
        if (i == NX) break;
    }
    int j;
    for (j=i-1; j >= 0; j--) {
        u = 0.5*( u + asin(kvalues[j]*sin(u)) );
    }
    return u;
}

/*
 cn --
     Compute the elliptic function cn

 Arguments:
     u           Independent variable
     k           Parameter
 Return value:
     cn(u,k)
 Note:
     If k == 1, then the iteration does not stop
*/

static double cn(double u, double k) {
    if (k == 0.0) return cos(u);
    if (k >= 1.0) return 1.0/cosh(u);

    u = IterateUK(u, k);
    return cos(u);

}

static double cn2(double u, double k)
{
   double a[NX+1];
   double c[NX+1];
   double b;
   double phi;
   int i, i2;

   if (k < EPS) return cos(u);
   if (k >= 1.0-EPS) return 1.0/cosh(u);

   i2 = 1.0;
   b = sqrt(1.0 - k * k);
   a[0] = 1.0;
   c[0] = k;

   for (i = 0; i < NX; i++) {
      if (abs(a[i] - b) < (a[i] * EPS)) break;
      i2 += i2;
      a[i+1] = 0.5 * (a[i] + b);
      c[i+1] = 0.5 * (a[i] - b);
      b = sqrt(a[i] * b);
   }

   phi = i2 * u * a[i];

fprintf(stderr, "phi=%.6f\n", phi);

   for (; i > 0; i--) phi = 0.5 * ( phi + asin( c[i] * sin(phi) / a[i]) );

   return cos(phi);
}



/*
 ::mapproj::ellRF --

	Computes the Carlson incomplete elliptic integral of the
	first kind:

	RF(x, y, z) = 1/2 * integral_0^inf dt/sqrt((t+x)*(t+y)*(t+z))

 Parameters:
	x, y, z -- Interchangeable parameters of the integral

 Results:
	Returns the value of RF
*/
static double ellRF (double x, double y, double z) {
    if ((x < 0.0) || (y < 0.0) || (z < 0.0)) {
        return 0.0;
    }
    double delx = 1.0;
    double dely = 1.0;
    double delz = 1.0;
    double mean;
    while((abs(delx) > 0.0025) || (abs(dely) > 0.0025) || (abs(delz) > 0.0025)) {
        double sx = sqrt(x);
        double sy = sqrt(y);
        double sz = sqrt(z);
        double len = sx * (sy + sz) + sy * sz;
        x = 0.25 * (x + len);
        y = 0.25 * (y + len);
        z = 0.25 * (z + len);
        mean = (x + y + z) / 3.0;
        delx = (mean - x) / mean;
        dely = (mean - y) / mean;
        delz = (mean - z) / mean;
    }
    double e2 = delx * dely - delz * delz;
    double e3 = delx * dely * delz;
    return (1.0 + (e2 / 24.0 - 0.1 - 3.0 * e3 / 44.0) * e2 + e3 / 14.0) / sqrt(mean);
}

/*
 ::mapproj::ellFaux -

	Computes the Legendre incomplete elliptic integral of the
	first kind when circular functions of the 'phi' argument
	are already available.

 Parameters:
	cos_phi - Cosine of the argument
	sin_phi - Sine of the argument
	k - Parameter

 Results:
	Returns F(atan(sin_phi/cos_phi), k)
*/
static double ellFaux (double cos_phi, double sin_phi, double k) {
    double rf = ellRF(cos_phi * cos_phi, 1.0 - k*k * sin_phi*sin_phi, 1.0);
    return sin_phi * rf;
}

FORWARD(s_forward); /* spheroid */
    double lam = lp.lam;
    double phi = lp.phi;

    /* Compute the auxiliary quantities 'm' and 'n'. Set 'm' to match */
    /* the sign of 'lambda' and 'n' to be positive if |lambda| > pi/2 */

    double cos_phiosqrt2 = 0.5*SQRT2*cos(phi);
    double cos_lambda = cos(lam);
    double sin_lambda = sin(lam);
    double cos_a = cos_phiosqrt2 * (sin(lam) + cos(lam));
    double cos_b = cos_phiosqrt2 * (sin(lam) - cos(lam));
    double sin_a = sqrt(1.0 - cos_a * cos_a);
    double sin_b = sqrt(1.0 - cos_b * cos_b);
    double cos_a_cos_b = cos_a * cos_b;
    double sin_a_sin_b = sin_a * sin_b;
    double sin2_m = 1.0 + cos_a_cos_b - sin_a_sin_b;
    double sin2_n = 1.0 - cos_a_cos_b - sin_a_sin_b;
    if (sin2_m < 0.0) sin2_m = 0.0;
    double sin_m = sqrt(sin2_m);
    if (sin2_m > 1.0) sin2_m = 1.0;
    double cos_m = sqrt(1.0 - sin2_m);
    if (sin_lambda < 0.0) sin_m = -sin_m;

    if (sin2_n < 0.0) sin2_n = 0.0;
    double sin_n = sqrt(sin2_n);
    if (sin2_n > 1.0) sin2_n = 1.0;
    double cos_n = sqrt(1.0 - sin2_n);
    if (cos_lambda > 0.0) sin_n = -sin_n;

    /* Compute elliptic integrals to map the disc to the square */

    double x = ellFaux(cos_m, sin_m, 0.5*SQRT2);
    double y = ellFaux(cos_n, sin_n, 0.5*SQRT2);

    /* Reflect the Southern Hemisphere outward */

    if (phi < 0.0) {
        if (lam < -3.0*FORTPI) {
            y = PeirceQuincuncialScale - y;
        }
        else if (lam < -FORTPI) {
            x = -PeirceQuincuncialScale - x;
        }
        else if (lam < FORTPI) {
            y = -PeirceQuincuncialScale - y;
        }
        else if (lam < 3.0*FORTPI) {
            x = PeirceQuincuncialScale - x;
        }
        else {
            y = PeirceQuincuncialScale - y;
        }
    }

    /* Rotate the square by 45 degrees to fit the screen better */

    double xo = (x - y) * 0.5*SQRT2;
    double yo = (x + y) * 0.5*SQRT2;

    /* custom map rotation */
    if (fabs(P->sinrot) < TOL) {
        xy.x = xo;
        xy.y = yo;
    }
    else {
        xy.x = xo * P->cosrot + yo * P->sinrot;
        xy.y = yo * P->cosrot - xo * P->sinrot;
    }

    return (xy);
}
INVERSE(s_inverse); /* spheroid */

    double xo, yo;

    /* custom map rotation */
    if (fabs(P->sinrot) < TOL) {
        xo = xy.x;
        yo = xy.y;
    }
    else {
        xo = xy.x * P->cosrot - xy.y * P->sinrot;
        yo = xy.y * P->cosrot + xy.x * P->sinrot;
    }

    /* handle large scale periodicity */
    /*
       FIXME: this can certainly be done more efficiently
    */
    while ((xo < -PeirceQuincuncialLimit*SQRT2) || (xo >= PeirceQuincuncialLimit*SQRT2)) {
        if (xo < -PeirceQuincuncialLimit*SQRT2)
            xo = -PeirceQuincuncialLimit*2.0*SQRT2 - xo;
        else
            xo = PeirceQuincuncialLimit*2.0*SQRT2 - xo;
        yo = -yo;
    }
    while ((yo < -PeirceQuincuncialLimit*SQRT2) || (yo >= PeirceQuincuncialLimit*SQRT2)) {
        if (yo < -PeirceQuincuncialLimit*SQRT2)
            yo = -PeirceQuincuncialLimit*2.0*SQRT2 - yo;
        else
            yo = PeirceQuincuncialLimit*2.0*SQRT2 - yo;
        xo = -xo;
    }

    /* Rotate x and y 45 degrees */
    double x = (xo + yo) * 0.5*SQRT2;
    double y = (yo - xo) * 0.5*SQRT2;

    /* Reflect Southern Hemisphere into the Northern */
    int southern = 0;
    if (x < -PeirceQuincuncialLimit) {
        x = -PeirceQuincuncialScale - x;
        southern = 1;
    }
    else if (x > PeirceQuincuncialLimit) {
        x = PeirceQuincuncialScale - x;
        southern = 1;
    }
    else if (y < -PeirceQuincuncialLimit) {
        y = -PeirceQuincuncialScale - y;
        southern = 1;
    }
    else if (y > PeirceQuincuncialLimit) {
        y = PeirceQuincuncialScale - y;
        southern = 1;
    }

    /* Now we know that latitude will be positive.  If X is negative, then */
    /* longitude will be negative; reflect the Western Hemisphere into the */
    /* Eastern. */
    int western = 0;
    if (x < 0.0) {
        western = 1;
        x = -x;
    }

    /* If Y is positive, the point is in the back hemisphere.  Reflect */
    /* it to the front. */
    int back = 0;
    if (y > 0.0) {
        back = 1;
        y = -y;
    }

    /* Finally, constrain longitude to be less than pi/4, by reflecting across */
    /* the 45 degree meridian. */
    int complement = 0;
    if (x >  -y) {
        complement = 1;
        double t = -x;
        x = -y;
        y = t;
    }

    /* Compute the elliptic functions to map the plane onto the sphere */
    double cnx = cn(x, 0.5*SQRT2);
    double cny = cn(y, 0.5*SQRT2);

    /* Undo the mapping to latitude and longitude */
    double a1 = acos(-cnx * cnx);
    double a2 = acos(cny * cny);
    double b = 0.5 * (a1 + a2);
    double a = 0.5 * (a1 - a2);
    double cos_a = cos(a);
    double cos_b = -cos(b);
    double lam = FORTPI - atan2(cos_b, cos_a);
    double phi = acos(sqrt(cos_b*cos_b+ cos_a*cos_a));

    /* Undo the reflections that were done above, to get correct latitude */
    /* and longitude */
    if (complement) lam = HALFPI - lam;
    if (back) lam = PI - lam;
    if (western) lam = -lam;
    if (southern) phi = -phi;

    lp.lam = lam;
    lp.phi = phi;

    return (lp);

}
FREEUP; if (P) pj_dalloc(P); }
ENTRY0(pq)
	double alpha = 0.0;
        if (pj_param(P->ctx, P->params, "talpha").i != 0)
		alpha = pj_param(P->ctx, P->params, "ralpha").f;
	P->es = 0.;
	P->inv = s_inverse; P->fwd = s_forward;
	P->sinrot = sin(alpha);
	P->cosrot = cos(alpha);
ENDENTRY(P)
