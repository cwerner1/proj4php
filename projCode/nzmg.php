<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/*
 * *****************************************************************************
  NAME                            NEW ZEALAND MAP GRID

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the New Zealand Map Grid projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.


  ALGORITHM REFERENCES

  1.  Department of Land and Survey Technical Circular 1973/32
  http://www.linz.govt.nz/docs/miscellaneous/nz-map-definition.pdf

  2.  OSG Technical Report 4.1
  http://www.linz.govt.nz/docs/miscellaneous/nzmg.pdf


  IMPLEMENTATION NOTES

  The two references use different symbols for the calculated values. This
  implementation uses the variable names similar to the symbols in reference [1].

  The alogrithm uses different units for delta latitude and delta longitude.
  The delta latitude is assumed to be in units of seconds of arc x 10^-5.
  The delta longitude is the usual radians. Look out for these conversions.

  The algorithm is described using complex arithmetic. There were three
  options:
 * find and use a Javascript library for complex arithmetic
 * write my own complex library
 * expand the complex arithmetic by hand to simple arithmetic

  This implementation has expanded the complex multiplication operations
  into parallel simple arithmetic operations for the real and imaginary parts.
  The imaginary part is way over to the right of the display; this probably
  violates every coding standard in the world, but, to me, it makes it much
  more obvious what is going on.

  The following complex operations are used:
  - addition
  - multiplication
  - division
  - complex number raised to integer power
  - summation

  A summary of complex arithmetic operations:
  (from http://en.wikipedia.org/wiki/Complex_arithmetic)
  addition:       (a + bi) + (c + di) = (a + c) + (b + d)i
  subtraction:    (a + bi) - (c + di) = (a - c) + (b - d)i
  multiplication: (a + bi) x (c + di) = (ac - bd) + (bc + ad)i
  division:       (a + bi) / (c + di) = [(ac + bd)/(cc + dd)] + [(bc - ad)/(cc + dd)]i

  The algorithm needs to calculate summations of simple and complex numbers. This is
  implemented using a for-loop, pre-loading the summed value to zero.

  The algorithm needs to calculate theta^2, theta^3, etc while doing a summation.
  There are three possible implementations:
  - use pow in the summation loop - except for complex numbers
  - precalculate the values before running the loop
  - calculate theta^n = theta^(n-1) * theta during the loop
  This implementation uses the third option for both real and complex arithmetic.

  For example
  psi_n = 1;
  sum = 0;
  for (n = 1; n <=6; n++) {
  psi_nOne = psi_n * psi;       // calculate psi^(n+1)
  psi_n = psi_nOne;
  sum = sum + A[n] * psi_n;
  }


  TEST VECTORS

  NZMG E, N:         2487100.638      6751049.719     metres
  NZGD49 long, lat:      172.739194       -34.444066  degrees

  NZMG E, N:         2486533.395      6077263.661     metres
  NZGD49 long, lat:      172.723106       -40.512409  degrees

  NZMG E, N:         2216746.425      5388508.765     metres
  NZGD49 long, lat:      169.172062       -46.651295  degrees

  Note that these test vectors convert from NZMG metres to lat/long referenced
  to NZGD49, not the more usual WGS84. The difference is about 70m N/S and about
  10m E/W.

  These test vectors are provided in reference [1]. Many more test
  vectors are available in
  http://www.linz.govt.nz/docs/topography/topographicdata/placenamesdatabase/nznamesmar08.zip
  which is a catalog of names on the 260-series maps.


  EPSG CODES

  NZMG     EPSG:27200
  NZGD49   EPSG:4272

  http://spatialreference.org/ defines these as
  ProjFourphp.defs["EPSG:4272"] = "+proj=longlat +ellps=intl +datum=nzgd49 +no_defs ";
  ProjFourphp.defs["EPSG:27200"] = "+proj=nzmg +lat_0=-41 +lon_0=173 +x_0=2510000 +y_0=6023150 +ellps=intl +datum=nzgd49 +units=m +no_defs ";


  LICENSE
  Copyright: Stephen Irons 2008
  Released under terms of the LGPL as per: http://www.gnu.org/copyleft/lesser.html

 * ***************************************************************************** */

/**
  Initialize New Zealand Map Grip projection
 */
class ProjFourphp_ProjNzmg
{

    /**
     * iterations: Number of iterations to refine inverse transform.
     *     0 -> km accuracy
     *     1 -> m accuracy -- suitable for most mapping applications
     *     2 -> mm accuracy
     */
    protected $iterations = 1;

    /**
     * 
     */
    public function init()
    {
        $this->A = array();
        $this->A[1] = +0.6399175073;
        $this->A[2] = -0.1358797613;
        $this->A[3] = +0.063294409;
        $this->A[4] = -0.02526853;
        $this->A[5] = +0.0117879;
        $this->A[6] = -0.0055161;
        $this->A[7] = +0.0026906;
        $this->A[8] = -0.001333;
        $this->A[9] = +0.00067;
        $this->A[10] = -0.00034;

        $this->bRe = array();
        $this->bIm = array();
        $this->bRe[1] = +0.7557853228;
        $this->bIm[1] = 0.0;
        $this->bRe[2] = +0.249204646;
        $this->bIm[2] = +0.003371507;
        $this->bRe[3] = -0.001541739;
        $this->bIm[3] = +0.041058560;
        $this->bRe[4] = -0.10162907;
        $this->bIm[4] = +0.01727609;
        $this->bRe[5] = -0.26623489;
        $this->bIm[5] = -0.36249218;
        $this->bRe[6] = -0.6870983;
        $this->bIm[6] = -1.1651967;

        $this->cRe = array();
        $this->cIm = array();
        $this->cRe[1] = +1.3231270439;
        $this->cIm[1] = 0.0;
        $this->cRe[2] = -0.577245789;
        $this->cIm[2] = -0.007809598;
        $this->cRe[3] = +0.508307513;
        $this->cIm[3] = -0.112208952;
        $this->cRe[4] = -0.15094762;
        $this->cIm[4] = +0.18200602;
        $this->cRe[5] = +1.01418179;
        $this->cIm[5] = +1.64497696;
        $this->cRe[6] = +1.9660549;
        $this->cIm[6] = +2.5127645;

        $this->D = array();
        $this->D[1] = +1.5627014243;
        $this->D[2] = +0.5185406398;
        $this->D[3] = -0.03333098;
        $this->D[4] = -0.1052906;
        $this->D[5] = -0.0368594;
        $this->D[6] = +0.007317;
        $this->D[7] = +0.01220;
        $this->D[8] = +0.00394;
        $this->D[9] = -0.0013;
    }

    /**
      New Zealand Map Grid Forward  - long/lat to x/y
      long/lat in radians
     */
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $deltaLat = $lat - $this->latZero;
        $deltaLon = $lon - $this->longZero;

        // 1. Calculate d_phi and d_psi    ...                          // and d_lambda
        // For this algorithm, delta_latitude is in seconds of arc x 10-5, so we need to scale to those units. Longitude is radians.
        $dPhi = $deltaLat / ProjFourphp::$common->secToRad * 1E-5;
        $dLambda = $deltaLon;
        $dPhiN = 1;  // d_phi^0

        $dPsi = 0;
        for ($n = 1; $n <= 10; $n++) {
            $dPhiN = $dPhiN * $dPhi;
            $dPsi = $dPsi + $this->A[$n] * $dPhiN;
        }

        // 2. Calculate theta
        $thRe = $dPsi;
        $thIm = $dLambda;

        // 3. Calculate z
        $thNRe = 1;
        $thNIm = 0;  // theta^0
        #$th_n_reOne;
        #$th_n_im1;

        $zRe = 0;
        $zIm = 0;
        for ($n = 1; $n <= 6; $n++) {
            $thNReOne = $thNRe * $thRe - $thNIm * $thIm;
            $thNImOne = $thNIm * $thRe + $thNRe * $thIm;
            $thNRe = $thNReOne;
            $thNIm = $thNImOne;
            $zRe = $zRe + $this->bRe[$n] * $thNRe - $this->bIm[$n] * $thNIm;
            $zIm = $zIm + $this->bIm[$n] * $thNRe + $this->bRe[$n] * $thNIm;
        }

        // 4. Calculate easting and northing
        $p->x = ($zIm * $this->a) + $this->xZero;
        $p->y = ($zRe * $this->a) + $this->yZero;

        return $p;
    }

    /**
      New Zealand Map Grid Inverse  -  x/y to long/lat
     */
    public function inverse($p)
    {

        $x = $p->x;
        $y = $p->y;

        $deltaX = $x - $this->xZero;
        $deltaY = $y - $this->yZero;

        // 1. Calculate z
        $zRe = $deltaY / $this->a;
        $zIm = $deltaX / $this->a;

        // 2a. Calculate theta - first approximation gives km accuracy
        $zNRe = 1;
        $zNIm = 0;  // z^0
        $zNReOne;
        $zNIm1;

        $thRe = 0;
        $thIm = 0;
        for ($n = 1; $n <= 6; $n++) {
            $zNReOne = $zNRe * $zRe - $zNIm * $zIm;
            $zNIm1 = $zNIm * $zRe + $zNRe * $zIm;
            $zNRe = $zNReOne;
            $zNIm = $zNIm1;
            $thRe = $thRe + $this->cRe[$n] * $zNRe - $this->cIm[$n] * $zNIm;
            $thIm = $thIm + $this->cIm[$n] * $zNRe + $this->cRe[$n] * $zNIm;
        }

        // 2b. Iterate to refine the accuracy of the calculation
        //        0 iterations gives km accuracy
        //        1 iteration gives m accuracy -- good enough for most mapping applications
        //        2 iterations bives mm accuracy
        for ($i = 0; $i < $this->iterations; $i++) {
            $thNRe = $thRe;
            $thNIm = $thIm;
            $thNReOne;
            $thNImOne;

            $numRe = $zRe;
            $numIm = $zIm;
            for ($n = 2; $n <= 6; $n++) {
                $thNReOne = $thNRe * th_re - $thNIm * $thIm;
                $thNImOne = $thNIm * $thRe + $thNRe * $thIm;
                $thNRe = $thNReOne;
                $thNIm = $thNImOne;
                $numRe = $numRe + ($n - 1) * ($this->bRe[$n] * $thNRe - $this->bIm[$n] * $thNIm);
                $numIm = $numIm + (n - 1) * ($this->bIm[$n] * $thNRe + $this->bRe[$n] * $thNIm);
            }

            $thNRe = 1;
            $thNIm = 0;
            $denRe = $this->bRe[1];
            $denIm = $this->bIm[1];
            for ($n = 2; $n <= 6; $n++) {
                $thNReOne = $thNRe * $thRe - $thNIm * $thIm;
                $thNImOne = $thNIm * $thRe + $thNRe * $thIm;
                $thNRe = $thNReOne;
                $thNIm = $thNImOne;
                $denRe = $denRe + $n * ($this->bRe[$n] * $thNRe - $this->bIm[$n] * $thNIm);
                $denIm = $denIm + $n * ($this->bIm[n] * $thNRe + $this->bRe[$n] * $thNIm);
            }

            // Complex division
            $denTwo = $denRe * $denRe + $denIm * $denIm;
            $thRe = ($numRe * $denRe + $numIm * $denIm) / $denTwo;
            $thIm = ($numIm * $denRe - $numRe * $denIm) / $denTwo;
        }

        // 3. Calculate d_phi              ...                                    // and d_lambda
        $dPsi = $thRe;
        $dLambda = $thIm;
        $dPsiN = 1;  // d_psi^0

        $dPhi = 0;
        for ($n = 1; $n <= 9; $n++) {
            $dPsiN = $dPsiN * $dPsi;
            $dPhi = $dPhi + $this->D[$n] * $dPsiN;
        }

        // 4. Calculate latitude and longitude
        // d_phi is calcuated in second of arc * 10^-5, so we need to scale back to radians. d_lambda is in radians.
        $lat = $this->latZero + ($dPhi * ProjFourphp::$common->secToRad * 1E5);
        $lon = $this->longZero + $dLambda;

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['nzmg'] = new ProjFourphp_ProjNzmg();
