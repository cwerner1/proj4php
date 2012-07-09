<?php
/**
 * @package Proj4
 */

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * **************************************************************************
  NAME                            TRANSVERSE MERCATOR

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Transverse Mercator projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  ALGORITHM REFERENCES

  1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
  Printing Office, Washington D.C., 1989.
 * ****************************************************************************
 *  */

/**
  Initialize Transverse Mercator projection
 */
class ProjFourphp_ProjTmerc
{

    private $_eZero, $_eOne, $_eTwo, $_eThree, $_mlZero;

    /**
     * 
     */
    public function init()
    {

        $this->_eZero = ProjFourphp_Common::eZeroFn($this->es);
        $this->_eOne = ProjFourphp_Common::eOnefn($this->es);
        $this->_eTwo = ProjFourphp_Common::eTwofn($this->es);
        $this->_eThree = ProjFourphp_Common::eThreefn($this->es);
        $this->_mlZero = $this->a *
            ProjFourphp_Common::mlfn($this->_eZero, $this->_eOne, $this->_eTwo, $this->_eThree, $this->latZero);
    }

    /**
      Transverse Mercator Forward  - long/lat to x/y
      long/lat in radians
     */
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $deltaLon = ProjFourphp_Common::adjustLon($lon - $this->longZero);
        //// Delta longitude
        #$con = 0;    // cone constant
        #$x = 0;
        #$y = 0;
        $sinPhi   = sin($lat);
        $cosPhi   = cos($lat);

        if (isset($this->sphere) && $this->sphere === true) {
            /* spherical form */
            $b = $cosPhi * sin($deltaLon);
            if ((abs(abs($b) - 1.0)) < .0000000001) {
                ProjFourphp::reportError("tmerc:forward: Point projects into infinity");
                return(93);
            } else {
                $x   = .5 * $this->a * $this->kZero *
                    log((1.0 + $b) / (1.0 - $b));
                $con = acos($cosPhi * cos($deltaLon) / sqrt(1.0 - $b * $b));
                if ($lat < 0) $con = - $con;
                $y   = $this->a * $this->kZero * ($con - $this->latZero);
            }
        } else {
            $al  = $cosPhi * $deltaLon;
            $als = pow($al, 2);
            $c   = $this->epTwo * pow($cosPhi, 2);
            $tq  = tan($lat);
            $t   = pow($tq, 2);
            $con = 1.0 - $this->es * pow($sinPhi, 2);
            $n   = $this->a / sqrt($con);

            $ml = $this->a *
                ProjFourphp_Common::mlfn($this->_eZero, $this->_eOne, $this->_eTwo, $this->_eThree, $lat);

            $x = $this->kZero * $n * $al *
                (1.0 + $als / 6.0 *
                (1.0 - $t + $c + $als / 20.0 *
                (5.0 - 18.0 * $t + pow($t, 2) + 72.0 * $c
                - 58.0 * $this->epTwo))) +
                $this->xZero;
            $y = $this->kZero *
                ($ml - $this->_mlZero + $n * $tq *
                ($als * (0.5 + $als / 24.0 *
                (5.0 - $t + 9.0 * $c + 4.0 * pow($c, 2)
                + $als / 30.0 *
                (61.0 - 58.0 * $t + pow($t, 2) + 600.0 * $c - 330.0 *
                $this->epTwo))))) + $this->yZero;
        }

        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /**
      Transverse Mercator Inverse  -  x/y to long/lat
     */
    public function inverse($p)
    {

        #$phi;  /* temporary angles       */
        #$delta_phi; /* difference between longitudes    */
        $maxIter = 6;      /* maximun number of iterations */

        if (isset($this->sphere) && $this->sphere === true) {
            /* spherical form */
            $f    = exp($p->x / ($this->a * $this->kZero));
            $g    = .5 * ($f - 1 / $f);
            $temp = $this->latZero + $p->y / ($this->a * $this->kZero);
            $h    = cos($temp);
            $con  = sqrt((1.0 - $h * $h) / (1.0 + $g * $g));
            $lat  = ProjFourphp_Common::asinz($con);
            if ($temp < 0) $lat  = -$lat;
            if (($g == 0) && ($h == 0)) {
                $lon = $this->longZero;
            } else {
                $lon = ProjFourphp_Common::adjustLon(atan2($g, $h) + $this->longZero);
            }
        } else {    // ellipsoidal form
            $x = $p->x - $this->xZero;
            $y = $p->y - $this->yZero;

            $con = ($this->_mlZero + $y / $this->kZero) / $this->a;
            $phi = $con;

            for ($i = 0; true; $i++) {
                $deltaPhi = (($con + $this->_eOne * sin(2.0 * $phi) -
                    $this->_eTwo * sin(4.0 * $phi) + $this->_eThree *
                    sin(6.0 * $phi))
                    / $this->_eZero) - $phi;
                $phi += $deltaPhi;
                if (abs($deltaPhi) <= ProjFourphp_Common::$epsln) break;
                if ($i >= $maxIter) {
                    ProjFourphp::reportError("tmerc:inverse: Latitude failed to converge");
                    return(95);
                }
            } // for()
            if (abs($phi) < ProjFourphp_Common::$halfPi) {
                // sincos(phi, &sin_phi, &cos_phi);
                $sinPhi = sin($phi);
                $cosPhi = cos($phi);
                $tanPhi = tan($phi);
                $c      = $this->epTwo * pow($cosPhi, 2);
                $cs     = pow($c, 2);
                $t      = pow($tanPhi, 2);
                $ts     = pow($t, 2);
                $con    = 1.0 - $this->es * pow($sinPhi, 2);
                $n      = $this->a / sqrt($con);
                $r      = $n * (1.0 - $this->es) / $con;
                $d      = $x / ($n * $this->kZero);
                $ds     = pow($d, 2);
                $lat    = $phi -
                    ($n * $tanPhi * $ds / $r) *
                    (0.5 - $ds / 24.0 *
                    (5.0 + 3.0 * $t + 10.0 * $c - 4.0 * $cs -
                    9.0 * $this->epTwo - $ds / 30.0 *
                    (61.0 + 90.0 * $t + 298.0 * $c + 45.0 * $ts - 252.0 *
                    $this->epTwo - 3.0 * $cs)));
                $lon    = ProjFourphp_Common::adjustLon($this->longZero + ($d * (1.0 - $ds / 6.0 * (1.0 + 2.0 * $t + $c - $ds / 20.0 * (5.0 - 2.0 * $c + 28.0 * $t - 3.0 * $cs + 8.0 * $this->epTwo + 24.0 * $ts))) / $cosPhi));
            } else {
                $lat = ProjFourphp_Common::$halfPi * ProjFourphp::$common->sign($y);
                $lon = $this->longZero;
            }
        }

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['tmerc'] = new ProjFourphp_ProjTmerc();