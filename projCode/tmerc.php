<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by Proj4php from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
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
 * ***************************************************************************** */

/**
  Initialize Transverse Mercator projection
 */
class Proj4php_ProjTmerc
{

    private $eZero, $e1, $e2, $e3, $mlZero;

    /**
     * 
     */
    public function init()
    {

        $this->eZero = Proj4php::$common->eZerofn($this->es);
        $this->e1 = Proj4php::$common->e1fn($this->es);
        $this->e2 = Proj4php::$common->e2fn($this->es);
        $this->e3 = Proj4php::$common->e3fn($this->es);
        $this->mlZero = $this->a * Proj4php::$common->mlfn($this->eZero, $this->e1, $this->e2, $this->e3, $this->latZero);
    }

    /**
      Transverse Mercator Forward  - long/lat to x/y
      long/lat in radians
     */
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $deltaLon = Proj4php_Common::adjustLon($lon - $this->longZero); // Delta longitude
        #$con = 0;    // cone constant
        #$x = 0;
        #$y = 0;
        $sinPhi = sin($lat);
        $cosPhi = cos($lat);

        if (isset($this->sphere) && $this->sphere === true) { /* spherical form */
            $b = $cosPhi * sin($deltaLon);
            if ((abs(abs($b) - 1.0)) < .0000000001) {
                Proj4php::reportError("tmerc:forward: Point projects into infinity");
                return(93);
            } else {
                $x = .5 * $this->a * $this->kZero * log((1.0 + $b) / (1.0 - $b));
                $con = acos($cosPhi * cos($deltaLon) / sqrt(1.0 - $b * $b));
                if ($lat < 0)
                    $con = - $con;
                $y = $this->a * $this->kZero * ($con - $this->latZero);
            }
        } else {
            $al = $cosPhi * $deltaLon;
            $als = pow($al, 2);
            $c = $this->ep2 * pow($cosPhi, 2);
            $tq = tan($lat);
            $t = pow($tq, 2);
            $con = 1.0 - $this->es * pow($sinPhi, 2);
            $n = $this->a / sqrt($con);

            $ml = $this->a * Proj4php::$common->mlfn($this->eZero, $this->e1, $this->e2, $this->e3, $lat);

            $x = $this->kZero * $n * $al * (1.0 + $als / 6.0 * (1.0 - $t + $c + $als / 20.0 * (5.0 - 18.0 * $t + pow($t, 2) + 72.0 * $c - 58.0 * $this->ep2))) + $this->xZero;
            $y = $this->kZero * ($ml - $this->mlZero + $n * $tq * ($als * (0.5 + $als / 24.0 * (5.0 - $t + 9.0 * $c + 4.0 * pow($c, 2) + $als / 30.0 * (61.0 - 58.0 * $t + pow($t, 2) + 600.0 * $c - 330.0 * $this->ep2))))) + $this->yZero;
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

        if (isset($this->sphere) && $this->sphere === true) { /* spherical form */
            $f = exp($p->x / ($this->a * $this->kZero));
            $g = .5 * ($f - 1 / $f);
            $temp = $this->latZero + $p->y / ($this->a * $this->kZero);
            $h = cos($temp);
            $con = sqrt((1.0 - $h * $h) / (1.0 + $g * $g));
            $lat = Proj4php_Common::asinz($con);
            if ($temp < 0)
                $lat = -$lat;
            if (($g == 0) && ($h == 0)) {
                $lon = $this->longZero;
            } else {
                $lon = Proj4php_Common::adjustLon(atan2($g, $h) + $this->longZero);
            }
        } else {    // ellipsoidal form
            $x = $p->x - $this->xZero;
            $y = $p->y - $this->yZero;

            $con = ($this->mlZero + $y / $this->kZero) / $this->a;
            $phi = $con;

            for ($i = 0; true; $i++) {
                $deltaPhi = (($con + $this->e1 * sin(2.0 * $phi) - $this->e2 * sin(4.0 * $phi) + $this->e3 * sin(6.0 * $phi)) / $this->eZero) - $phi;
                $phi += $deltaPhi;
                if (abs($deltaPhi) <= Proj4php_Common::$epsln)
                    break;
                if ($i >= $maxIter) {
                    Proj4php::reportError("tmerc:inverse: Latitude failed to converge");
                    return(95);
                }
            } // for()
            if (abs($phi) < Proj4php_Common::$halfPi) {
                // sincos(phi, &sin_phi, &cos_phi);
                $sinPhi = sin($phi);
                $cosPhi = cos($phi);
                $tanPhi = tan($phi);
                $c = $this->ep2 * pow($cosPhi, 2);
                $cs = pow($c, 2);
                $t = pow($tanPhi, 2);
                $ts = pow($t, 2);
                $con = 1.0 - $this->es * pow($sinPhi, 2);
                $n = $this->a / sqrt($con);
                $r = $n * (1.0 - $this->es) / $con;
                $d = $x / ($n * $this->kZero);
                $ds = pow($d, 2);
                $lat = $phi - ($n * $tanPhi * $ds / $r) * (0.5 - $ds / 24.0 * (5.0 + 3.0 * $t + 10.0 * $c - 4.0 * $cs - 9.0 * $this->ep2 - $ds / 30.0 * (61.0 + 90.0 * $t + 298.0 * $c + 45.0 * $ts - 252.0 * $this->ep2 - 3.0 * $cs)));
                $lon = Proj4php_Common::adjustLon($this->longZero + ($d * (1.0 - $ds / 6.0 * (1.0 + 2.0 * $t + $c - $ds / 20.0 * (5.0 - 2.0 * $c + 28.0 * $t - 3.0 * $cs + 8.0 * $this->ep2 + 24.0 * $ts))) / $cosPhi));
            } else {
                $lat = Proj4php_Common::$halfPi * Proj4php::$common->sign($y);
                $lon = $this->longZero;
            }
        }

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

Proj4php::$proj['tmerc'] = new Proj4php_ProjTmerc();