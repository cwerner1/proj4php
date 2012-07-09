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
/* * *************************************************************************
  NAME                    VAN DER GRINTEN

  PURPOSE:	Transforms input Easting and Northing to longitude and
  latitude for the Van der Grinten projection.  The
  Easting and Northing must be in meters.  The longitude
  and latitude values will be returned in radians.

  PROGRAMMER              DATE
  ----------              ----
  T. Mittan		March, 1993

  This function was adapted from the Van Der Grinten projection code
  (FORTRAN) in the General Cartographic Transformation Package software
  which is available from the U.S. Geological Survey National Mapping Division.

  ALGORITHM REFERENCES

  1.  "New Equal-Area Map Projections for Noncircular Regions", John P. Snyder,
  The American Cartographer, Vol 15, No. 4, October 1988, pp. 341-355.

  2.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  3.  "Software Documentation for GCTP General Cartographic Transformation
  Package", U.S. Geological Survey National Mapping Division, May 1982.
 * ****************************************************************************
 */

class ProjFourphp_ProjVandg
{
    /* Initialize the Van Der Grinten projection
      ---------------------------------------- */

    /**
     * Radius of earth
     * @var type 
     */
    public static $radius = 6370997.0;

    public function init()
    {
        
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        /* Forward equations
          ----------------- */
        $dlon = ProjFourphp_Common::adjustLon($lon - $this->longZero);
        $x;
        $y;

        if (abs($lat) <= ProjFourphp_Common::$epsln) {
            $x     = $this->xZero + static::$radius * $dlon;
            $y     = $this->yZero;
        }
        $theta = ProjFourphp_common::asinz(2.0 * abs($lat / ProjFourphp_common::$pi));
        if ((abs($dlon) <= ProjFourphp_Common::$epsln) ||
            (abs(abs($lat) - ProjFourphp_Common::$halfPi) <= ProjFourphp_Common::$epsln)) {
            $x = $this->xZero;
            if ($lat >= 0) {
                $y = $this->yZero + ProjFourphp_common::$pi *
                    static::$radius * tan(.5 * $theta);
            } else {
                $y     = $this->yZero + ProjFourphp_common::$pi *
                    static::$radius * - tan(.5 * $theta);
            }
            //  return(OK);
        }
        $al    = .5 * abs((ProjFourphp_common::$pi / $dlon) - ($dlon / ProjFourphp_common::$pi));
        $asq   = $al * $al;
        $sinth = sin($theta);
        $costh = cos($theta);

        $g   = $costh / ($sinth + $costh - 1.0);
        $gsq = $g * $g;
        $m   = $g * (2.0 / $sinth - 1.0);
        $msq = $m * $m;
        $con = ProjFourphp_common::$pi * static::$radius * ($al * ($g - $msq) +
            sqrt($asq * ($g - $sq) * ($g - $msq) - ($msq + $asq) * ($gsq - $msq))
            ) / ($msq + $asq);
        if ($dlon < 0) {
            $con = -$con;
        }
        $x   = $this->xZero + $con;
        $con = abs($con / (ProjFourphp_common::$pi * static::$radius));
        if ($lat >= 0) {
            $y = $this->yZero + ProjFourphp_common::$pi * static::$radius *
                sqrt(1.0 - $con * $con - 2.0 * $al * $con);
        } else {
            $y = $this->yZero - ProjFourphp_common::$pi * static::$radius *
                sqrt(1.0 - $con * $con - 2.0 * $al * $con);
        }

        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /* Van Der Grinten inverse equations--mapping x,y to lat/long
      --------------------------------------------------------- */

    public function inverse($p)
    {

        /*
          $dlon;
          $xx;
          $yy;
          $xys;
          $c1;
          $c2;
          $c3;
          $al;
          $asq;
          $aOne;
          $m1;
          $con;
          $th1;
          $d;
         */

        /* inverse equations
          ----------------- */
        $p->x -= $this->xZero;
        $p->y -= $this->yZero;
        $con    = ProjFourphp_common::$pi * static::$radius;
        $xx     = $p->x / $con;
        $yy     = $p->y / $con;
        $xys    = $xx * $xx + $yy * $yy;
        $cOne   = -abs($yy) * (1.0 + $xys);
        $cTwo   = $cOne - 2.0 * $yy * $yy + $xx * $xx;
        $cThree = -2.0 * $cOne + 1.0 + 2.0 * $yy * $yy + $xys * $xys;
        $d      = $yy * $yy / $cThree +
            (2.0 * $cTwo * $cTwo * $cTwo / $cThree / $cThree / $cThree - 9.0 *
            $cOne * $cTwo / $cThree / $cThree) / 27.0;
        $aOne   = ($cOne - $cTwo * $cTwo / 3.0 / $cThree) / $cThree;
        $mOne   = 2.0 * sqrt(-$aOne / 3.0);
        $con    = ((3.0 * $d) / $aOne) / $mOne;
        if (abs($con) > 1.0) {
            if ($con >= 0.0) {
                $con = 1.0;
            } else {
                $con   = -1.0;
            }
        }
        $thOne = acos($con) / 3.0;
        if ($p->$y >= 0) {
            $lat = (-$mOne * cos($thOne + ProjFourphp_common::$pi / 3.0) -
                $cTwo / 3.0 / $cThree) * ProjFourphp_common::$pi;
        } else {
            $lat = -(-m1 * cos($thOne + ProjFourphp_common::$pi / 3.0) -
                $cTwo / 3.0 / $cThree) * ProjFourphp_common::$pi;
        }

        if (abs($xx) < ProjFourphp_Common::$epsln) {
            $lon    = $this->longZero;
        }
        $adjLon = $this->longZero + ProjFourphp_common::$pi *
            ($xys - 1.0 +
            sqrt(1.0 + 2.0 * ($xx * $xx - $yy * $yy) + $xys * $xys))
            / 2.0 / $xx;
        $lon    = ProjFourphp_Common::adjustLon($adjLon);

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['vandg'] = new ProjFourphp_ProjVandg();