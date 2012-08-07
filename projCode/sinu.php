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
/* * ***************************************************************************
  NAME                  		SINUSOIDAL

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Sinusoidal projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  D. Steinwand, EROS      May, 1991

  This function was adapted from the Sinusoidal projection code (FORTRAN) in the
  General Cartographic Transformation Package software which is available from
  the U.S. Geological Survey National Mapping Division.

  ALGORITHM REFERENCES

  1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  "Software Documentation for GCTP General Cartographic Transformation
  Package", U.S. Geological Survey National Mapping Division, May 1982.
 * ************************************************************************** */

class ProjFourphp_ProjSinu
{
    /* Initialize the Sinusoidal projection
      ------------------------------------ */

    public function init()
    {
        /* Place parameters in static storage for common use
          ------------------------------------------------- */
        #$this->R = 6370997.0; //Radius of earth

        if (!$this->sphere) {
            $this->en = ProjFourphp::$common->pjEnfn($this->es);
        } else {
            $this->n  = 1.;
            $this->m  = 0.;
            $this->es = 0;
            $this->cY = sqrt(($this->m + 1.) / $this->n);
            $this->cX = $this->cY / ($this->m + 1.);
        }
    }

    /* Sinusoidal forward equations--mapping lat,long to x,y
      ----------------------------------------------------- */

    public function forward($p)
    {

        #$x,y,delta_lon;
        $lon = $p->x;
        $lat = $p->y;

        /* Forward equations
          ----------------- */
        $lon = ProjFourphp_Common::adjustLon($lon - $this->longZero);

        if (isset($this->sphere)) {
            if (!$this->m) {
                $lat = $this->n != 1. ? asin($this->n * sin($lat)) : $lat;
            } else {
                $k = $this->n * sin($lat);
                for ($i = ProjFourphp_Common::$maxIter; $i; --$i) {
                    $v = ($this->m * $lat + sin($lat) - $k) /
                        ($this->m + cos($lat));
                    $lat -= $v;
                    if (abs($v) < ProjFourphp_Common::$epsln) break;
                }
            }
            $x = $this->a * $this->cX * $lon * ($this->m + cos($lat));
            $y = $this->a * $this->cY * $lat;
        } else {

            $s = sin($lat);
            $c = cos($lat);
            $y = $this->a * ProjFourphp_Common::jMlfn($lat, $s, $c, $this->en);
            $x = $this->a * $lon * $c / sqrt(1. - $this->es * $s * $s);
        }

        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function inverse($p)
    {
        #$lat;
        #$temp;
        #$lon;

        /* Inverse equations
          ----------------- */
        $p->x -= $this->xZero;
        $p->y -= $this->yZero;
        $lat = $p->y / $this->a;

        if (isset($this->sphere)) {

            $p->y /= $this->cY;
            if ($this->m) {
                $lat = asin(($this->m * $p->y + sin($p->y)) / $this->n);
            } else {
                $lat = $p->y;
            }
            $lon = $p->x / ($this->cX * ($this->m + cos($p->y)));
        } else {
            $lat =
                ProjFourphp::$common->pjInvMlfn($p->y / $this->a, $this->es, $this->en);
            $s   = abs($lat);

            if ($s < ProjFourphp_Common::$halfPi) {
                $s    = sin($lat);
                $temp =
                    $this->longZero + $p->x * sqrt(1. - $this->es * $s * $s) /
                    ($this->a * cos($lat));
                //temp = $this->longZero + $p->x / ($this->a * cos($lat));
                $lon  = ProjFourphp_Common::adjustLon($temp);
            } else if (($s - ProjFourphp_Common::$epsln)
                < ProjFourphp_Common::$halfPi) {
                $lon = $this->longZero;
            }
        }

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['sinu'] = new ProjFourphp_ProjSinu();