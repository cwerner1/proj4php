<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
  NAME                            MOLLWEIDE

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the MOllweide projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  D. Steinwand, EROS      May, 1991;  Updated Sept, 1992; Updated Feb, 1993
  S. Nelson, EDC		Jun, 2993;	Made corrections in precision and
  number of iterations.

  ALGORITHM REFERENCES

  1.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
  Printing Office, Washington D.C., 1989.

  2.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.
 * ***************************************************************************** */

class ProjFourphp_ProjMoll
{
    /* Initialize the Mollweide projection
      ------------------------------------ */

    public function init()
    {
        //no-op
    }

    /* Mollweide forward equations--mapping lat,long to x,y
      ---------------------------------------------------- */

    public function forward($p)
    {

        /* Forward equations
          ----------------- */
        $lon = $p->x;
        $lat = $p->y;

        $deltaLon = ProjFourphp_Common::adjustLon($lon - $this->longZero);
        $theta = $lat;
        $con = ProjFourphp_common::$pi * sin($lat);

        /* Iterate using the Newton-Raphson method to find theta
          ----------------------------------------------------- */
        for ($i = 0; true; ++$i) {
            $deltaTheta = -($theta + sin($theta) - $con) / (1.0 + cos($theta));
            $theta += $deltaTheta;
            if (abs($deltaTheta) < ProjFourphp_Common::$epsln)
                break;
            if ($i >= 50) {
                ProjFourphp::reportError("moll:Fwd:IterationError");
                //return(241);
            }
        }
        $theta /= 2.0;

        /* If the latitude is 90 deg, force the x coordinate to be "0 . false easting"
          this is done here because of precision problems with "cos(theta)"
          -------------------------------------------------------------------------- */
        if (ProjFourphp_common::$pi / 2 - abs($lat) < ProjFourphp_Common::$epsln)
            $deltaLon = 0;
        $x = 0.900316316158 * $this->a * $deltaLon * cos($theta) + $this->xZero;
        $y = 1.4142135623731 * $this->a * sin($theta) + $this->yZero;

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
        #$theta;
        #$arg;

        /* Inverse equations
          ----------------- */
        $p->x-= $this->xZero;
        //~ $p->y -= $this->yZero;
        $arg = $p->y / (1.4142135623731 * $this->a);

        /* Because of division by zero problems, 'arg' can not be 1.0.  Therefore
          a number very close to one is used instead.
          ------------------------------------------------------------------- */
        if (abs($arg) > 0.999999999999)
            $arg = 0.999999999999;
        $theta = asin($arg);
        $lon = ProjFourphp_Common::adjustLon($this->longZero + ($p->x / (0.900316316158 * $this->a * cos($theta))));
        if ($lon < (-ProjFourphp_common::$pi))
            $lon = -ProjFourphp_common::$pi;
        if ($lon > ProjFourphp_common::$pi)
            $lon = ProjFourphp_common::$pi;
        $arg = (2.0 * $theta + sin(2.0 * $theta)) / ProjFourphp_common::$pi;
        if (abs($arg) > 1.0)
            $arg = 1.0;
        $lat = asin($arg);
        //return(OK);

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['moll'] = new ProjFourphp_ProjMoll();