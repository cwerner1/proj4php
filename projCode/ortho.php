<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
  NAME                             ORTHOGRAPHIC

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Orthographic projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  T. Mittan		Mar, 1993

  ALGORITHM REFERENCES

  1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
  Printing Office, Washington D.C., 1989.
 * ***************************************************************************** */

class ProjFourphp_ProjOrtho
{
    /* Initialize the Orthographic projection
      ------------------------------------- */

    public function init($def)
    {
        //double temp;			/* temporary variable		*/

        /* Place parameters in static storage for common use
          ------------------------------------------------- */;
        $this->sinPOneFour = sin($this->latZero);
        $this->cosPOneFour = cos($this->latZero);
    }

    /* Orthographic forward equations--mapping lat,long to x,y
      --------------------------------------------------- */

    public function forward($p)
    {

        /*
          $sinphi;
          $cosphi; // sin and cos value
          $dlon;  // delta longitude value
          $coslon;  // cos of longitude
          $ksp;  // scale factor
          $g;
         */

        $lon = $p->x;
        $lat = $p->y;

        /* Forward equations
          ----------------- */
        $dlon = ProjFourphp_Common::adjustLon($lon - $this->longZero);

        $sinphi = sin($lat);
        $cosphi = cos($lat);

        $coslon = cos($dlon);
        $g = $this->sinPOneFour * sinphi + $this->cosPOneFour * $cosphi * $coslon;
        $ksp = 1.0;

        if (($g > 0) || (abs($g) <= ProjFourphp_Common::$epsln)) {
            $x = $this->a * $ksp * $cosphi * sin($dlon);
            $y = $this->yZero + $this->a * $ksp * ($this->cosPOneFour * $sinphi - $this->sinPOneFour * $cosphi * $coslon);
        } else {
            ProjFourphp::reportError("orthoFwdPointError");
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

        /*
          $rh;  // height above ellipsoid
          $z;  // angle
          $sinz;
          $cosz; // sin of z and cos of z
          $temp;
          $con;
          $lon;
          $lat;
         */

        /* Inverse equations
          ----------------- */
        $p->x -= $this->xZero;
        $p->y -= $this->yZero;
        $rh = sqrt($p->x * $p->x + $p->y * $p->y);
        if ($rh > $this->a + .0000001) {
            ProjFourphp::reportError("orthoInvDataError");
        }
        $z = ProjFourphp::$common . asinz($rh / $this->a);

        $sinz = sin($z);
        $cosz = cos($z);

        $lon = $this->longZero;
        if (abs($rh) <= ProjFourphp_Common::$epsln) {
            $lat = $this->latZero;
        }
        $lat = ProjFourphp::$common . asinz($cosz * $this->sinPOneFour + ($p->y * $sinz * $this->cosPOneFour) / $rh);
        $con = abs($this->latZero) - ProjFourphp_Common::$halfPi;
        if (abs(con) <= ProjFourphp_Common::$epsln) {
            if ($this->latZero >= 0) {
                $lon = ProjFourphp_Common::adjustLon($this->longZero + atan2($p->x, -$p->y));
            } else {
                $lon = ProjFourphp_Common::adjustLon($this->longZero - atan2(-$p->x, $p->y));
            }
        }
        $con = $cosz - $this->sinPOneFour * sin($lat);

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['ortho'] = new ProjFourphp_ProjOrtho();