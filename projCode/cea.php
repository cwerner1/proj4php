<?php

/* * *****************************************************************************
  NAME                    LAMBERT CYLINDRICAL EQUAL AREA

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Lambert Cylindrical Equal Area projection.
  This class of projection includes the Behrmann and
  Gall-Peters Projections.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  R. Marsden              August 2009
  Winwaed Software Tech LLC, http://www.winwaed.com

  This function was adapted from the Miller Cylindrical Projection in the ProjFourphp
  library.

  Note: This implementation assumes a Spherical Earth. The (commented) code
  has been included for the ellipsoidal forward transform, but derivation of
  the ellispoidal inverse transform is beyond me. Note that most of the
  ProjFourphp implementations do NOT currently support ellipsoidal figures.
  Therefore this is not seen as a problem - especially this lack of support
  is explicitly stated here.

  ALGORITHM REFERENCES

  1.  "Cartographic Projection Procedures for the UNIX Environment -
  A User's Manual" by Gerald I. Evenden, USGS Open File Report 90-284
  and Release 4 Interim Reports (2003)

  2.  Snyder, John P., "Flattening the Earth - Two Thousand Years of Map
  Projections", Univ. Chicago Press, 1993
 * ***************************************************************************** */

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
class ProjFourphp_ProjCea
{
    /* Initialize the Cylindrical Equal Area projection
      ------------------------------------------- */

    public function init()
    {
        //no-op
    }

    /* Cylindrical Equal Area forward equations--mapping lat,long to x,y
      ------------------------------------------------------------ */

    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        /* Forward equations
          ----------------- */
        $dlon =  ProjFourphp_Common::adjustLon($lon - $this->longZero);
        $x = $this->xZero + $this->a * $dlon * cos($this->latTs);
        $y = $this->yZero + $this->a * sin($lat) / cos($this->latTs);
        /* Elliptical Forward Transform
          Not implemented due to a lack of a matchign inverse function
          {
          $Sin_Lat = sin(lat);
          $Rn = $this->a * (sqrt(1.0e0 - $this->es * Sin_Lat * Sin_Lat ));
          x = $this->xZero + $this->a * dlon * cos($this->latTs);
          y = $this->yZero + Rn * sin(lat) / cos($this->latTs);
          }
         */

        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /**
     * Cylindrical Equal Area inverse equations--mapping x,y to lat/long
     * 
     * @param type $p
     * @return type 
     */
    public function inverse($p)
    {
        $p->x -= $this->xZero;
        $p->y -= $this->yZero;

        $p->x =  ProjFourphp_Common::adjustLon($this->longZero + ($p->x / $this->a) / cos($this->latTs));
        $p->y = asin(($p->y / $this->a) * cos($this->latTs));

        return $p;
    }

}

ProjFourphp::$proj['cea'] = new ProjFourphp_ProjCea();