<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
  NAME                    MILLER CYLINDRICAL

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Miller Cylindrical projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  T. Mittan		March, 1993

  This function was adapted from the Lambert Azimuthal Equal Area projection
  code (FORTRAN) in the General Cartographic Transformation Package software
  which is available from the U.S. Geological Survey National Mapping Division.

  ALGORITHM REFERENCES

  1.  "New Equal-Area Map Projections for Noncircular Regions", John P. Snyder,
  The American Cartographer, Vol 15, No. 4, October 1988, pp. 341-355.

  2.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  3.  "Software Documentation for GCTP General Cartographic Transformation
  Package", U.S. Geological Survey National Mapping Division, May 1982.
 * ***************************************************************************** */

class ProjFourphp_ProjMill
{
    /* Initialize the Miller Cylindrical projection
      ------------------------------------------- */

    public function init()
    {
        //no-op
    }

    /* Miller Cylindrical forward equations--mapping lat,long to x,y
      ------------------------------------------------------------ */

    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        /* Forward equations
          ----------------- */
        $dlon = ProjFourphp_Common::adjustLon($lon - $this->longZero);
        $x = $this->xZero + $this->a * $dlon;
        $y = $this->yZero + $this->a * log(tan((ProjFourphp::$common->pi / 4.0) + ($lat / 2.5))) * 1.25;

        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /* Miller Cylindrical inverse equations--mapping x,y to lat/long
      ------------------------------------------------------------ */

    public function inverse($p)
    {

        $p->x -= $this->xZero;
        $p->y -= $this->yZero;

        $lon = ProjFourphp_Common::adjustLon($this->longZero + $p->x / $this->a);
        $lat = 2.5 * (atan(exp(0.8 * $p->y / $this->a)) - ProjFourphp::$common->pi / 4.0);

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['mill'] = new ProjFourphp_ProjMill();