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
  NAME                             EQUIRECTANGULAR

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Equirectangular projection.  The
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
 * ************************************************************************** */
class ProjFourphp_ProjEqui
{

    public function init()
    {
        if (!$this->xZero) $this->xZero    = 0;
        if (!$this->yZero) $this->yZero    = 0;
        if (!$this->latZero) $this->latZero  = 0;
        if (!$this->longZero) $this->longZero = 0;
        ///$this->tTwo;
    }

    /* Equirectangular forward equations--mapping lat,long to x,y
      --------------------------------------------------------- */

    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $dlon = ProjFourphp_Common::adjustLon($lon - $this->longZero);
        $x    = $this->xZero + $this->a * $dlon * cos($this->latZero);
        $y    = $this->yZero + $this->a * $lat;

        $this->tOne = $x;
        $this->tTwo = cos($this->latZero);
        $p->x       = $x;
        $p->y       = $y;
        return $p;
    }

    /* Equirectangular inverse equations--mapping x,y to lat/long
      --------------------------------------------------------- */

    public function inverse($p)
    {

        $p->x -= $this->xZero;
        $p->y -= $this->yZero;
        $lat = $p->y / $this->a;

        if (abs($lat) > ProjFourphp_Common::$halfPi) {
            ProjFourphp::reportError("equi:Inv:DataError");
        }
        $lon  =
            ProjFourphp_Common::adjustLon($this->longZero + $p->x / ($this->a * cos($this->latZero)));
        $p->x = $lon;
        $p->y = $lat;
    }

}

ProjFourphp::$proj['equi'] = new ProjFourphp_ProjEqui();