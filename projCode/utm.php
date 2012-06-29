<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
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
class ProjFourphp_ProjUtm
{

    public $dependsOn = 'tmerc';
    public $utmSouth = false; // UTM north/south

    /**
     *
     * @return void 
     */

    public function init()
    {

        if (!isset($this->zone)) {
            ProjFourphp::reportError("utm:init: zone must be specified for UTM");
            return;
        }

        $this->latZero = 0.0;
        $this->longZero = ((6 * abs($this->zone)) - 183) * ProjFourphp::$common->dToR;
        $this->xZero = 500000.0;
        $this->yZero = $this->utmSouth ? 10000000.0 : 0.0;
        $this->kZero = 0.9996;
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function forward($p)
    {
        return ProjFourphp::$proj['tmerc']->forward($p);
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function inverse($p)
    {
        return ProjFourphp::$proj['tmerc']->inverse($p);
    }

}

ProjFourphp::$proj['utm'] = new ProjFourphp_ProjUtm();