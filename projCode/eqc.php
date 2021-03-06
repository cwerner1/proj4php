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
/* similar to equi.js FIXME projFour uses eqc */
class ProjFourphp_ProjEqc
{

    public function init()
    {

        if (!$this->xZero) $this->xZero    = 0;
        if (!$this->yZero) $this->yZero    = 0;
        if (!$this->latZero) $this->latZero  = 0;
        if (!$this->longZero) $this->longZero = 0;
        if (!$this->latTs) {
            $this->latTs = 0;
        }
        if (!$this->title) {
            $this->title = "Equidistant Cylindrical (Plate Carre)";
        }

        $this->rc = cos($this->latTs);
    }

    // forward equations--mapping lat,long to x,y
    // -----------------------------------------------------------------
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $dlon = ProjFourphp_Common::adjustLon($lon - $this->longZero);
        $dlat = ProjFourphp_Common::adjust_lat($lat - $this->latZero);
        $p->x = $this->xZero + ($this->a * $dlon * $this->rc);
        $p->y = $this->yZero + ($this->a * $dlat );
        return $p;
    }

    // inverse equations--mapping x,y to lat/long
    // -----------------------------------------------------------------
    public function inverse($p)
    {

        $x = $p->x;
        $y = $p->y;

        $p->x =
            ProjFourphp_Common::adjustLon($this->longZero + (($x - $this->xZero) / ($this->a * $this->rc)));
        $p->y =
            ProjFourphp_Common::adjustLat($this->latZero + (($y - $this->yZero) / ($this->a )));
        return $p;
    }

}

ProjFourphp::$proj['eqc'] = new ProjFourphp_ProjEqc();