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
  NAME                            EQUIDISTANT CONIC

  PURPOSE:	Transforms input longitude and latitude to Easting and Northing
  for the Equidistant Conic projection.  The longitude and
  latitude must be in radians.  The Easting and Northing values
  will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  T. Mittan		Mar, 1993

  ALGORITHM REFERENCES

  1.  Snyder, John $p->, "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  Snyder, John $p-> and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
  Printing Office, Washington D.C., 1989.
 * ************************************************************************** */

/* Variables common to all subroutines in this code file
  ----------------------------------------------------- */

class ProjFourphp_ProjEqdc
{
    /* Initialize the Equidistant Conic projection
      ------------------------------------------ */

    public function init()
    {

        /* Place parameters in static storage for common use
          ------------------------------------------------- */
        if (!$this->mode) $this->mode   = 0; //chosen default mode
        $this->temp   = $this->b / $this->a;
        $this->es     = 1.0 - pow($this->temp, 2);
        $this->e      = sqrt($this->es);
        $this->eZero  = ProjFourphp_Common::eZeroFn($this->es);
        $this->eOne   = ProjFourphp_Common::eOnefn($this->es);
        $this->eTwo   = ProjFourphp_Common::eTwofn($this->es);
        $this->eThree = ProjFourphp_Common::eThreefn($this->es);

        $this->sinphi = sin($this->latOne);
        $this->cosphi = cos($this->latOne);

        $this->msOne =
            ProjFourphp_Common::msfnz($this->e, $this->sinphi, $this->cosphi);
        $this->mlOne =
            ProjFourphp_Common::mlfn($this->eZero, $this->eOne, $this->eTwo, $this->eThree, $this->latOne);

        /* format B
          --------- */
        if ($this->mode != 0) {
            if (abs($this->latOne + $this->latTwo) < ProjFourphp_Common::$epsln) {
                ProjFourphp::reportError("eqdc:Init:EqualLatitudes");
                //return(81);
            }
            $this->sinphi = sin($this->latTwo);
            $this->cosphi = cos($this->latTwo);

            $this->msTwo =
                ProjFourphp_Common::msfnz($this->e, $this->sinphi, $this->cosphi);
            $this->mlTwo =
                ProjFourphp_Common::mlfn($this->eZero, $this->eOne, $this->eTwo, $this->eThree, $this->latTwo);
            if (abs($this->latOne - $this->latTwo)
                >= ProjFourphp_Common::$epsln) {
                $this->ns =
                    ($this->msOne - $this->msTwo)
                    / ($this->mlTwo - $this->mlOne);
            } else {
                $this->ns = $this->sinphi;
            }
        } else {
            $this->ns     = $this->sinphi;
        }
        $this->g      = $this->mlOne + $this->msOne / $this->ns;
        $this->mlZero =
            ProjFourphp_Common::mlfn($this->eZero, $this->eOne, $this->eTwo, $this->eThree, $this->latZero);
        $this->rh     = $this->a * ($this->g - $this->mlZero);
    }

    /* Equidistant Conic forward equations--mapping lat,long to x,y
      ----------------------------------------------------------- */

    public function forward($p)
    {
        $lon = $p->x;
        $lat = $p->y;

        /* Forward equations
          ----------------- */
        $ml    =
            ProjFourphp_Common::mlfn($this->eZero, $this->eOne, $this->eTwo, $this->eThree, $lat);
        $rhOne = $this->a * ($this->g - $ml);
        $theta =
            $this->ns * ProjFourphp_Common::adjustLon($lon - $this->longZero);

        $x    = $this->xZero + $rhOne * sin($theta);
        $y    = $this->yZero + $this->rh - $rhOne * cos($theta);
        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /* Inverse equations
      ----------------- */

    public function inverse($p)
    {

        $p->x -= $this->xZero;
        $p->y = $this->rh - $p->y + $this->yZero;

        if ($this->ns >= 0) {
            $rhOne = sqrt($p->x * $p->x + $p->y * $p->y);
            $con   = 1.0;
        } else {
            $rhOne = -sqrt($p->x * $p->x + $p->y * $p->y);
            $con   = -1.0;
        }
        $theta = 0.0;
        if ($rhOne != 0.0) $theta = atan2($con * $p->x, $con * $p->y);
        $ml    = $this->g - $rhOne / $this->a;
        $lat   =
            static::phi3z($ml, $this->eZero, $this->eOne, $this->eTwo, $this->eThree);
        $lon   =
            ProjFourphp_Common::adjustLon($this->longZero + $theta / $this->ns);

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

    /* Function to compute latitude, phi3, for the inverse of the Equidistant
      Conic projection.
      ----------------------------------------------------------------- */

    public static function phi3z($ml, $eZero, $eOne, $eTwo, $eThree)
    {

        $phi = $ml;
        for ($i = 0; $i < 15; $i++) {
            $dphi = ($ml + $eOne * sin(2.0 * $phi) - $eTwo * sin(4.0 * $phi) + $eThree * sin(6.0 * $phi)) / $eZero - $phi;
            $phi += $dphi;
            if (abs($dphi) <= .0000000001) {
                return $phi;
            }
        }

        ProjFourphp::reportError("PHI3Z-CONV:Latitude failed to converge after 15 iterations");

        return null;
    }

}

ProjFourphp::$proj['eqdc'] = new ProjFourphp_ProjEqdc();