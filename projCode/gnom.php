<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * ***************************************************************************
  NAME                             GNOMONIC

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Gnomonic Projection.
  Implementation based on the existing sterea and ortho
  implementations.

  PROGRAMMER              DATE
  ----------              ----
  Richard Marsden         November 2009

  ALGORITHM REFERENCES

  1.  Snyder, John P., "Flattening the Earth - Two Thousand Years of Map
  Projections", University of Chicago Press 1993

  2.  Wolfram Mathworld "Gnomonic Projection"
  http://mathworld.wolfram.com/GnomonicProjection.html
  Accessed: 12th November 2009
 * **************************************************************************** */

class ProjFourphp_ProjGnom
{

    /**
     * Initialize the Gnomonic projection
     * 
     * @todo $def not used in context...?
     * @param type $def 
     */
    public function init($def)
    {

        /* Place parameters in static storage for common use
          ------------------------------------------------- */
        $this->sinPOneFour = sin($this->latZero);
        $this->cosPOneFour = cos($this->latZero);

        // Approximation for projecting points to the horizon (infinity)
        $this->infinity_dist = 1000 * $this->a;
        $this->rc = 1;
    }

    /* Gnomonic forward equations--mapping lat,long to x,y
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
        $dlon =  ProjFourphp_Common::adjustLon($lon - $this->longZero);

        $sinphi = sin($lat);
        $cosphi = cos($lat);

        $coslon = cos($dlon);
        $g = $this->sinPOneFour * $sinphi + $this->cosPOneFour * $cosphi * $coslon;
        $ksp = 1.0;

        if ((g > 0) || (abs(g) <= ProjFourphp_Common::$epsln)) {
            $x = $this->xZero + $this->a * $ksp * $cosphi * sin($dlon) / $g;
            $y = $this->yZero + $this->a * $ksp * ($this->cosPOneFour * $sinphi - $this->sinPOneFour * $cosphi * $coslon) / $g;
        } else {
            ProjFourphp::reportError("orthoFwdPointError");

            // Point is in the opposing hemisphere and is unprojectable
            // We still need to return a reasonable point, so we project 
            // to infinity, on a bearing 
            // equivalent to the northern hemisphere equivalent
            // This is a reasonable approximation for short shapes and lines that 
            // straddle the horizon.

            $x = $this->xZero + $this->infinity_dist * $cosphi * sin($dlon);
            $y = $this->yZero + $this->infinity_dist * ($this->cosPOneFour * $sinphi - $this->sinPOneFour * $cosphi * $coslon);
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
          $rh;  // Rho
          $z;  // angle
          $sinc;
          $cosc;
          $c;
          $lon;
          $lat;
         */

        /* Inverse equations
          ----------------- */
        $p->x = ($p->x - $this->xZero) / $this->a;
        $p->y = ($p->y - $this->yZero) / $this->a;

        $p->x /= $this->kZero;
        $p->y /= $this->kZero;

        if (($rh = sqrt($p->x * $p->x + $p->y * $p->y))) {
            $c = atanTwo($rh, $this->rc);
            $sinc = sin($c);
            $cosc = cos($c);

            $lat = ProjFourphp_Common::asinz($cosc * $this->sinPOneFour + ($p->y * $sinc * $this->cosPOneFour) / $rh);
            $lon = atanTwo($p->x * sinc, rh * $this->cosPOneFour * $cosc - $p->y * $this->sinPOneFour * $sinc);
            $lon =  ProjFourphp_Common::adjustLon($this->longZero + $lon);
        } else {
            $lat = $this->phic0;
            $lon = 0.0;
        }

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['gnom'] = new ProjFourphp_ProjGnom();