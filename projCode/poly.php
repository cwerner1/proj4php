<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* Function to compute, phi4, the latitude for the inverse of the
  Polyconic projection.
  ------------------------------------------------------------ */
function phi4z($eccent, $eZero, $eOne, $eTwo, $eThree, $a, $b, &$c, $phi)
{
    /*
      $sinphi;
      $sinTwoph;
      $tanph;
      $ml;
      $mlp;
      $conOne;
      $conTwo;
      $con3;
      $dphi;
      $i;
     */

    $phi = $a;
    for ($i   = 1; $i <= 15; $i++) {
        $sinphi   = sin($phi);
        $tanphi   = tan($phi);
        $c        = $tanphi * sqrt(1.0 - $eccent * $sinphi * $sinphi);
        $sinTwoph = sin(2.0 * $phi);
        /*
          ml = e0 * *phi - eOne * sinTwoph + eTwo * sin (4.0 *  *phi);
          mlp = e0 - 2.0 * eOne * cos (2.0 *  *phi) + 4.0 * eTwo *  cos (4.0 *  *phi);
         */
        $ml       = $eZero * $phi - $eOne * $sinTwoph + $eTwo * sin(4.0 * $phi)
            - $eThree * sin(6.0 * $phi);
        $mlp      = $eZero - 2.0 * $eOne * cos(2.0 * $phi) + 4.0 * $eTwo *
            cos(4.0 * $phi) - 6.0 * $eThree * cos(6.0 * $phi);
        $conOne   = 2.0 * $ml + $c * ($ml * $ml + $b) - 2.0 * $a *
            ($c * $ml + 1.0);
        $conTwo   = $eccent * $sinTwoph * ($ml * $ml + $b - 2.0 * $a * $ml) /
            (2.0 * $c);
        $conThree = 2.0 * ($a - $ml) *
            ($c * $mlp - 2.0 / $sinTwoph) - 2.0 * $mlp;
        $dphi     = $conOne / ($conTwo + $conThree);
        $phi += $dphi;
        if (abs($dphi) <= .0000000001) return($phi);
    }

    ProjFourphp::reportError("phi4z: No convergence");

    return null;
}

/* Function to compute the constant e4 from the input of the eccentricity
  of the spheroid, x.  This constant is used in the Polar Stereographic
  projection.
  -------------------------------------------------------------------- */

function e4fn($x)
{
    #$con;
    #$com;
    $con = 1.0 + $x;
    $com = 1.0 - $x;
    return (sqrt((pow($con, $con)) * (pow($com, $com))));
}

/* * **************************************************************************
  NAME                             POLYCONIC

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Polyconic projection.  The
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
 * ************************************************************************* */

class ProjFourphp_ProjPoly
{
    /* Initialize the POLYCONIC projection
      ---------------------------------- */

    public function init()
    {
        #$temp;   /* temporary variable		 */
        if ($this->latZero == 0) $this->latZero = 90; //$this->latZero ca

        /* Place parameters in static storage for common use
          ------------------------------------------------- */
        $this->temp = $this->b / $this->a;
        $this->es = 1.0 - pow($this->temp, 2); // devait etre dans tmerc.js
        // mais n y est pas donc je commente sinon retour de valeurs nulles 
        $this->e = sqrt($this->es);
        $this->eZero = ProjFourphp_Common::eZeroFn($this->es);
        $this->eOne = ProjFourphp_Common::eOnefn($this->es);
        $this->eTwo = ProjFourphp_Common::eTwofn($this->es);
        $this->eThree = ProjFourphp_Common::eThreefn($this->es);
        $this->mlZero = ProjFourphp_Common::mlfn($this->eZero, $this->eOne, $this->eTwo, $this->eThree, $this->latZero); //si que des zeros le calcul ne se fait pas
        //if (!$this->mlZero) {$this->mlZero=0;}
    }

    /* Polyconic forward equations--mapping lat,long to x,y
      --------------------------------------------------- */

    public function forward($p)
    {

        /*
          $sinphi;
          $cosphi; // sin and cos value
          $al;    // temporary values
          $c;    // temporary values
          $con;
          $ml;  // cone constant, small m
          $ms;    // small m
          $x;
          $y;
         */

        $lon = $p->x;
        $lat = $p->y;

        $con = ProjFourphp_Common::adjustLon($lon - $this->longZero);

        if (abs($lat) <= .0000001) {
            $x = $this->xZero + $this->a * $con;
            $y = $this->yZero - $this->a * $this->mlZero;
        } else {
            $sinphi = sin($lat);
            $cosphi = cos($lat);

            $ml = ProjFourphp_Common::mlfn($this->eZero, $this->eOne, $this->eTwo, $this->eThree, $lat);
            $ms = ProjFourphp_Common::msfnz($this->e, $sinphi, $cosphi);

            $x = $this->xZero + $this->a * $ms * sin($sinphi) / $sinphi;
            $y = $this->yZero + $this->a * ($ml - $this->mlZero + $ms * (1.0 - cos($sinphi)) / $sinphi);
        }

        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /* Inverse equations
      ----------------- */

    public function inverse($p)
    {

        /*
          $sin_phi;
          $cos_phi; // sin and cos values
          $al;     // temporary values
          $b;     // temporary values
          $c;     // temporary values
          $con;
          $ml;   // cone constant, small m
          $iflg;    // error flag
          $lon;
          $lat;
         */

        $p->x -= $this->xZero;
        $p->y -= $this->yZero;
        $al   = $this->mlZero + $p->y / $this->a;
        $iflg = 0;

        if (abs($al) <= .0000001) {
            $lon = $p->x / $this->a + $this->longZero;
            $lat = 0.0;
        } else {
            $b    = $al * $al + ($p->x / $this->a) * ($p->x / $this->a);
            $iflg = phi4z($this->es, $this->eZero, $this->eOne, $this->eTwo, $this->eThree, $this->al, $b, $c, $lat);
            if ($iflg != 1) return($iflg);
            $lon  = ProjFourphp_Common::adjustLon((ProjFourphp_Common::asinz($p->x * $c / $this->a) / sin($lat)) + $this->longZero);
        }

        $p->x = $lon;
        $p->y = $lat;
        return $p;
    }

}

ProjFourphp::$proj['poly'] = new ProjFourphp_ProjPoly();