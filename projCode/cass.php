<?php

/* * *****************************************************************************
  NAME                            CASSINI

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Cassini projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.
  Ported from PROJ.4.


  ALGORITHM REFERENCES

  1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
 * ***************************************************************************** */

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
//ProjFourphp.defs["EPSG:28191"] = "+proj=cass +lat_0=31.73409694444445 +lon_0=35.21208055555556 +x_0=170251.555 +y_0=126867.909 +a=6378300.789 +b=6356566.435 +towgs84=-275.722,94.7824,340.894,-8.001,-4.42,-11.821,1 +units=m +no_defs";
// Initialize the Cassini projection
// -----------------------------------------------------------------

class ProjFourphp_ProjCass
{

    public function init()
    {
        if (!$this->sphere) {
            $this->en = ProjFourphp::$common->pjEnfn($this->es);
            $this->mZero = ProjFourphp_Common::jMlfn($this->latZero, sin($this->latZero), cos($this->latZero), $this->en);
        }
    }

    protected $C1 = .16666666666666666666;
    protected $C2 = .00833333333333333333;
    protected $C3 = .04166666666666666666;
    protected $C4 = .33333333333333333333;
    protected $C5 = .06666666666666666666;

    /* Cassini forward equations--mapping lat,long to x,y
      ----------------------------------------------------------------------- */

    public function forward($p)
    {

        /* Forward equations
          ----------------- */
        #$x;
        #$y;
        $lam = $p->x;
        $phi = $p->y;
        $lam =  ProjFourphp_Common::adjustLon($lam - $this->longZero);

        if ($this->sphere) {
            $x = asin(cos($phi) * sin($lam));
            $y = atan2(tan($phi), cos($lam)) - $this->phiZero;
        } else {
            //ellipsoid
            $this->n = sin($phi);
            $this->c = cos($phi);
            $y = $this->pj_mlfn($phi, $this->n, $this->c, $this->en);
            $this->n = 1. / sqrt(1. - $this->es * $this->n * $this->n);
            $this->tn = tan($phi);
            $this->t = $this->tn * $this->tn;
            $this->aOne = $lam * $this->c;
            $this->c *= $this->es * $this->c / (1 - $this->es);
            $this->aTwo = $this->aOne * $this->aOne;
            $x = $this->n * $this->aOne * (1. - $this->aTwo * $this->t * ($this->C1 - (8. - $this->t + 8. * $this->c) * $this->aTwo * $this->C2));
            $y -= $this->mZero - $this->n * $this->tn * $this->aTwo * (.5 + (5. - $this->t + 6. * $this->c) * $this->aTwo * $this->C3);
        }

        $p->x = $this->a * $x + $this->xZero;
        $p->y = $this->a * $y + $this->yZero;

        return $p;
    }

    /* Inverse equations
      ----------------- */

    public function inverse($p)
    {
        $p->x -= $this->xZero;
        $p->y -= $this->yZero;
        $x = $p->x / $this->a;
        $y = $p->y / $this->a;

        if ($this->sphere) {
            $this->dd = $y + $this->latZero;
            $phi = asin(sin($this->dd) * cos($x));
            $lam = atan2(tan($x), cos($this->dd));
        } else {
            /* ellipsoid */
            $phOne = ProjFourphp::$common->pjInvMlfn($this->mZero + $y, $this->es, $this->en);
            $this->tn = tan($phOne);
            $this->t = $this->tn * $this->tn;
            $this->n = sin($phOne);
            $this->r = 1. / (1. - $this->es * $this->n * $this->n);
            $this->n = sqrt($this->r);
            $this->r *= (1. - $this->es) * $this->n;
            $this->dd = $x / $this->n;
            $this->dTwo = $this->dd * $this->dd;
            $phi = $phOne - ($this->n * $this->tn / $this->r) * $this->dTwo * (.5 - (1. + 3. * $this->t) * $this->dTwo * $this->C3);
            $lam = $this->dd * (1. + $this->t * $this->dTwo * (-$this->C4 + (1. + 3. * $this->t) * $this->dTwo * $this->C5)) / cos($phOne);
        }
        $p->x =  ProjFourphp_Common::adjustLon($this->longZero + $lam);
        $p->y = $phi;

        return $p;
    }

}

ProjFourphp::$proj['cass'] = new ProjFourphp_ProjCass();