<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
  NAME                       OBLIQUE MERCATOR (HOTINE)

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Oblique Mercator projection.  The
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

class ProjFourphp_ProjOmerc
{
    /* Initialize the Oblique Mercator  projection
      ------------------------------------------ */

    public function init()
    {
        if (!$this->mode) {
            $this->mode = 0;
        }
        if (!$this->lon1) {
            $this->lon1 = 0;
            $this->mode = 1;
        }
        if (!$this->lon2) {
            $this->lon2 = 0;
        }
        if (!$this->latTwo) {
            $this->latTwo = 0;
        }

        /* Place parameters in static storage for common use
          ------------------------------------------------- */
        $temp = $this->b / $this->a;
        $es = 1.0 - pow($temp, 2);
        $e = sqrt($es);

        $this->sinP20 = sin($this->latZero);
        $this->cosP20 = cos($this->latZero);

        $this->con = 1.0 - $this->es * $this->sinP20 * $this->sinP20;
        $this->com = sqrt(1.0 - $es);
        $this->bl = sqrt(1.0 + $this->es * pow($this->cosP20, 4.0) / (1.0 - $es));
        $this->al = $this->a * $this->bl * $this->kZero * $this->com / $this->con;
        if (abs($this->latZero) < ProjFourphp_Common::$epsln) {
            $this->ts = 1.0;
            $this->d = 1.0;
            $this->el = 1.0;
        } else {
            $this->ts = ProjFourphp::$common->tsfnz($this->e, $this->latZero, $this->sinP20);
            $this->con = sqrt($this->con);
            $this->d = $this->bl * $this->com / ($this->cosP20 * $this->con);
            if (($this->d * $this->d - 1.0) > 0.0) {
                if ($this->latZero >= 0.0) {
                    $this->f = $this->d + sqrt($this->d * $this->d - 1.0);
                } else {
                    $this->f = $this->d - sqrt($this->d * $this->d - 1.0);
                }
            } else {
                $this->f = $this->d;
            }
            $this->el = $this->f * pow($this->ts, $this->bl);
        }

        //$this->longc=52.60353916666667;

        if ($this->mode != 0) {
            $this->g = .5 * ($this->f - 1.0 / $this->f);
            $this->gama = ProjFourphp_Common::asinz(sin($this->alpha) / $this->d);
            $this->longc = $this->longc - ProjFourphp_Common::asinz($this->g * tan($this->gama)) / $this->bl;

            /* Report parameters common to format B
              ------------------------------------- */
            //genrpt(azimuth * R2D,"Azimuth of Central Line:    ");
            //cenlon(lon_origin);
            // cenlat(lat_origin);

            $this->con = abs($this->latZero);
            if (($this->con > ProjFourphp_Common::$epsln) && (abs($this->con - ProjFourphp_Common::$halfPi) > ProjFourphp_Common::$epsln)) {
                $this->singam = sin($this->gama);
                $this->cosgam = cos($this->gama);

                $this->sinaz = sin($this->alpha);
                $this->cosaz = cos($this->alpha);

                if ($this->latZero >= 0) {
                    $this->u = ($this->al / $this->bl) * atan(sqrt($this->d * $this->d - 1.0) / $this->cosaz);
                } else {
                    $this->u = -($this->al / $this->bl) * atan(sqrt($this->d * $this->d - 1.0) / $this->cosaz);
                }
            } else {
                ProjFourphp::reportError("omerc:Init:DataError");
            }
        } else {
            $this->sinphi = sin($this->atOne);
            $this->ts1 = ProjFourphp::$common->tsfnz($this->e, $this->latOne, $this->sinphi);
            $this->sinphi = sin($this->latTwo);
            $this->ts2 = ProjFourphp::$common->tsfnz($this->e, $this->latTwo, $this->sinphi);
            $this->h = pow($this->ts1, $this->bl);
            $this->l = pow($this->ts2, $this->bl);
            $this->f = $this->el / $this->h;
            $this->g = .5 * ($this->f - 1.0 / $this->f);
            $this->j = ($this->el * $this->el - $this->l * $this->h) / ($this->el * $this->el + $this->l * $this->h);
            $this->p = ($this->l - $this->h) / ($this->l + $this->h);
            $this->dlon = $this->lon1 - $this->lon2;
            if ($this->dlon < -ProjFourphp::$common->pi)
                $this->lon2 = $this->lon2 - 2.0 * ProjFourphp::$common->pi;
            if ($this->dlon > ProjFourphp::$common->pi)
                $this->lon2 = $this->lon2 + 2.0 * ProjFourphp::$common->pi;
            $this->dlon = $this->lon1 - $this->lon2;
            $this->longc = .5 * ($this->lon1 + $this->lon2) - atan($this->j * tan(.5 * $this->bl * $this->dlon) / $this->p) / $this->bl;
            $this->dlon = ProjFourphp_Common::adjustLon($this->lon1 - $this->longc);
            $this->gama = atan(sin($this->bl * $this->dlon) / $this->g);
            $this->alpha = ProjFourphp_Common::asinz($this->d * sin($this->gama));

            /* Report parameters common to format A
              ------------------------------------- */
            if (abs($this->latOne - $this->latTwo) <= ProjFourphp_Common::$epsln) {
                ProjFourphp::reportError("omercInitDataError");
                //return(202);
            } else {
                $this->con = abs($this->latOne);
            }
            if (($this->con <= ProjFourphp_Common::$epsln) || (abs($this->con - ProjFourphp_Common::$halfPi) <= ProjFourphp_Common::$epsln)) {
                ProjFourphp::reportError("omercInitDataError");
                //return(202);
            } else {
                if (abs(abs($this->latZero) - ProjFourphp_Common::$halfPi) <= ProjFourphp_Common::$epsln) {
                    ProjFourphp::reportError("omercInitDataError");
                    //return(202);
                }
            }

            $this->singam = sin($this->gam);
            $this->cosgam = cos($this->gam);

            $this->sinaz = sin($this->alpha);
            $this->cosaz = cos($this->alpha);


            if ($this->latZero >= 0) {
                $this->u = ($this->al / $this->bl) * atan(sqrt($this->d * $this->d - 1.0) / $this->cosaz);
            } else {
                $this->u = -($this->al / $this->bl) * atan(sqrt($this->d * $this->d - 1.0) / $this->cosaz);
            }
        }
    }

    /* Oblique Mercator forward equations--mapping lat,long to x,y
      ---------------------------------------------------------- */

    public function forward($p)
    {

        /*
          $theta;   // angle
          $sin_phi;
          $cos_phi;  // sin and cos value
          $b;  // temporary values
          $c;
          $t;
          $tq; // temporary values
          $con;
          $n;
          $ml; // cone constant, small m
          $q;
          $us;
          $vl;
          $ul;
          $vs;
          $s;
          $dlon;
          $ts1;
         */

        $lon = $p->x;
        $lat = $p->y;

        /* Forward equations
          ----------------- */
        $sin_phi = sin($lat);
        $dlon = ProjFourphp_Common::adjustLon($lon - $this->longc);
        $vl = sin($this->bl * $dlon);
        if (abs(abs($lat) - ProjFourphp_Common::$halfPi) > ProjFourphp_Common::$epsln) {
            $ts1 = ProjFourphp::$common->tsfnz($this->e, $lat, $sin_phi);
            $q = $this->el / (pow($ts1, $this->bl));
            $s = .5 * ($q - 1.0 / $q);
            $t = .5 * ($q + 1.0 / $q);
            $ul = ($s * $this->singam - $vl * $this->cosgam) / $t;
            $con = cos($this->bl * $dlon);
            if (abs(con) < .0000001) {
                $us = $this->al * $this->bl * $dlon;
            } else {
                $us = $this->al * atan(($s * $this->cosgam + $vl * $this->singam) / $con) / $this->bl;
                if ($con < 0)
                    $us = $us + ProjFourphp::$common->pi * $this->al / $this->bl;
            }
        } else {
            if ($lat >= 0) {
                $ul = $this->singam;
            } else {
                $ul = -$this->singam;
            }
            $us = $this->al * $lat / $this->bl;
        }
        if (abs(abs($ul) - 1.0) <= ProjFourphp_Common::$epsln) {
            //alert("Point projects into infinity","omer-for");
            ProjFourphp::reportError("omercFwdInfinity");
            //return(205);
        }
        $vs = .5 * $this->al * log((1.0 - $ul) / (1.0 + $ul)) / $this->bl;
        $us = $us - $this->u;
        $p->x = $this->xZero + $vs * $this->cosaz + $us * $this->sinaz;
        $p->y = $this->yZero + $us * $this->cosaz - $vs * $this->sinaz;

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
          $deltaLon; /* Delta longitude (Given longitude - center
          $theta;  /* angle
          $delta_theta; /* adjusted longitude
          $sin_phi;
          $cos_phi; /* sin and cos value
          $b;  /* temporary values
          $c;
          $t;
          $tq; /* temporary values
          $con;
          $n;
          $ml; /* cone constant, small m
          $vs;
          $us;
          $q;
          $s;
          $ts1;
          $vl;
          $ul;
          $bs;
          $dlon;
          $flag;
         */

        /* Inverse equations
          ----------------- */
        $p->x -= $this->xZero;
        $p->y -= $this->yZero;
        #$flag = 0;
        $vs = $p->x * $this->cosaz - $p->y * $this->sinaz;
        $us = $p->y * $this->cosaz + $p->x * $this->sinaz;
        $us = $us + $this->u;
        $q = exp(-$this->bl * $vs / $this->al);
        $s = .5 * ($q - 1.0 / $q);
        $t = .5 * ($q + 1.0 / $q);
        $vl = sin($this->bl * $us / $this->al);
        $ul = ($vl * $this->cosgam + $s * $this->singam) / $t;
        if (abs(abs($ul) - 1.0) <= ProjFourphp_Common::$epsln) {
            $lon = $this->longc;
            if (ul >= 0.0) {
                $lat = ProjFourphp_Common::$halfPi;
            } else {
                $lat = -ProjFourphp_Common::$halfPi;
            }
        } else {
            $con = 1.0 / $this->bl;
            $ts1 = pow(($this->el / sqrt((1.0 + $ul) / (1.0 - $ul))), $con);
            $lat = ProjFourphp::$common->phi2z($this->e, $ts1);
            //if ($flag != 0)
            //return($flag);
            //~ con = cos($this->bl * us /al);
            $theta = $this->longc - atan2(($s * $this->cosgam - $vl * $this->singam), $con) / $this->bl;
            $lon = ProjFourphp_Common::adjustLon($theta);
        }
        $p->x = $lon;
        $p->y = $lat;
        return $p;
    }

}

ProjFourphp::$proj['omerc'] = new ProjFourphp_ProjOmerc();