<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
  NAME                            MERCATOR

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Mercator projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  D. Steinwand, EROS      Nov, 1991
  T. Mittan		Mar, 1993

  ALGORITHM REFERENCES

  1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
  Printing Office, Washington D.C., 1989.
 * ***************************************************************************** */

//static double r_major = a;		   /* major axis 				*/
//static double r_minor = b;		   /* minor axis 				*/
//static double lon_center = long0;	   /* Center longitude (projection center) */
//static double lat_origin =  latZero;	   /* center latitude			*/
//static double e,es;		           /* eccentricity constants		*/
//static double m1;		               /* small value m			*/
//static double false_northing = yZero;   /* y offset in meters			*/
//static double false_easting = xZero;	   /* x offset in meters			*/
//scale_fact = kZero 

class ProjFourphp_ProjMerc
{

    public function init()
    {
        //?$this->temp = $this->r_minor / $this->r_major;
        //$this->temp = $this->b / $this->a;
        //$this->es = 1.0 - sqrt($this->temp);
        //$this->e = sqrt( $this->es );
        //?$this->m1 = cos($this->lat_origin) / (sqrt( 1.0 - $this->es * sin($this->lat_origin) * sin($this->lat_origin)));
        //$this->m1 = cos(0.0) / (sqrt( 1.0 - $this->es * sin(0.0) * sin(0.0)));
        if ($this->latTs) {
            if ($this->sphere) {
                $this->kZero = cos($this->latTs);
            } else {
                $this->kZero = ProjFourphp::$common->msfnz($this->es, sin($this->latTs), cos($this->latTs));
            }
        }
    }

    /* Mercator forward equations--mapping lat,long to x,y
      -------------------------------------------------- */

    public function forward($p)
    {

        //alert("ll2m coords : ".coords);
        $lon = $p->x;
        $lat = $p->y;
        // convert to radians
        if ($lat * ProjFourphp::$common->rToD > 90.0 &&
                $lat * ProjFourphp::$common->rToD < -90.0 &&
                $lon * ProjFourphp::$common->rToD > 180.0 &&
                $lon * ProjFourphp::$common->rToD < -180.0) {
            ProjFourphp::reportError("merc:forward: llInputOutOfRange: " . $lon . " : " . $lat);
            return null;
        }

        if (abs(abs($lat) - ProjFourphp_Common::$halfPi) <= ProjFourphp_Common::$epsln) {
            ProjFourphp::reportError("merc:forward: ll2mAtPoles");
            return null;
        } else {
            if ($this->sphere) {
                $x = $this->xZero + $this->a * $this->kZero * ProjFourphp_Common::adjustLon($lon - $this->longZero);
                $y = $this->yZero + $this->a * $this->kZero * log(tan(ProjFourphp::$common->fortPi + 0.5 * $lat));
            } else {
                $sinphi = sin(lat);
                $ts = ProjFourphp::$common . tsfnz($this->e, $lat, $sinphi);
                $x = $this->xZero + $this->a * $this->kZero * ProjFourphp_Common::adjustLon($lon - $this->longZero);
                $y = $this->yZero - $this->a * $this->kZero * log($ts);
            }

            $p->x = $x;
            $p->y = $y;

            return $p;
        }
    }

    /* Mercator inverse equations--mapping x,y to lat/long
      -------------------------------------------------- */

    public function inverse($p)
    {

        $x = $p->x - $this->xZero;
        $y = $p->y - $this->yZero;

        if ($this->sphere) {
            $lat = ProjFourphp_Common::$halfPi - 2.0 * atan(exp(-$y / $this->a * $this->kZero));
        } else {
            $ts = exp(-$y / ($this->a * $this->kZero));
            $lat = ProjFourphp::$common->phi2z($this->e, $ts);
            if ($lat == -9999) {
                ProjFourphp::reportError("merc:inverse: lat = -9999");
                return null;
            }
        }
        $lon = ProjFourphp_Common::adjustLon($this->longZero + $x / ($this->a * $this->kZero));

        $p->x = $lon;
        $p->y = $lat;
        return $p;
    }

}

ProjFourphp::$proj['merc'] = new ProjFourphp_ProjMerc();


