<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by Proj4php from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
  NAME                            LAMBERT CONFORMAL CONIC

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Lambert Conformal Conic projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.


  ALGORITHM REFERENCES

  1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
 * ***************************************************************************** */


//<2104> +proj=lcc +lat_1=10.16666666666667 +lat_0=10.16666666666667 +lon_0=-71.60561777777777 +k_0=1 +xZero=-17044 +xZero=-23139.97 +ellps=intl +units=m +no_defs  no_defs
// Initialize the Lambert Conformal conic projection
// -----------------------------------------------------------------
//class Proj4phpProjlcc = Class.create();
class Proj4php_ProjLcc
{

    public function init()
    {
        // array of:  r_maj,r_min,latOne,latTwo,c_lon,c_lat,false_east,false_north
        //double c_lat;                   /* center latitude                      */
        //double c_lon;                   /* center longitude                     */
        //double latOne;                    /* first standard parallel              */
        //double latTwo;                    /* second standard parallel             */
        //double r_maj;                   /* major axis                           */
        //double r_min;                   /* minor axis                           */
        //double false_east;              /* x offset in meters                   */
        //double false_north;             /* y offset in meters                   */
        //if latTwo is not defined
        if (!isset($this->latTwo)) {
            $this->latTwo = $this->latZero;
        }

        //if kZero is not defined
        if (!isset($this->kZero))
            $this->kZero = 1.0;

        // Standard Parallels cannot be equal and on opposite sides of the equator
        if (abs($this->latOne + $this->latTwo) < Proj4php_Common::$epsln) {
            Proj4php::reportError("lcc:init: Equal Latitudes");
            return;
        }

        $temp = $this->b / $this->a;
        $this->e = sqrt(1.0 - $temp * $temp);

        $sin1 = sin($this->latOne);
        $cos1 = cos($this->latOne);
        $msOne = Proj4php::$common->msfnz($this->e, $sin1, $cos1);
        $ts1 = Proj4php::$common->tsfnz($this->e, $this->latOne, $sin1);

        $sin2 = sin($this->latTwo);
        $cos2 = cos($this->latTwo);
        $msTwo = Proj4php::$common->msfnz($this->e, $sin2, $cos2);
        $ts2 = Proj4php::$common->tsfnz($this->e, $this->latTwo, $sin2);

        $tsZero = Proj4php::$common->tsfnz($this->e, $this->latZero, sin($this->latZero));

        if (abs($this->latOne - $this->latTwo) > Proj4php_Common::$epsln) {
            $this->ns = log($msOne / $msTwo) / log($ts1 / $ts2);
        } else {
            $this->ns = $sin1;
        }
        $this->f0 = $msOne / ($this->ns * pow($ts1, $this->ns));
        $this->rh = $this->a * $this->f0 * pow($tsZero, $this->ns);

        if (!isset($this->title))
            $this->title = "Lambert Conformal Conic";
    }

    // Lambert Conformal conic forward equations--mapping lat,long to x,y
    // -----------------------------------------------------------------
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        // convert to radians
        if ($lat <= 90.0 && $lat >= -90.0 && $lon <= 180.0 && $lon >= -180.0) {
            //lon = lon * Proj4php::$common.D2R;
            //lat = lat * Proj4php::$common.D2R;
        } else {
            Proj4php::reportError("lcc:forward: llInputOutOfRange: " . $lon . " : " . $lat);
            return null;
        }

        $con = abs(abs($lat) - Proj4php_Common::$halfPi);

        if ($con > Proj4php_Common::$epsln) {
            $ts = Proj4php::$common->tsfnz($this->e, $lat, sin($lat));
            $rhOne = $this->a * $this->f0 * pow($ts, $this->ns);
        } else {
            $con = $lat * $this->ns;
            if ($con <= 0) {
                Proj4php::reportError("lcc:forward: No Projection");
                return null;
            }
            $rhOne = 0;
        }

        $theta = $this->ns * Proj4php_Common::adjustLon($lon - $this->longZero);
        $p->x = $this->kZero * ($rhOne * sin($theta)) + $this->xZero;
        $p->y = $this->kZero * ($this->rh - $rhOne * cos($theta)) + $this->yZero;

        return $p;
    }

    /**
     * Lambert Conformal Conic inverse equations--mapping x,y to lat/long
     * 
     * @param type $p
     * @return null 
     */
    public function inverse($p)
    {

        $x = ($p->x - $this->xZero) / $this->kZero;
        $y = ($this->rh - ($p->y - $this->yZero) / $this->kZero);
        if ($this->ns > 0) {
            $rhOne = sqrt($x * $x + $y * $y);
            $con = 1.0;
        } else {
            $rhOne = -sqrt($x * $x + $y * $y);
            $con = -1.0;
        }
        $theta = 0.0;
        if ($rhOne != 0) {
            $theta = atan2(($con * $x), ($con * $y));
        }
        if (($rhOne != 0) || ($this->ns > 0.0)) {
            $con = 1.0 / $this->ns;
            $ts = pow(($rhOne / ($this->a * $this->f0)), $con);
            $lat = Proj4php::$common->phi2z($this->e, $ts);
            if ($lat == -9999)
                return null;
        } else {
            $lat = -Proj4php_Common::$halfPi;
        }
        $lon = Proj4php_Common::adjustLon($theta / $this->ns + $this->longZero);

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

Proj4php::$proj['lcc'] = new Proj4php_ProjLcc();


