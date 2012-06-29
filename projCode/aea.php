<?php

/* * *****************************************************************************
  NAME                     ALBERS CONICAL EQUAL AREA

  PURPOSE:	Transforms input longitude and latitude to Easting and Northing
  for the Albers Conical Equal Area projection.  The longitude
  and latitude must be in radians.  The Easting and Northing
  values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  T. Mittan,       	Feb, 1992

  ALGORITHM REFERENCES

  1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
  U.S. Geological Survey Professional Paper 1453 , United State Government
  Printing Office, Washington D.C., 1989.
 * ***************************************************************************** */

/**
 * Author : Julien Moquet
 * 
 * Inspired by Proj4php from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
class Proj4php_ProjAea
{

    /**
     *
     * @return void 
     */
    public function init()
    {

        if (abs($this->latOne + $this->latTwo) < Proj4php_Common::$epsln) {
            Proj4php::reportError("aeaInitEqualLatitudes");
            return;
        }
        $this->temp = $this->b / $this->a;
        $this->es = 1.0 - pow($this->temp, 2);
        $this->e3 = sqrt($this->es);

        $this->sinPo = sin($this->latOne);
        $this->cosPo = cos($this->latOne);
        $this->tOne = $this->sinPo;
        $this->con = $this->sinPo;
        $this->msOne = Proj4php_Common::msfnz($this->e3, $this->sinPo, $this->cosPo);
        $this->qsOne = Proj4php_Common::qsfnz($this->e3, $this->sinPo, $this->cosPo);

        $this->sinPo = sin($this->latTwo);
        $this->cosPo = cos($this->latTwo);
        $this->tTwo = $this->sinPo;
        $this->msTwo = Proj4php_Common::msfnz($this->e3, $this->sinPo, $this->cosPo);
        $this->qsTwo = Proj4php_Common::qsfnz($this->e3, $this->sinPo, $this->cosPo);

        $this->sinPo = sin($this->latZero);
        $this->cosPo = cos($this->latZero);
        $this->t3 = $this->sinPo;
        $this->qsZero = Proj4php_Common::qsfnz($this->e3, $this->sinPo, $this->cosPo);

        if (abs($this->latOne - $this->latTwo) > Proj4php_Common::$epsln) {
            $this->nsZero = ($this->msOne * $this->msOne - $this->msTwo * $this->msTwo) / ($this->qsTwo - $this->qsOne);
        } else {
            $this->nsZero = $this->con;
        }

        $this->c = $this->msOne * $this->msOne + $this->nsZero * $this->qsOne;
        $this->rh = $this->a * sqrt($this->c - $this->nsZero * $this->qsZero) / $this->nsZero;
    }

    /**
     * Albers Conical Equal Area forward equations--mapping lat,long to x,y
     *
     * @param Point $p
     * @return Point $p 
     */
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $this->sinPhi = sin($lat);
        $this->cosPhi = cos($lat);

        $qs = Proj4php_Common::qsfnz($this->e3, $this->sinPhi, $this->cosPhi);
        $rhOne = $this->a * sqrt($this->c - $this->nsZero * $qs) / $this->nsZero;
        $theta = $this->nsZero * Proj4php_Common::adjustLon($lon - $this->longZero);
        $x = rhOne * sin($theta) + $this->xZero;
        $y = $this->rh - $rhOne * cos($theta) + $this->yZero;

        $p->x = $x;
        $p->y = $y;

        return $p;
    }

    /**
     *
     * @param Point $p
     * @return Point $p
     */
    public function inverse($p)
    {

        $p->x -= $this->xZero;
        $p->y = $this->rh - $p->y + $this->yZero;

        if ($this->nsZero >= 0) {
            $rhOne = sqrt($p->x * $p->x + $p->y * $p->y);
            $con = 1.0;
        } else {
            $rhOne = -sqrt($p->x * $p->x + $p->y * $p->y);
            $con = -1.0;
        }

        $theta = 0.0;
        if ($rhOne != 0.0) {
            $theta = atan2($con * $p->x, $con * $p->y);
        }

        $con = $rhOne * $this->nsZero / $this->a;
        $qs = ($this->c - $con * $con) / $this->nsZero;

        if ($this->e3 >= 1e-10) {
            $con = 1 - .5 * (1.0 - $this->es) * log((1.0 - $this->e3) / (1.0 + $this->e3)) / $this->e3;
            if (abs(abs($con) - abs($qs)) > .0000000001) {
                $lat = $this->phi1z($this->e3, $qs);
            } else {
                if ($qs >= 0) {
                    $lat = .5 * Proj4php_Common::$pi;
                } else {
                    $lat = -.5 * Proj4php_Common::$pi;
                }
            }
        } else {
            $lat = $this->phi1z($this->e3, $qs);
        }

        $lon = Proj4php_Common::adjustLon($theta / $this->nsZero + $this->longZero);

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

    /**
     * Function to compute phi1, the latitude for the inverse of the Albers Conical Equal-Area projection.
     *
     * @param type $eccent
     * @param type $qs
     * @return $phi or null on Convergence error
     */
    public function phi1z($eccent, $qs)
    {

        $phi = Proj4php_Common::asinz(.5 * $qs);

        if ($eccent < Proj4php_Common::$epsln)
            return $phi;

        $eccnts = $eccent * $eccent;
        for ($i = 1; $i <= 25; ++$i) {
            $sinphi = sin($phi);
            $cosphi = cos($phi);
            $con = $eccent * $sinphi;
            $com = 1.0 - $con * $con;
            $dphi = .5 * $com * $com / $cosphi * ($qs / (1.0 - $eccnts) - $sinphi / $com + .5 / $eccent * log((1.0 - $con) / (1.0 + $con)));
            $phi = $phi + $dphi;
            if (abs($dphi) <= 1e-7)
                return $phi;
        }

        Proj4php::reportError("aea:phi1z:Convergence error");

        return null;
    }

}

Proj4php::$proj['aea'] = new Proj4php_ProjAea();