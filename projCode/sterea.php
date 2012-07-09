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
class ProjFourphp_ProjSterea
{

    protected $_dependsOn = 'gauss';

    /**
     *
     * @return void 
     */
    public function init()
    {

        if (!$this->rc) {
            ProjFourphp::reportError("sterea:init:E_ERROR_0");
            return;
        }

        $this->sincZero = sin($this->phicZero);
        $this->coscZero = cos($this->phicZero);
        $this->rTwo = 2.0 * $this->rc;

        if (!$this->title)
            $this->title = "Oblique Stereographic Alternative";
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function forward($p)
    {

        $p->x = ProjFourphp_Common::adjustLon($p->x - $this->longZero); /* adjust del longitude */
        $p = ProjFourphp::$proj['gauss']->forward($p);
        $sinc = sin($p->y);
        $cosc = cos($p->y);
        $cosl = cos($p->x);
        $k = $this->kZero * $this->rTwo / (1.0 + $this->sincZero * $sinc + $this->coscZero * $cosc * $cosl);

        $p->x = $k * $cosc * sin($p->x);
        $p->y = $k * ($this->coscZero * sinc - $this->sincZero * $cosc * $cosl);

        $p->x = $this->a * $p->x + $this->xZero;
        $p->y = $this->a * $p->y + $this->yZero;

        return $p;
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function inverse($p)
    {

        #$lon;
        #$lat;
        $p->x = ($p->x - $this->xZero) / $this->a; /* descale and de-offset */
        $p->y = ($p->y - $this->yZero) / $this->a;

        $p->x /= $this->kZero;
        $p->y /= $this->kZero;

        if (($rho = sqrt($p->x * $p->x + $p->y * $p->y))) {
            $c = 2.0 * atan2($rho, $this->rTwo);
            $sinc = sin($c);
            $cosc = cos($c);
            $lat = asin($cosc * $this->sincZero + $p->y * $sinc * $this->coscZero / $rho);
            $lon = atan2($p->x * $sinc, $rho * $this->coscZero * $cosc - $p->y * $this->sincZero * $sinc);
        } else {
            $lat = $this->phicZero;
            $lon = 0.;
        }

        $p->x = $lon;
        $p->y = $lat;
        $p = ProjFourphp::$proj['gauss']->inverse($p);
        $p->x = ProjFourphp_Common::adjustLon($p->x + $this->longZero); /* adjust longitude to CM */

        return $p;
    }

}

ProjFourphp::$proj['sterea'] = new ProjFourphp_ProjSterea();