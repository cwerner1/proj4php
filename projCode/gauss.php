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
class ProjFourphp_ProjGauss
{

    /**
     * 
     */
    public function init()
    {

        $sphi = sin($this->latZero);
        $cphi = cos($this->latZero);
        $cphi *= $cphi;
        $this->rc = sqrt(1.0 - $this->es) / (1.0 - $this->es * $sphi * $sphi);
        $this->c = sqrt(1.0 + $this->es * $cphi * $cphi / (1.0 - $this->es));
        $this->phicZero = asin($sphi / $this->c);
        $this->ratexp = 0.5 * $this->c * $this->e;
        $this->k = tan(0.5 * $this->phicZero + ProjFourphp_Common::$fortPi) / (pow(tan(0.5 * $this->latZero + ProjFourphp::$common->fortPi), $this->c) * ProjFourphp::$common->srat($this->e * $sphi, $this->ratexp));
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function forward($p)
    {
        $lon = $p->x;
        $lat = $p->y;

        $p->y = 2.0 *
            atan($this->k * pow(tan(0.5 * $lat + ProjFourphp_Common::$fortPi), $this->c * ProjFourphp_Common::srat($this->e * sin($lat), $this->ratexp)) - ProjFourphp_Common::$halfPi);
        $p->x = $this->c * $lon;

        return $p;
    }

    /**
     *
     * @param type $p
     * @return null 
     */
    public function inverse($p)
    {

        $delTol = 1e-14;
        $lon    = $p->x / $this->c;
        $lat    = $p->y;
        $num    = pow(tan(0.5 * $lat + ProjFourphp_Common::$fortPi) / $this->k, 1. / $this->c);

        for ($i = ProjFourphp_Common::$maxIter; $i > 0; --$i) {
            $lat = 2.0 *
                atan($num * ProjFourphp_Common::srat($this->e * sin($p->y), -0.5 * $this->e)) - ProjFourphp_Common::$halfPi;
            if (abs($lat - $p->y) < $delTol) break;
            $p->y = $lat;
        }

        /* convergence failed */
        if (!$i) {
            ProjFourphp::reportError("gauss:inverse:convergence failed");
            return null;
        }

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['gauss'] = new ProjFourphp_ProjGauss();