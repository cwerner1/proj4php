<?php

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
        $this->C = sqrt(1.0 + $this->es * $cphi * $cphi / (1.0 - $this->es));
        $this->phic0 = asin($sphi / $this->C);
        $this->ratexp = 0.5 * $this->C * $this->e;
        $this->K = tan(0.5 * $this->phic0 + ProjFourphp::$common->fortPi) / (pow(tan(0.5 * $this->latZero + ProjFourphp::$common->fortPi), $this->C) * ProjFourphp::$common->srat($this->e * $sphi, $this->ratexp));
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

        $p->y = 2.0 * atan($this->K * pow(tan(0.5 * $lat + ProjFourphp::$common->fortPi), $this->C) * ProjFourphp::$common->srat($this->e * sin($lat), $this->ratexp)) - ProjFourphp_Common::$halfPi;
        $p->x = $this->C * $lon;

        return $p;
    }

    /**
     *
     * @param type $p
     * @return null 
     */
    public function inverse($p)
    {

        $DEL_TOL = 1e-14;
        $lon = $p->x / $this->C;
        $lat = $p->y;
        $num = pow(tan(0.5 * $lat + ProjFourphp::$common . FORTPI) / $this->K, 1. / $this->C);

        for ($i = ProjFourphp::$common . MAX_ITER; $i > 0; --$i) {
            $lat = 2.0 * atan($num * ProjFourphp::$common->srat($this->e * sin($p->y), -0.5 * $this->e)) - ProjFourphp_Common::$halfPi;
            if (abs($lat - $p->y) < $DEL_TOL)
                break;
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