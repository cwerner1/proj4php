<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by Proj4php from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
// Initialize the Stereographic projection
class Proj4php_ProjStere
{

    static protected $tol = 1.e-8;
    static protected $niter = 8;
    static protected $conv = 1.e-10;
    static protected $sPole = 0;
    static protected $nPole = 1;
    static protected $obliq = 2;
    static protected $equit = 3;

    /**
     *
     * @param type $phit
     * @param type $sinphi
     * @param type $eccen
     * @return type 
     */
    public function ssfn_($phit, $sinphi, $eccen)
    {
        $sinphi *= $eccen;
        return (tan(.5 * (Proj4php_Common::$halfPi + $phit)) * pow((1. - $sinphi) / (1. + $sinphi), .5 * $eccen));
    }

    /**
     * 
     */
    public function init()
    {
        $this->phits = $this->latTs ? $this->latTs : Proj4php_Common::$halfPi;
        $t = abs($this->latZero);
        if ((abs($t) - Proj4php_Common::$halfPi) < Proj4php_Common::$epsln) {
            $this->mode = $this->latZero < 0. ? self::$sPole : self::$nPole;
        } else {
            $this->mode = $t > Proj4php_Common::$epsln ? self::$obliq : self::$equit;
        }
        $this->phits = abs($this->phits);
        if ($this->es) {
            #$X;

            switch ($this->mode) {
                case self::$nPole:
                case self::$sPole:
                    if (abs($this->phits - Proj4php_Common::$halfPi) < Proj4php_Common::$epsln) {
                        $this->akmOne = 2. * $this->kZero / sqrt(pow(1 + $this->e, 1 + $this->e) * pow(1 - $this->e, 1 - $this->e));
                    } else {
                        $t = sin($this->phits);
                        $this->akmOne = cos($this->phits) / Proj4php::$common->tsfnz($this->e, $this->phits, $t);
                        $t *= $this->e;
                        $this->akmOne /= sqrt(1. - $t * $t);
                    }
                    break;
                case self::$equit:
                    $this->akmOne = 2. * $this->kZero;
                    break;
                case self::$obliq:
                    $t = sin($this->latZero);
                    $X = 2. * atan($this->ssfn_($this->latZero, $t, $this->e)) - Proj4php_Common::$halfPi;
                    $t *= $this->e;
                    $this->akmOne = 2. * $this->kZero * cos($this->latZero) / sqrt(1. - $t * $t);
                    $this->sinX1 = sin($X);
                    $this->cosX1 = cos($X);
                    break;
            }
        } else {
            switch ($this->mode) {
                case self::$obliq:
                    $this->sinphZero = sin($this->latZero);
                    $this->cosphZero = cos($this->latZero);
                case self::$equit:
                    $this->akmOne = 2. * $this->kZero;
                    break;
                case self::$sPole:
                case self::$nPole:
                    $this->akmOne = abs($this->phits - Proj4php_Common::$halfPi) >= Proj4php_Common::$epsln ?
                            cos($this->phits) / tan(Proj4php::$common->fortPi - .5 * $this->phits) :
                            2. * $this->kZero;
                    break;
            }
        }
    }

    /**
     * Stereographic forward equations--mapping lat,long to x,y
     *
     * @param type $p
     * @return type 
     */
    public function forward($p)
    {

        $lon = $p->x;
        $lon = Proj4php_Common::adjustLon($lon - $this->longZero);
        $lat = $p->y;
        #$x;
        #$y;

        if ($this->sphere) {
            /*
              $sinphi;
              $cosphi;
              $coslam;
              $sinlam;
             */

            $sinphi = sin($lat);
            $cosphi = cos($lat);
            $coslam = cos($lon);
            $sinlam = sin($lon);
            switch ($this->mode) {
                case self::$equit:
                    $y = 1. + $cosphi * $coslam;
                    if (y <= Proj4php_Common::$epsln) {
                        Proj4php::reportError("stere:forward:Equit");
                    }
                    $y = $this->akmOne / $y;
                    $x = $y * $cosphi * $sinlam;
                    $y *= $sinphi;
                    break;
                case self::$obliq:
                    $y = 1. + $this->sinphZero * $sinphi + $this->cosphZero * $cosphi * $coslam;
                    if ($y <= Proj4php_Common::$epsln) {
                        Proj4php::reportError("stere:forward:Obliq");
                    }
                    $y = $this->akmOne / $y;
                    $x = $y * $cosphi * $sinlam;
                    $y *= $this->cosphZero * $sinphi - $this->sinphZero * $cosphi * $coslam;
                    break;
                case self::$nPole:
                    $coslam = -$coslam;
                    $lat = -$lat;
                //Note  no break here so it conitnues through S_POLE
                case self::$sPole:
                    if (abs($lat - Proj4php_Common::$halfPi) < self::$tol) {
                        Proj4php::reportError("stere:forward:S_POLE");
                    }
                    $y = $this->akmOne * tan(Proj4php::$common->fortPi + .5 * $lat);
                    $x = $sinlam * $y;
                    $y *= $coslam;
                    break;
            }
        } else {
            $coslam = cos($lon);
            $sinlam = sin($lon);
            $sinphi = sin($lat);
            if ($this->mode == self::$obliq || $this->mode == self::$equit) {
                $Xt = 2. * atan($this->ssfn_($lat, $sinphi, $this->e));
                $sinX = sin($Xt - Proj4php_Common::$halfPi);
                $cosX = cos($Xt);
            }
            switch ($this->mode) {
                case self::$obliq:
                    $A = $this->akmOne / ($this->cosX1 * (1. + $this->sinX1 * $sinX + $this->cosX1 * $cosX * $coslam));
                    $y = $A * ($this->cosX1 * $sinX - $this->sinX1 * $cosX * $coslam);
                    $x = $A * $cosX;
                    break;
                case self::$equit:
                    $A = 2. * $this->akmOne / (1. + $cosX * $coslam);
                    $y = $A * $sinX;
                    $x = $A * $cosX;
                    break;
                case self::$sPole:
                    $lat = -$lat;
                    $coslam = - $coslam;
                    $sinphi = -$sinphi;
                case self::$nPole:
                    $x = $this->akmOne * Proj4php::$common->tsfnz($this->e, $lat, $sinphi);
                    $y = - $x * $coslam;
                    break;
            }
            $x = $x * $sinlam;
        }

        $p->x = $x * $this->a + $this->xZero;
        $p->y = $y * $this->a + $this->yZero;

        return $p;
    }

    /**
     * Stereographic inverse equations--mapping x,y to lat/long
     *
     * @param type $p
     * @return type 
     */
    public function inverse($p)
    {
        $x = ($p->x - $this->xZero) / $this->a;   /* descale and de-offset */
        $y = ($p->y - $this->yZero) / $this->a;
        /*
          $lon;
          $lat;
          $cosphi;
          $sinphi;
          $rho;
          $tp = 0.0;
          $phi_l = 0.0;
          $i;
         */
        $halfe = 0.0;
        $pi2 = 0.0;

        if ($this->sphere) {
            /*
              $c;
              $rh;
              $sinc;
              $cosc;
             */

            $rh = sqrt($x * $x + $y * $y);
            $c = 2. * atan($rh / $this->akmOne);
            $sinc = sin($c);
            $cosc = cos($c);
            $lon = 0.;

            switch ($this->mode) {
                case self::$equit:
                    if (abs($rh) <= Proj4php_Common::$epsln) {
                        $lat = 0.;
                    } else {
                        $lat = asin($y * $sinc / $rh);
                    }
                    if ($cosc != 0. || $x != 0.)
                        $lon = atan2($x * $sinc, $cosc * $rh);
                    break;
                case self::$obliq:
                    if (abs($rh) <= Proj4php_Common::$epsln) {
                        $lat = $this->phi0;
                    } else {
                        $lat = asin($cosc * $this->sinphZero + $y * $sinc * $this->cosphZero / $rh);
                    }
                    $c = $cosc - $this->sinphZero * sin($lat);
                    if ($c != 0. || $x != 0.) {
                        $lon = atan2($x * $sinc * $this->cosphZero, $c * $rh);
                    }
                    break;
                case self::$nPole:
                    $y = -$y;
                case self::$sPole:
                    if (abs($rh) <= Proj4php_Common::$epsln) {
                        $lat = $this->phi0;
                    } else {
                        $lat = asin($this->mode == self::$sPole ? -$cosc : $cosc );
                    }
                    $lon = ($x == 0. && $y == 0.) ? 0. : atan2($x, $y);
                    break;
            }
            $p->x = Proj4php_Common::adjustLon($lon + $this->longZero);
            $p->y = $lat;
        } else {
            $rho = sqrt($x * $x + $y * $y);
            switch ($this->mode) {
                case self::$obliq:
                case self::$equit:
                    $tp = 2. * atan2($rho * $this->cosX1, $this->akmOne);
                    $cosphi = cos($tp);
                    $sinphi = sin($tp);
                    if ($rho == 0.0) {
                        $phi_l = asin($cosphi * $this->sinX1);
                    } else {
                        $phi_l = asin($cosphi * $this->sinX1 + ($y * $sinphi * $this->cosX1 / $rho));
                    }

                    $tp = tan(.5 * (Proj4php_Common::$halfPi + $phi_l));
                    $x *= $sinphi;
                    $y = $rho * $this->cosX1 * $cosphi - $y * $this->sinX1 * $sinphi;
                    $pi2 = Proj4php_Common::$halfPi;
                    $halfe = .5 * $this->e;
                    break;
                case self::$nPole:
                    $y = -$y;
                case self::$sPole:
                    $tp = - $rho / $this->akmOne;
                    $phi_l = Proj4php_Common::$halfPi - 2. * atan($tp);
                    $pi2 = -Proj4php_Common::$halfPi;
                    $halfe = -.5 * $this->e;
                    break;
            }
            for ($i = $this->NITER; $i--; $phi_l = $lat) { //check this
                $sinphi = $this->e * sin($phi_l);
                $lat = 2. * atan($tp * pow((1. + $sinphi) / (1. - $sinphi), $halfe)) - $pi2;
                if (abs(phi_l - lat) < $this->CONV) {
                    if ($this->mode == self::$sPole)
                        $lat = -$lat;
                    $lon = ($x == 0. && $y == 0.) ? 0. : atan2($x, $y);
                    $p->x = Proj4php_Common::adjustLon($lon + $this->longZero);
                    $p->y = $lat;
                    return $p;
                }
            }
        }
    }

}

Proj4php::$proj['stere'] = new Proj4php_ProjStere();