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
// Initialize the Stereographic projection
class ProjFourphp_ProjStere
{

    static protected $_tol   = 1.e-8;
    static protected $_niter = 8;
    static protected $_conv  = 1.e-10;
    static protected $_sPole = 0;
    static protected $_nPole = 1;
    static protected $_obliq = 2;
    static protected $_equit = 3;

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
        return (tan(.5
                * (ProjFourphp_Common::$halfPi + $phit))
            * pow((1. - $sinphi) / (1. + $sinphi), .5 * $eccen));
    }

    /**
     * 
     */
    public function init()
    {
        if ($this->latTs) {
            $this->phits = $this->latTs;
        } else {
            $this->phits = ProjFourphp_Common::$halfPi;
        }
        $t           = abs($this->latZero);
        if ((abs($t) - ProjFourphp_Common::$halfPi)
            < ProjFourphp_Common::$epsln) {
            $this->mode = $this->latZero < 0. ? self::$_sPole : self::$_nPole;
        } else {
            if ($t > ProjFourphp_Common::$epsln) {
                $this->mode = self::$_obliq;
            } else {
                $this->mode  = self::$_equit;
            }
        }
        $this->phits = abs($this->phits);
        if ($this->es) {
            #$X;

            switch ($this->mode) {
                case self::$_nPole:
                case self::$_sPole:
                    if (abs($this->phits - ProjFourphp_Common::$halfPi)
                        < ProjFourphp_Common::$epsln) {
                        $this->akmOne = 2. *
                            $this->kZero / sqrt(pow(1 + $this->e, 1 + $this->e)
                                * pow(1 - $this->e, 1 - $this->e));
                    } else {
                        $t             = sin($this->phits);
                        $this->akmOne  =
                            cos($this->phits) /
                            ProjFourphp_Common::tsfnz($this->e, $this->phits, $t);
                        $t *= $this->e;
                        $this->akmOne /= sqrt(1. - $t * $t);
                    }
                    break;
                case self::$_equit:
                    $this->akmOne  = 2. * $this->kZero;
                    break;
                case self::$_obliq:
                    $t             = sin($this->latZero);
                    $x             =
                        2. * atan($this->ssfn_($this->latZero, $t, $this->e))
                        - ProjFourphp_Common::$halfPi;
                    $t *= $this->e;
                    $this->akmOne  = 2. * $this->kZero * cos($this->latZero)
                        / sqrt(1. - $t * $t);
                    $this->sinXOne = sin($x);
                    $this->cosXOne = cos($x);
                    break;
            }
        } else {
            switch ($this->mode) {
                case self::$_obliq:
                    $this->sinphZero = sin($this->latZero);
                    $this->cosphZero = cos($this->latZero);
                case self::$_equit:
                    $this->akmOne    = 2. * $this->kZero;
                    break;
                case self::$_sPole:
                case self::$_nPole:
                    $this->akmOne    =
                        abs($this->phits - ProjFourphp_Common::$halfPi)
                        >= ProjFourphp_Common::$epsln ?
                        cos($this->phits)
                        / tan(ProjFourphp::$common->fortPi - .5 * $this->phits) :
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
        $lon = ProjFourphp_Common::adjustLon($lon - $this->longZero);
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
                case self::$_equit:
                    $y = 1. + $cosphi * $coslam;
                    if (y <= ProjFourphp_Common::$epsln) {
                        ProjFourphp::reportError("stere:forward:Equit");
                    }
                    $y = $this->akmOne / $y;
                    $x = $y * $cosphi * $sinlam;
                    $y *= $sinphi;
                    break;
                case self::$_obliq:
                    $y = 1. + $this->sinphZero * $sinphi + $this->cosphZero *
                        $cosphi * $coslam;
                    if ($y <= ProjFourphp_Common::$epsln) {
                        ProjFourphp::reportError("stere:forward:Obliq");
                    }
                    $y      = $this->akmOne / $y;
                    $x      = $y * $cosphi * $sinlam;
                    $y *= $this->cosphZero * $sinphi - $this->sinphZero *
                        $cosphi * $coslam;
                    break;
                case self::$_nPole:
                    $coslam = -$coslam;
                    $lat    = -$lat;
                //Note  no break here so it conitnues through S_POLE
                case self::$_sPole:
                    if (abs($lat - ProjFourphp_Common::$halfPi) < self::$_tol) {
                        ProjFourphp::reportError("stere:forward:S_POLE");
                    }
                    $y = $this->akmOne *
                        tan(ProjFourphp::$common->fortPi + .5 * $lat);
                    $x = $sinlam * $y;
                    $y *= $coslam;
                    break;
            }
        } else {
            $coslam = cos($lon);
            $sinlam = sin($lon);
            $sinphi = sin($lat);
            if ($this->mode == self::$_obliq
                || $this->mode == self::$_equit) {
                $xt   = 2. * atan($this->ssfn_($lat, $sinphi, $this->e));
                $sinX = sin($xt - ProjFourphp_Common::$halfPi);
                $cosX = cos($xt);
            }
            switch ($this->mode) {
                case self::$_obliq:
                    $a      = $this->akmOne / ($this->cosXOne *
                        (1. + $this->sinXOne * $sinX +
                        $this->cosXOne * $cosX * $coslam));
                    $y      = $a * ($this->cosXOne * $sinX -
                        $this->sinXOne * $cosX * $coslam);
                    $x      = $a * $cosX;
                    break;
                case self::$_equit:
                    $a      = 2. * $this->akmOne / (1. + $cosX * $coslam);
                    $y      = $a * $sinX;
                    $x      = $a * $cosX;
                    break;
                case self::$_sPole:
                    $lat    = -$lat;
                    $coslam = - $coslam;
                    $sinphi = -$sinphi;
                case self::$_nPole:
                    $x      = $this->akmOne *
                        ProjFourphp_Common::tsfnz($this->e, $lat, $sinphi);
                    $y      = - $x * $coslam;
                    break;
            }
            $x      = $x * $sinlam;
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
        $x     = ($p->x - $this->xZero) / $this->a;   /* descale and de-offset */
        $y     = ($p->y - $this->yZero) / $this->a;
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
        $piTwo = 0.0;

        if ($this->sphere) {
            /*
              $c;
              $rh;
              $sinc;
              $cosc;
             */

            $rh   = sqrt($x * $x + $y * $y);
            $c    = 2. * atan($rh / $this->akmOne);
            $sinc = sin($c);
            $cosc = cos($c);
            $lon  = 0.;

            switch ($this->mode) {
                case self::$_equit:
                    if (abs($rh) <= ProjFourphp_Common::$epsln) {
                        $lat = 0.;
                    } else {
                        $lat = asin($y * $sinc / $rh);
                    }
                    if ($cosc != 0. || $x != 0.) {
                        $lon = atan2($x * $sinc, $cosc * $rh);
                    }
                    break;
                case self::$_obliq:
                    if (abs($rh) <= ProjFourphp_Common::$epsln) {
                        $lat = $this->phiZero;
                    } else {
                        $lat = asin($cosc * $this->sinphZero + $y * $sinc * $this->cosphZero / $rh);
                    }
                    $c   = $cosc - $this->sinphZero * sin($lat);
                    if ($c != 0. || $x != 0.) {
                        $lon = atan2($x * $sinc * $this->cosphZero, $c * $rh);
                    }
                    break;
                case self::$_nPole:
                    $y   = -$y;
                case self::$_sPole:
                    if (abs($rh) <= ProjFourphp_Common::$epsln) {
                        $lat = $this->phiZero;
                    } else {
                        $lat  = asin($this->mode == self::$_sPole ? -$cosc : $cosc);
                    }
                    $lon  = ($x == 0. && $y == 0.) ? 0. : atan2($x, $y);
                    break;
            }
            $p->x = ProjFourphp_Common::adjustLon($lon + $this->longZero);
            $p->y = $lat;
        } else {
            $rho = sqrt($x * $x + $y * $y);
            switch ($this->mode) {
                case self::$_obliq:
                case self::$_equit:
                    $tp     = 2. * atan2($rho * $this->cosXOne, $this->akmOne);
                    $cosphi = cos($tp);
                    $sinPhi = sin($tp);
                    if ($rho == 0.0) {
                        $phiL = asin($cosphi * $this->sinXOne);
                    } else {
                        $phiL = asin($cosphi * $this->sinXOne + ($y * $sinPhi * $this->cosXOne / $rho));
                    }

                    $tp    = tan(.5 * (ProjFourphp_Common::$halfPi + $phiL));
                    $x *= $sinPhi;
                    $y     = $rho * $this->cosXOne * $cosphi - $y *
                        $this->sinXOne * $sinPhi;
                    $piTwo = ProjFourphp_Common::$halfPi;
                    $halfe = .5 * $this->e;
                    break;
                case self::$_nPole:
                    $y     = -$y;
                case self::$_sPole:
                    $tp    = - $rho / $this->akmOne;
                    $phiL  = ProjFourphp_Common::$halfPi - 2. * atan($tp);
                    $piTwo = -ProjFourphp_Common::$halfPi;
                    $halfe = -.5 * $this->e;
                    break;
            }
            for ($i = static::$_niter; $i--; $phiL = $lat) { //check this
                $sinPhi = $this->e * sin($phiL);
                $lat    = 2. *
                    atan($tp * pow((1. + $sinPhi) / (1. - $sinPhi), $halfe))
                    - $piTwo;
                if (abs($phiL - $lat) < static::$_conv) {
                    if ($this->mode == self::$_sPole) $lat  = -$lat;
                    $lon  = ($x == 0. && $y == 0.) ? 0. : atan2($x, $y);
                    $p->x =
                        ProjFourphp_Common::adjustLon($lon + $this->longZero);
                    $p->y = $lat;
                    return $p;
                }
            }
        }
    }

}

ProjFourphp::$proj['stere'] = new ProjFourphp_ProjStere();