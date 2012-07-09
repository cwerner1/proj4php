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
/* * **************************************************************************
  NAME                  LAMBERT AZIMUTHAL EQUAL-AREA

  PURPOSE:	Transforms input longitude and latitude to Easting and
  Northing for the Lambert Azimuthal Equal-Area projection.  The
  longitude and latitude must be in radians.  The Easting
  and Northing values will be returned in meters.

  PROGRAMMER              DATE
  ----------              ----
  D. Steinwand, EROS      March, 1991

  This function was adapted from the Lambert Azimuthal Equal Area projection
  code (FORTRAN) in the General Cartographic Transformation Package software
  which is available from the U.S. Geological Survey National Mapping Division.

  ALGORITHM REFERENCES

  1.  "New Equal-Area Map Projections for Noncircular Regions", John P. Snyder,
  The American Cartographer, Vol 15, No. 4, October 1988, pp. 341-355.

  2.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
  Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
  State Government Printing Office, Washington D.C., 1987.

  3.  "Software Documentation for GCTP General Cartographic Transformation
  Package", U.S. Geological Survey National Mapping Division, May 1982.
 * ************************************************************************* */

class ProjFourphp_ProjLaea
{

    static protected $_sPole     = 1;
    static protected $_nPole     = 2;
    static protected $_equit     = 3;
    static protected $_obliq     = 4;
    static protected $_pZeroZero = .33333333333333333333;
    static protected $_pZeroOne  = .17222222222222222222;
    static protected $_pZeroTwo  = .10257936507936507936;
    static protected $_pOneZero  = .06388888888888888888;
    static protected $_pOneOne   = .06640211640211640211;
    static protected $_pTwoZero  = .01641501294219154443;

    /* Initialize the Lambert Azimuthal Equal Area projection
      ------------------------------------------------------ */

    public function init()
    {
        $t = abs($this->latZero);
        if (abs($t - ProjFourphp_Common::$halfPi) < ProjFourphp_Common::$epsln) {
            $this->mode = $this->latZero < 0. ? self::$_sPole : self::$_nPole;
        } else if (abs($t) < ProjFourphp_Common::$epsln) {
            $this->mode = self::$_equit;
        } else {
            $this->mode = self::$_obliq;
        }
        if ($this->es > 0) {
            #$sinphi;

            $this->qp = ProjFourphp::$common->qsfnz($this->e, 1.0);
            $this->mmf = .5 / (1. - $this->es);
            $this->apa = $this->authset($this->es);
            switch ($this->mode) {
                case self::$_nPole:
                case self::$_sPole:
                    $this->dd = 1.;
                    break;
                case self::$_equit:
                    $this->rq = sqrt(.5 * $this->qp);
                    $this->dd = 1. / $this->rq;
                    $this->xmf = 1.;
                    $this->ymf = .5 * $this->qp;
                    break;
                case self::$_obliq:
                    $this->rq = sqrt(.5 * $this->qp);
                    $sinphi = sin($this->latZero);
                    $this->sinbOne = ProjFourphp::$common->qsfnz($this->e, $sinphi) / $this->qp;
                    $this->cosbOne = sqrt(1. - $this->sinbOne * $this->sinbOne);
                    $this->dd = cos($this->latZero) / (sqrt(1. - $this->es * $sinphi * $sinphi) * $this->rq * $this->cosbOne);
                    $this->ymf = ($this->xmf = $this->rq) / $this->dd;
                    $this->xmf *= $this->dd;
                    break;
            }
        } else {
            if ($this->mode == self::$_obliq) {
                $this->sinphZero = sin($this->latZero);
                $this->cosphZero = cos($this->latZero);
            }
        }
    }

    /* Lambert Azimuthal Equal Area forward equations--mapping lat,long to x,y
      ----------------------------------------------------------------------- */

    public function forward($p)
    {

        /* Forward equations
          ----------------- */
        #$x;
        #$y;
        $lam = $p->x;
        $phi = $p->y;
        $lam = ProjFourphp_Common::adjustLon($lam - $this->longZero);

        if ($this->sphere) {
            /*
              $coslam;
              $cosphi;
              $sinphi;
             */

            $sinphi = sin($phi);
            $cosphi = cos($phi);
            $coslam = cos($lam);
            switch ($this->mode) {
                case self::$_obliq:
                case self::$_equit:
                    $y = ($this->mode == self::$_equit) ? 1. + $cosphi * $coslam : 1. + $this->sinphZero * $sinphi + $this->cosphZero * $cosphi * $coslam;
                    if (y <= ProjFourphp_Common::$epsln) {
                        ProjFourphp::reportError("laea:fwd:y less than eps");
                        return null;
                    }
                    $y      = sqrt(2. / $y);
                    $x      = $y * cosphi * sin($lam);
                    $y *= ($this->mode == self::$_equit) ? $sinphi : $this->cosphZero * $sinphi - $this->sinphZero * $cosphi * $coslam;
                    break;
                case self::$_nPole:
                    $coslam = -$coslam;
                case self::$_sPole:
                    if (abs($phi + $this->phiZero) < ProjFourphp_Common::$epsln) {
                        ProjFourphp::reportError("laea:fwd:phi < eps");
                        return null;
                    }
                    $y = ProjFourphp::$common->fortPi - $phi * .5;
                    $y = 2. * (($this->mode == self::$_sPole) ? cos($y) : sin($y));
                    $x = $y * sin($lam);
                    $y *= $coslam;
                    break;
            }
        } else {
            /*
              $coslam;
              $sinlam;
              $sinphi;
              $q;
             */
            $sinb = 0.0;
            $cosb = 0.0;
            $b    = 0.0;

            $coslam = cos($lam);
            $sinlam = sin($lam);
            $sinphi = sin($phi);
            $q      = ProjFourphp::$common->qsfnz($this->e, $sinphi);
            if ($this->mode == self::$_obliq || $this->mode == self::$_equit) {
                $sinb = $q / $this->qp;
                $cosb = sqrt(1. - $sinb * $sinb);
            }
            switch ($this->mode) {
                case self::$_obliq:
                    $b = 1. + $this->sinbOne * $sinb + $this->cosbOne * $cosb * $coslam;
                    break;
                case self::$_equit:
                    $b = 1. + $cosb * $coslam;
                    break;
                case self::$_nPole:
                    $b = ProjFourphp_Common::$halfPi + $phi;
                    $q = $this->qp - $q;
                    break;
                case self::$_sPole:
                    $b = $phi - ProjFourphp_Common::$halfPi;
                    $q = $this->qp + $q;
                    break;
            }
            if (abs($b) < ProjFourphp_Common::$epsln) {
                ProjFourphp::reportError("laea:fwd:b < eps");
                return null;
            }
            switch ($this->mode) {
                case self::$_obliq:
                case self::$_equit:
                    $b = sqrt(2. / $b);
                    if ($this->mode == self::$_obliq) {
                        $y = $this->ymf * $b * ($this->cosbOne * $sinb - $this->sinbOne * $cosb * $coslam);
                    } else {
                        $y = ($b = sqrt(2. / (1. + $cosb * $coslam))) * $sinb * $this->ymf;
                    }
                    $x = $this->xmf * $b * $cosb * $sinlam;
                    break;
                case self::$_nPole:
                case self::$_sPole:
                    if (q >= 0.) {
                        $x = ($b = sqrt($q)) * $sinlam;
                        $y = $coslam * (($this->mode == self::$_sPole) ? $b : -$b);
                    } else {
                        $x = $y = 0.;
                    }
                    break;
            }
        }

        //v 1.0
        /*
          $sin_lat=sin(lat);
          $cos_lat=cos(lat);

          $sin_delta_lon=sin(delta_lon);
          $cos_delta_lon=cos(delta_lon);

          $g =$this->sin_lat_o * sin_lat +$this->cos_lat_o * cos_lat * cos_delta_lon;
          if (g == -1.0) {
          ProjFourphp::reportError("laea:fwd:Point projects to a circle of radius "+ 2.0 * R);
          return null;
          }
          $ksp = $this->a * sqrt(2.0 / (1.0 + g));
          $x = ksp * cos_lat * sin_delta_lon + $this->xZero;
          $y = ksp * ($this->cos_lat_o * sin_lat - $this->sin_lat_o * cos_lat * cos_delta_lon) + $this->yZero;
         */
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
            $cosz = 0.0;
            #$rh;
            $sinz = 0.0;

            $rh  = sqrt($x * $x + $y * $y);
            $phi = $rh * .5;
            if ($phi > 1.) {
                ProjFourphp::reportError("laea:Inv:DataError");
                return null;
            }
            $phi = 2. * asin($phi);
            if ($this->mode == self::$_obliq || $this->mode == self::$_equit) {
                $sinz = sin($phi);
                $cosz = cos($phi);
            }
            switch ($this->mode) {
                case self::$_equit:
                    $phi = (abs($rh) <= ProjFourphp_Common::$epsln) ? 0. : asin($y * $sinz / $rh);
                    $x *= $sinz;
                    $y   = $cosz * $rh;
                    break;
                case self::$_obliq:
                    $phi = (abs($rh) <= ProjFourphp_Common::$epsln) ? $this->phiZero : asin($cosz * $this->sinphZero + $y * $sinz * $this->cosphZero / $rh);
                    $x *= $sinz * $this->cosphZero;
                    $y   = ($cosz - sin($phi) * $this->sinphZero) * $rh;
                    break;
                case self::$_nPole:
                    $y   = -$y;
                    $phi = ProjFourphp_Common::$halfPi - $phi;
                    break;
                case self::$_sPole:
                    $phi -= ProjFourphp_Common::$halfPi;
                    break;
            }
            $lam = ($y == 0. && ($this->mode == self::$_equit || $this->mode == self::$_obliq)) ? 0. : atan2($x, $y);
        } else {
            /*
              $cCe;
              $sCe;
              $q;
              $rho;
             */
            $ab = 0.0;

            switch ($this->mode) {
                case self::$_equit:
                case self::$_obliq:
                    $x /= $this->dd;
                    $y *= $this->dd;
                    $rho = sqrt($x * $x + $y * $y);
                    if ($rho < ProjFourphp_Common::$epsln) {
                        $p->x = 0.;
                        $p->y = $this->phiZero;
                        return $p;
                    }
                    $sCe = 2. * asin(.5 * $rho / $this->rq);
                    $cCe = cos($sCe);
                    $x *= ($sCe = sin($sCe));
                    if ($this->mode == self::$_obliq) {
                        $ab = $cCe * $this->sinbOne + $y * $sCe * $this->cosbOne / $rho;
                        $q  = $this->qp * $ab;
                        $y  = $rho * $this->cosbOne * $cCe - $y * $this->sinbOne * $sCe;
                    } else {
                        $ab = $y * $sCe / $rho;
                        $q  = $this->qp * $ab;
                        $y  = $rho * $cCe;
                    }
                    break;
                case self::$_nPole:
                    $y  = -$y;
                case self::$_sPole:
                    $q  = ($x * $x + $y * $y);
                    if (!$q) {
                        $p->x = 0.;
                        $p->y = $this->phiZero;
                        return $p;
                    }
                    /*
                      q = $this->qp - q;
                     */
                    $ab = 1. - $q / $this->qp;
                    if ($this->mode == self::$_sPole) {
                        $ab  = - $ab;
                    }
                    break;
            }
            $lam = atan2($x, $y);
            $phi = $this->authlat(asin($ab), $this->apa);
        }

        /*
          $Rh = sqrt($p->x *$p->x +$p->y * $p->y);
          $temp = Rh / (2.0 * $this->a);

          if (temp > 1) {
          ProjFourphp::reportError("laea:Inv:DataError");
          return null;
          }

          $z = 2.0 * ProjFourphp::$common.asinz(temp);
          $sin_z=sin(z);
          $cos_z=cos(z);

          $lon =$this->longZero;
          if (abs(Rh) > ProjFourphp::$common->EPSLN) {
          $lat = ProjFourphp::$common.asinz($this->sin_lat_o * cos_z +$this-> cos_lat_o * sin_z *$p->y / Rh);
          $temp =abs($this->latZero) - ProjFourphp::$common->HALF_PI;
          if (abs(temp) > ProjFourphp::$common->EPSLN) {
          temp = cos_z -$this->sin_lat_o * sin(lat);
          if(temp!=0.0) lon=ProjFourphp::$common->adjust_lon($this->longZero+atan2($p->x*sin_z*$this->cos_lat_o,temp*Rh));
          } else if ($this->latZero < 0.0) {
          lon = ProjFourphp::$common->adjust_lon($this->longZero - atan2(-$p->x,$p->y));
          } else {
          lon = ProjFourphp::$common->adjust_lon($this->longZero + atan2($p->x, -$p->y));
          }
          } else {
          lat = $this->latZero;
          }
         */
        //return(OK);
        $p->x = ProjFourphp_Common::adjustLon($this->longZero + $lam);
        $p->y = $phi;
        return $p;
    }

    /**
     * determine latitude from authalic latitude
     * 
     * @param type $es
     * @return type 
     */
    public function authset($es)
    {
        #$t;
        $apa = array();
        $apa[0] = $es * static::$_pZeroZero;
        $t      = $es * $es;
        $apa[0] += $t * static::$_pZeroOne;
        $apa[1] = $t * static::$_pOneZero;
        $t *= $es;
        $apa[0] += $t * static::$_pTwoZero;
        $apa[1] += $t * static::$_pOneOne;
        $apa[2] = $t * static::$_pTwoZero;
        return $apa;
    }

    /**
     *
     * @param type $beta
     * @param type $apa
     * @return type 
     */
    public function authlat($beta, $apa)
    {
        $t = $beta + $beta;
        return $beta + $apa[0] * sin($t) + $apa[1] * sin($t + $t) + $apa[2] *
            sin($t + $t + $t);
    }

}

ProjFourphp::$proj['laea'] = new ProjFourphp_ProjLaea();