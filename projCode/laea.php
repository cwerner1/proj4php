<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *****************************************************************************
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
 * ***************************************************************************** */

class ProjFourphp_ProjLaea
{

    protected $sPole = 1;
    protected $nPole = 2;
    protected $equit = 3;
    protected $obliq = 4;
    protected $P00 = .33333333333333333333;
    protected $P01 = .17222222222222222222;
    protected $P02 = .10257936507936507936;
    protected $P10 = .06388888888888888888;
    protected $P11 = .06640211640211640211;
    protected $P20 = .01641501294219154443;

    /* Initialize the Lambert Azimuthal Equal Area projection
      ------------------------------------------------------ */

    public function init()
    {
        $t = abs($this->latZero);
        if (abs($t - ProjFourphp_Common::$halfPi) < ProjFourphp_Common::$epsln) {
            $this->mode = $this->latZero < 0. ? self::$sPole : self::$nPole;
        } else if (abs($t) < ProjFourphp_Common::$epsln) {
            $this->mode = self::$equit;
        } else {
            $this->mode = self::$obliq;
        }
        if ($this->es > 0) {
            #$sinphi;

            $this->qp = ProjFourphp::$common->qsfnz($this->e, 1.0);
            $this->mmf = .5 / (1. - $this->es);
            $this->apa = $this->authset($this->es);
            switch ($this->mode) {
                case self::$nPole:
                case self::$sPole:
                    $this->dd = 1.;
                    break;
                case self::$equit:
                    $this->rq = sqrt(.5 * $this->qp);
                    $this->dd = 1. / $this->rq;
                    $this->xmf = 1.;
                    $this->ymf = .5 * $this->qp;
                    break;
                case self::$obliq:
                    $this->rq = sqrt(.5 * $this->qp);
                    $sinphi = sin($this->latZero);
                    $this->sinb1 = ProjFourphp::$common->qsfnz($this->e, $sinphi) / $this->qp;
                    $this->cosb1 = sqrt(1. - $this->sinb1 * $this->sinb1);
                    $this->dd = cos($this->latZero) / (sqrt(1. - $this->es * $sinphi * $sinphi) * $this->rq * $this->cosb1);
                    $this->ymf = ($this->xmf = $this->rq) / $this->dd;
                    $this->xmf *= $this->dd;
                    break;
            }
        } else {
            if ($this->mode == self::$obliq) {
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
                case self::$obliq:
                case self::$equit:
                    $y = ($this->mode == self::$equit) ? 1. + $cosphi * $coslam : 1. + $this->sinphZero * $sinphi + $this->cosphZero * $cosphi * $coslam;
                    if (y <= ProjFourphp_Common::$epsln) {
                        ProjFourphp::reportError("laea:fwd:y less than eps");
                        return null;
                    }
                    $y = sqrt(2. / $y);
                    $x = $y * cosphi * sin($lam);
                    $y *= ($this->mode == self::$equit) ? $sinphi : $this->cosphZero * $sinphi - $this->sinphZero * $cosphi * $coslam;
                    break;
                case self::$nPole:
                    $coslam = -$coslam;
                case self::$sPole:
                    if (abs($phi + $this->phiZero) < ProjFourphp_Common::$epsln) {
                        ProjFourphp::reportError("laea:fwd:phi < eps");
                        return null;
                    }
                    $y = ProjFourphp::$common->fortPi - $phi * .5;
                    $y = 2. * (($this->mode == self::$sPole) ? cos($y) : sin($y));
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
            $b = 0.0;

            $coslam = cos($lam);
            $sinlam = sin($lam);
            $sinphi = sin($phi);
            $q = ProjFourphp::$common->qsfnz($this->e, $sinphi);
            if ($this->mode == self::$obliq || $this->mode == self::$equit) {
                $sinb = $q / $this->qp;
                $cosb = sqrt(1. - $sinb * $sinb);
            }
            switch ($this->mode) {
                case self::$obliq:
                    $b = 1. + $this->sinb1 * $sinb + $this->cosb1 * $cosb * $coslam;
                    break;
                case self::$equit:
                    $b = 1. + $cosb * $coslam;
                    break;
                case self::$nPole:
                    $b = ProjFourphp_Common::$halfPi + $phi;
                    $q = $this->qp - $q;
                    break;
                case self::$sPole:
                    $b = $phi - ProjFourphp_Common::$halfPi;
                    $q = $this->qp + $q;
                    break;
            }
            if (abs($b) < ProjFourphp_Common::$epsln) {
                ProjFourphp::reportError("laea:fwd:b < eps");
                return null;
            }
            switch ($this->mode) {
                case self::$obliq:
                case self::$equit:
                    $b = sqrt(2. / $b);
                    if ($this->mode == self::$obliq) {
                        $y = $this->ymf * $b * ($this->cosb1 * $sinb - $this->sinb1 * $cosb * $coslam);
                    } else {
                        $y = ($b = sqrt(2. / (1. + $cosb * $coslam))) * $sinb * $this->ymf;
                    }
                    $x = $this->xmf * $b * $cosb * $sinlam;
                    break;
                case self::$nPole:
                case self::$sPole:
                    if (q >= 0.) {
                        $x = ($b = sqrt($q)) * $sinlam;
                        $y = $coslam * (($this->mode == self::$sPole) ? $b : -$b);
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

            $rh = sqrt($x * $x + $y * $y);
            $phi = $rh * .5;
            if ($phi > 1.) {
                ProjFourphp::reportError("laea:Inv:DataError");
                return null;
            }
            $phi = 2. * asin($phi);
            if ($this->mode == self::$obliq || $this->mode == self::$equit) {
                $sinz = sin($phi);
                $cosz = cos($phi);
            }
            switch ($this->mode) {
                case self::$equit:
                    $phi = (abs($rh) <= ProjFourphp_Common::$epsln) ? 0. : asin($y * $sinz / $rh);
                    $x *= $sinz;
                    $y = $cosz * $rh;
                    break;
                case self::$obliq:
                    $phi = (abs($rh) <= ProjFourphp_Common::$epsln) ? $this->phiZero : asin($cosz * $this->sinphZero + $y * $sinz * $this->cosphZero / $rh);
                    $x *= $sinz * $this->cosphZero;
                    $y = ($cosz - sin($phi) * $this->sinphZero) * $rh;
                    break;
                case self::$nPole:
                    $y = -$y;
                    $phi = ProjFourphp_Common::$halfPi - $phi;
                    break;
                case self::$sPole:
                    $phi -= ProjFourphp_Common::$halfPi;
                    break;
            }
            $lam = ($y == 0. && ($this->mode == self::$equit || $this->mode == self::$obliq)) ? 0. : atan2($x, $y);
        } else {
            /*
              $cCe;
              $sCe;
              $q;
              $rho;
             */
            $ab = 0.0;

            switch ($this->mode) {
                case self::$equit:
                case self::$obliq:
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
                    if ($this->mode == self::$obliq) {
                        $ab = $cCe * $this->sinb1 + $y * $sCe * $this->cosb1 / $rho;
                        $q = $this->qp * $ab;
                        $y = $rho * $this->cosb1 * $cCe - $y * $this->sinb1 * $sCe;
                    } else {
                        $ab = $y * $sCe / $rho;
                        $q = $this->qp * $ab;
                        $y = $rho * $cCe;
                    }
                    break;
                case self::$nPole:
                    $y = -$y;
                case self::$sPole:
                    $q = ($x * $x + $y * $y);
                    if (!$q) {
                        $p->x = 0.;
                        $p->y = $this->phiZero;
                        return $p;
                    }
                    /*
                      q = $this->qp - q;
                     */
                    $ab = 1. - $q / $this->qp;
                    if ($this->mode == self::$sPole) {
                        $ab = - $ab;
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
        $APA = array();
        $APA[0] = $es * $this->P00;
        $t = $es * $es;
        $APA[0] += $t * $this->P01;
        $APA[1] = $t * $this->P10;
        $t *= $es;
        $APA[0] += $t * $this->P02;
        $APA[1] += $t * $this->P11;
        $APA[2] = $t * $this->P20;
        return $APA;
    }

    /**
     *
     * @param type $beta
     * @param type $APA
     * @return type 
     */
    public function authlat($beta, $APA)
    {
        $t = $beta + $beta;
        return($beta + $APA[0] * sin($t) + $APA[1] * sin($t + $t) + $APA[2] * sin($t + $t + $t));
    }

}

ProjFourphp::$proj['laea'] = new ProjFourphp_ProjLaea();