<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by Proj4js from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */

/** datum object
 */
class Proj4php_Datum
{

    public $datumType;
    public $datumParams;

    /**
     *
     * @param type $proj 
     */
    public function __construct($proj)
    {

        $this->datumType = Proj4php_Common::$pjdWgs84;   //default setting

        if (isset($proj->datumCode) && $proj->datumCode == 'none') {
            $this->datumType = Proj4php_Common::$pjdNodatum;
        }

        if (isset($proj->datumParams)) {

            for ($i = 0; $i < sizeof($proj->datumParams); $i++) {
                $proj->datumParams[$i] = floatval($proj->datumParams[$i]);
            }

            if ($proj->datumParams[0] != 0 || $proj->datumParams[1] != 0 || $proj->datumParams[2] != 0) {
                $this->datumType = Proj4php_Common::$pjdThreeParam;
            }

            if (sizeof($proj->datumParams) > 3) {
                if ($proj->datumParams[3] != 0 || $proj->datumParams[4] != 0 ||
                        $proj->datumParams[5] != 0 || $proj->datumParams[6] != 0) {

                    $this->datumType = Proj4php_Common::$pjd7Param;
                    $proj->datumParams[3] *= Proj4php_Common::$secToRad;
                    $proj->datumParams[4] *= Proj4php_Common::$secToRad;
                    $proj->datumParams[5] *= Proj4php_Common::$secToRad;
                    $proj->datumParams[6] = ($proj->datumParams[6] / 1000000.0) + 1.0;
                }
            }

            $this->datumParams = $proj->datumParams;
        }
        if (isset($proj)) {
            $this->a = $proj->a;    //datum object also uses these values
            $this->b = $proj->b;
            $this->es = $proj->es;
            $this->ep2 = $proj->ep2;
            #$this->datum_params = $proj->datum_params;
        }
    }

    /**
     *
     * @param type $dest
     * @return boolean Returns TRUE if the two datums match, otherwise FALSE.
     * @throws type 
     */
    public function compare_datums($dest)
    {

        if ($this->datumType != $dest->datumType) {
            return false; // false, datums are not equal
        } else if ($this->a != $dest->a || abs($this->es - $dest->es) > 0.000000000050) {
            // the tolerence for es is to ensure that GRS80 and WGS84
            // are considered identical
            return false;
        } else if ($this->datumType == Proj4php_Common::$pjdThreeParam) {
            return ($this->datumParams[0] == $dest->datumParams[0]
                    && $this->datumParams[1] == $dest->datumParams[1]
                    && $this->datumParams[2] == $dest->datumParams[2]);
        } else if ($this->datumType == Proj4php_Common::$pjd7Param) {
            return ($this->datumParams[0] == $dest->datumParams[0]
                    && $this->datumParams[1] == $dest->datumParams[1]
                    && $this->datumParams[2] == $dest->datumParams[2]
                    && $this->datumParams[3] == $dest->datumParams[3]
                    && $this->datumParams[4] == $dest->datumParams[4]
                    && $this->datumParams[5] == $dest->datumParams[5]
                    && $this->datumParams[6] == $dest->datumParams[6]);
        } else if ($this->datumType == Proj4php_Common::$pjdGridshift ||
                $dest->datumType == Proj4php_Common::$pjdGridshift) {
            throw(new Exception("ERROR: Grid shift transformations are not implemented."));
            return false;
        }

        return true; // datums are equal
    }

    /*
     * The function Convert_Geodetic_To_Geocentric converts geodetic coordinates
     * (latitude, longitude, and height) to geocentric coordinates (X, Y, Z),
     * according to the current ellipsoid parameters.
     *
     *    Latitude  : Geodetic latitude in radians                     (input)
     *    Longitude : Geodetic longitude in radians                    (input)
     *    Height    : Geodetic height, in meters                       (input)
     *    X         : Calculated Geocentric X coordinate, in meters    (output)
     *    Y         : Calculated Geocentric Y coordinate, in meters    (output)
     *    Z         : Calculated Geocentric Z coordinate, in meters    (output)
     *
     */

    public function GeodeticToGeocentric($p)
    {

        $longitude = $p->x;
        $latitude = $p->y;
        $height = isset($p->z) ? $p->z : 0;   //Z value not always supplied
        $errorCode = 0;  //  GEOCENT_NO_ERROR;

        /*
         * * Don't blow up if Latitude is just a little out of the value
         * * range as it may just be a rounding issue.  Also removed longitude
         * * test, it should be wrapped by cos() and sin().  NFW for PROJ.4, Sep/2001.
         */
        if ($latitude < - Proj4php_Common::$halfPi && $latitude > -1.001 * Proj4php_Common::$halfPi) {
            $latitude = - Proj4php_Common::$halfPi;
        } else if ($latitude > Proj4php_Common::$halfPi && $latitude < 1.001 * Proj4php_Common::$halfPi) {
            $latitude = Proj4php_Common::$halfPi;
        } else if (($latitude < - Proj4php_Common::$halfPi) || ($latitude > Proj4php_Common::$halfPi)) {
            /* Latitude out of range */
            Proj4php::reportError('geocent:lat out of range:' . $latitude);
            return null;
        }

        if ($longitude > Proj4php_Common::$pi)
            $longitude -= (2 * Proj4php_Common::$pi);

        $sinLat = sin($latitude); /*  sin(Latitude)  */
        $cosLat = cos($latitude); /*  cos(Latitude)  */
        $sinTwoLat = $sinLat * $sinLat; /*  Square of sin(Latitude)  */
        $rn = $this->a / (sqrt(1.0e0 - $this->es * $sinTwoLat)); /*  Earth radius at location  */
        $p->x = ($rn + $height) * $cosLat * cos($longitude);
        $p->y = ($rn + $height) * $cosLat * sin($longitude);
        $p->z = (($rn * (1 - $this->es)) + $height) * $sinLat;

        return $errorCode;
    }

    /**
     *
     * @param object $p
     * @return type 
     */
    public function geocentricToGeodetic($pC)
    {

        /* local defintions and variables */
        /* end-criterium of loop, accuracy of sin(Latitude) */
        $genau = 1.E-12;
        $genau2 = ($genau * $genau);
        $maxiter = 30;
        $x = $pC->x;
        $y = $pC->y;
        $z = $pC->z ? $pC->z : 0.0;   //Z value not always supplied

        /*
          $P;        // distance between semi-minor axis and location
          $RR;       // distance between center and location
          $CT;       // sin of geocentric latitude
          $ST;       // cos of geocentric latitude
          $RX;
          $RK;
          $RN;       // Earth radius at location
          $CPHI0;    // cos of start or old geodetic latitude in iterations
          $SPHI0;    // sin of start or old geodetic latitude in iterations
          $CPHI;     // cos of searched geodetic latitude
          $SPHI;     // sin of searched geodetic latitude
          $SDPHI;    // end-criterium: addition-theorem of sin(Latitude(iter)
         * -Latitude(iter-1))
          $At_Pole;     // indicates location is in polar region
          $iter;        // of continous iteration, max. 30 is always enough (s.a.)
          $Longitude;
          $Latitude;
          $height;
         */

        $atPole = false;
        $p = sqrt($x * $x + $y * $y);
        $rr = sqrt($x * $x + $y * $y + $z * $z);

        /*      special cases for latitude and longitude */
        if ($p / $this->a < $genau) {

            /*  special case, if P=0. (X=0., Y=0.) */
            $atPole = true;
            $longitude = 0.0;

            /*  if (X,Y,Z)=(0.,0.,0.) then Height becomes semi-minor axis
             *  of ellipsoid (=center of mass), Latitude becomes PI/2 */
            if ($rr / $this->a < $genau) {
                $latitude = Proj4php_Common::$halfPi;
                $height = -$this->b;
                return;
            }
        } else {
            /*  ellipsoidal (geodetic) longitude
             *  interval: -PI < Longitude <= +PI */
            $longitude = atan2($y, $x);
        }

        /* --------------------------------------------------------------
         * Following iterative algorithm was developped by
         * "Institut für Erdmessung", University of Hannover, July 1988.
         * Internet: www.ife.uni-hannover.de
         * Iterative computation of CPHI,SPHI and Height.
         * Iteration of CPHI and SPHI to 10**-12 radian res$p->
         * 2*10**-7 arcsec.
         * --------------------------------------------------------------
         */
        $ct = $z / $rr;
        $st = $p / $rr;
        $rx = 1.0 / sqrt(1.0 - $this->es * (2.0 - $this->es) * $st * $st);
        $cphio = $st * (1.0 - $this->es) * $rx;
        $sphio = $ct * $rx;
        $iter = 0;

        /* loop to find sin(Latitude) res$p-> Latitude
         * until |sin(Latitude(iter)-Latitude(iter-1))| < genau */
        do {
            ++$iter;
            $rn = $this->a / sqrt(1.0 - $this->es * $sphio * $sphio);

            /*  ellipsoidal (geodetic) height */
            $height = $p * $cphio + $z * $sphio - $rn * (1.0 - $this->es * $sphio * $sphio);

            $rk = $this->es * $rn / ($rn + $height);
            $rx = 1.0 / sqrt(1.0 - $rk * (2.0 - $rk) * $st * $st);
            $cphi = $st * (1.0 - $rk) * $rx;
            $sphi = $ct * $rx;
            $sdphi = $sphi * $cphio - $cphi * $sphio;
            $cphio = $cphi;
            $sphio = $sphi;
        } while ($sdphi * $sdphi > $genau2 && $iter < $maxiter);

        /*      ellipsoidal (geodetic) latitude */
        $latitude = atan($sphi / abs($cphi));

        $pC->x = $longitude;
        $pC->y = $latitude;
        $pC->z = $height;

        return $pC;
    }

    /**
     * Convert_Geocentric_To_Geodetic
     * The method used here is derived from 'An Improved Algorithm for
     * Geocentric to Geodetic Coordinate Conversion', by Ralph Toms, Feb 1996
     * 
     * @param object Point $p
     * @return object Point $p
     */
    public function geocentricToGeodeticNoniter($p)
    {

        /*
          $Longitude;
          $Latitude;
          $height;

          $W;        // distance from Z axis
          $W2;       // square of distance from Z axis
          $T0;       // initial estimate of vertical component
          $T1;       // corrected estimate of vertical component
          $S0;       // initial estimate of horizontal component
          $S1;       // corrected estimate of horizontal component
          $Sin_B0;   // sin(B0), B0 is estimate of Bowring aux variable
          $Sin3_B0;  // cube of sin(B0)
          $Cos_B0;   // cos(B0)
          $Sin_p1;   // sin(phi1), phi1 is estimated latitude
          $Cos_p1;   // cos(phi1)
          $Rn;       // Earth radius at location
          $Sum;      // numerator of cos(phi1)
          $At_Pole;  // indicates location is in polar region
         */

        $x = floatval($p->x);  // cast from string to float
        $y = floatval($p->y);
        $z = floatval($p->z ? $p->z : 0 );

        $atPole = false;
        if ($x <> 0.0) {
            $longitude = atan2($y, $x);
        } else {
            if ($y > 0) {
                $longitude = Proj4php_Common::$halfPi;
            } else if (Y < 0) {
                $longitude = - Proj4php_Common::$halfPi;
            } else {
                $atPole = true;
                $longitude = 0.0;
                if ($z > 0.0) { /* north pole */
                    $latitude = Proj4php_Common::$halfPi;
                } else if ($z< 0.0) { /* south pole */
                    $latitude = - Proj4php_Common::$halfPi;
                } else { /* center of earth */
                    $latitude = Proj4php_Common::$halfPi;
                    $height = -$this->b;
                    return;
                }
            }
        }
        $wTwo = $x * $x + $y * $y;
        $w = sqrt($wTwo);
        $to = $z * Proj4php_Common::$adC;
        $sZero = sqrt($to * $to + $wTwo);
        $sinBZero = $to / $sZero;
        $cosBzero = $w / $sZero;
        $sinThreeBZero = $sinBZero * $sinBZero * $sinBZero;
        $tOne = $z + $this->b * $this->ep2 * $sinThreeBZero;
        $sum = $w - $this->a * $this->es * $cosBzero * $cosBzero * $cosBzero;
        $sOne = sqrt($tOne * $tOne + $sum * $sum);
        $sinPOne = $tOne / $sOne;
        $cosPOne = $sum / $sOne;
        $rn = $this->a / sqrt(1.0 - $this->es * $sinPOne * $sinPOne);
        if ($cosPOne >= Proj4php_Common::$cosOf67P5) {
            $height = $w / $cosPOne - $rn;
        } else if ($cosPOne <= - Proj4php_Common::$cosOf67P5) {
            $height = $w / -$cosPOne - $rn;
        } else {
            $height = $z / $sinPOne + $rn * ($this->es - 1.0);
        }
        if ($atPole == false) {
            $latitude = atan($sinPOne / $cosPOne);
        }

        $p->x = $longitude;
        $p->y = $latitude;
        $p->z = $height;

        return $p;
    }

    /*     * ************************************************************ */

    // pj_geocentic_to_wgs84( p )
    //  p = point to transform in geocentric coordinates (x,y,z)
    public function geocentricToWgs84($p)
    {

        if ($this->datumType == Proj4php_Common::$pjdThreeParam) {
            // if( x[io] == HUGE_VAL )
            //    continue;
            $p->x += $this->datumParams[0];
            $p->y += $this->datumParams[1];
            $p->z += $this->datumParams[2];
        } else if ($this->datumType == Proj4php_Common::$pjd7Param) {
            $dxBf = $this->datumParams[0];
            $dyBf = $this->datumParams[1];
            $dzBf = $this->datumParams[2];
            $rxBf = $this->datumParams[3];
            $ryBf = $this->datumParams[4];
            $rzBf = $this->datumParams[5];
            $mBf = $this->datumParams[6];
            // if( x[io] == HUGE_VAL )
            //    continue;
            $p->x = $mBf * ( $p->x - $rzBf * $p->y + $ryBf * $p->z) + $dxBf;
            $p->y = $mBf * ( $rzBf * $p->x + $p->y - $rxBf * $p->z) + $dyBf;
            $p->z = $mBf * (-$ryBf * $p->x + $rxBf * $p->y + $p->z) + $dzBf;
        }
    }

    /*     * ************************************************************* */

    // pj_geocentic_from_wgs84()
    //  coordinate system definition,
    //  point to transform in geocentric coordinates (x,y,z)
    public function geocentric_from_wgs84($p)
    {

        if ($this->datumType == Proj4php_Common::$pjdThreeParam) {
            //if( x[io] == HUGE_VAL )
            //    continue;
            $p->x -= $this->datumParams[0];
            $p->y -= $this->datumParams[1];
            $p->z -= $this->datumParams[2];
        } else if ($this->datumType == Proj4php_Common::$pjd7Param) {
            $dxBf = $this->datumParams[0];
            $dyBf = $this->datumParams[1];
            $dzBf = $this->datumParams[2];
            $rxBf = $this->datumParams[3];
            $ryBf = $this->datumParams[4];
            $rzBf = $this->datumParams[5];
            $mBf = $this->datumParams[6];
            $xTmp = ($p->x - $dxBf) / $mBf;
            $yTmp = ($p->y - $dyBf) / $mBf;
            $zTmp = ($p->z - $dzBf) / $mBf;
            //if( x[io] == HUGE_VAL )
            //    continue;

            $p->x = $xTmp + $rzBf * $yTmp - $ryBf * $zTmp;
            $p->y = -$rzBf * $xTmp + $yTmp + $rxBf * $zTmp;
            $p->z = $ryBf * $xTmp - $rxBf * $yTmp + $zTmp;
        } //cs_geocentric_from_wgs84()
    }

}

