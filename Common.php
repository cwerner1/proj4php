<?php
/**
 * @package Proj4
 */

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourjs from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodmap.com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
class ProjFourphp_Common
{

    public static $pi        = 3.141592653589793238; //Math.PI,
    public static $halfPi    = 1.570796326794896619; //Math.PI*0.5,
    public static $twoPi     = 6.283185307179586477; //Math.PI*2,
    public static $fortPi    = 0.78539816339744833;
    public static $rToD      = 57.29577951308232088;
    public static $dToR      = 0.01745329251994329577;
    public static $secToRad  = 4.84813681109535993589914102357e-6;
    /* SEC_TO_RAD = Pi/180/3600 */
    public static $epsln     = 1.0e-10;
    public static $maxIter   = 20;
    // following constants from geocent.c
    public static $cosOf67P5 = 0.38268343236508977;  /* cosine of 67.5 degrees */
    public static $adC       = 1.0026000;
    /* Toms region 1 constant */

    /* datum_type values */
    public static $pjdUnknown      = 0;
    public static $pjdThreeParam   = 1;
    public static $pjdSevenParam   = 2;
    public static $pjdGridshift    = 3;
    public static $pjdWgsEightFour = 4;   // WGS84 or equivalent
    public static $pjdNodatum      = 5;   // WGS84 or equivalent

    const SRS_WGS84_SEMIMAJOR = 6378137.0;  // only used in grid shift 

    //transforms
    // ellipoid pj_set_ell.c

    public static $sixth  = .1666666666666666667; /* 1/6 */
    public static $raFour = .04722222222222222222; /* 17/360 */
    public static $raSix  = .02215608465608465608; /* 67/3024 */
    public static $rvFour = .06944444444444444444; /* 5/72 */
    public static $rvSix  = .04243827160493827160; /* 55/1296 */


    /* meridinal distance for ellipsoid and inverse
     * *	8th degree - accurate to < 1e-5 meters when used in conjuction
     * *		with typical major axis values.
     * *	Inverse determines phi to EPS (1e-11) radians, about 1e-6 seconds.
     */
    protected static $_cZeroZero   = 1.0;
    protected static $_cZeroTwo    = .25;
    protected static $_cZeroFour   = .046875;
    protected static $_cZeroSix    = .01953125;
    protected static $_cZeroEight  = .01068115234375;
    protected static $_cTwoTwo     = .75;
    protected static $_cFourFour   = .46875;
    protected static $_cFourSix    = .01302083333333333333;
    protected static $_cFourEight  = .00712076822916666666;
    protected static $_cSixSix     = .36458333333333333333;
    protected static $_cSixEight   = .00569661458333333333;
    protected static $_cEightEight = .3076171875;

    /*
     * Function to compute the constant small m which is the radius of
     *   a parallel of latitude, phi, divided by the semimajor axis.
     */

    public static function msfnz($eccent, $sinphi, $cosphi)
    {
        $con = $eccent * $sinphi;
        return $cosphi / (sqrt(1.0 - $con * $con));
    }

    /**
     * Function to compute the constant small t for use in the forward
     *  computations in the Lambert Conformal Conic and the Polar
     *   Stereographic projections.
     */
    public static function tsfnz($eccent, $phi, $sinphi)
    {
        $con = $eccent * $sinphi;
        $com = 0.5 * $eccent;
        $con = pow(((1.0 - $con) / (1.0 + $con)), $com);
        return (tan(.5 * (self::$halfPi - $phi)) / $con);
    }

    /** Function to compute the latitude angle, phi2, for the inverse of the
      //   Lambert Conformal Conic and Polar Stereographic projections.
      //
      // rise up an assertion if there is no convergence.
      // ----------------------------------------------------------------
     */
    public static function phi2z($eccent, $ts)
    {
        $eccnth = .5 * $eccent;
        $phi    = self::$halfPi - 2 * atan($ts);
        for ($i      = 0; $i <= 15; $i++) {
            $con  = $eccent * sin($phi);
            $dphi = self::$halfPi -
                2 * atan($ts * (pow(((1.0 - $con) / (1.0 + $con)), $eccnth)))
                - $phi;
            $phi += $dphi;
            if (abs($dphi) <= .0000000001) return $phi;
        }
        assert("false; /* phi2z has NoConvergence */");
        return (-9999);
    }

    /* Function to compute constant small q which is the radius of a 
      parallel of latitude, phi, divided by the semimajor axis.
      ------------------------------------------------------------ */

    public static function qsfnz($eccent, $sinphi)
    {
        if ($eccent > 1.0e-7) {
            $con = $eccent * $sinphi;
            return (( 1.0 - $eccent * $eccent) * ($sinphi / (1.0 - $con * $con) - (.5 / $eccent) * log((1.0 - $con) / (1.0 + $con))));
        } else {
            return(2.0 * $sinphi);
        }
    }

    /* Function to eliminate roundoff errors in asin
      ---------------------------------------------- */

    public static function asinz($x)
    {
        if (abs($x) > 1.0) {
            $x = ($x > 1.0) ? 1.0 : -1.0;
        }
        return asin($x);
    }

// following functions from gctpc cproj.c for transverse mercator projections
    public static function eZerofn($x)
    {
        return(1.0 - 0.25 * $x * (1.0 + $x / 16.0 * (3.0 + 1.25 * $x)));
    }

    public static function eOnefn($x)
    {
        return(0.375 * $x * (1.0 + 0.25 * $x * (1.0 + 0.46875 * $x)));
    }

    public static function eTwofn($x)
    {
        return(0.05859375 * $x * $x * (1.0 + 0.75 * $x));
    }

    public static function eThreefn($x)
    {
        return($x * $x * $x * (35.0 / 3072.0));
    }

    public static function mlfn($e0, $eOne, $eTwo, $eThree, $phi)
    {
        return($e0 * $phi - $eOne * sin(2.0 * $phi) + $eTwo * sin(4.0 * $phi) - $eThree * sin(6.0 * $phi));
    }

    public static function srat($esinp, $exp)
    {
        return(pow((1.0 - $esinp) / (1.0 + $esinp), $exp));
    }

// Function to return the sign of an argument
    public static function sign($x)
    {
        if ($x < 0.0) {
            return(-1);
        } else {
            return(1);
        }
    }

// Function to adjust longitude to -180 to 180; input in radians
    public static function adjustLon($x)
    {
        $x = (abs($x) < self::$pi) ? $x : ($x - (self::$sign($x) * self::$twoPi) );
        return $x;
    }

// IGNF - DGR : algorithms used by IGN France
// Function to adjust latitude to -90 to 90; input in radians
    public static function adjustLat($x)
    {
        $x = (abs($x) < self::$halfPi) ? $x : ($x - (self::$sign($x) * self::$pi) );
        return $x;
    }

// Latitude Isometrique - close to tsfnz ...
    public static function latiso($eccent, $phi, $sinphi)
    {
        if (abs($phi) > self::$halfPi) return +NaN;
        if ($phi == self::$halfPi) return INF;
        if ($phi == -1.0 * self::$halfPi) return -1.0 * INF;

        $con = $eccent * $sinphi;
        return log(tan((self::$halfPi + $phi) / 2.0)) + $eccent * log((1.0 - $con) / (1.0 + $con)) / 2.0;
    }

    public static function fL($x, $l)
    {
        return 2.0 * atan($x * exp($l)) - self::$halfPi;
    }

// Inverse Latitude Isometrique - close to phTwoz
    public static function invlatiso($eccent, $ts)
    {
        $phi  = self::$fL(1.0, $ts);
        $iPhi = 0.0;
        $con  = 0.0;
        do {
            $iPhi = $phi;
            $con  = $eccent * sin($iPhi);
            $phi  = self::$fL(exp($eccent * log((1.0 + $con) / (1.0 - $con)) / 2.0), $ts);
        } while (abs($phi - $iPhi) > 1.0e-12);
        return $phi;
    }

    // Grande Normale
    public static function gN($a, $e, $sinphi)
    {
        $temp = $e * $sinphi;
        return $a / sqrt(1.0 - $temp * $temp);
    }

    //code from the PROJ.4 pj_mlfn.c file;  this may be useful for other projections
    public static function pjEnfn($es)
    {

        $en = array();
        $en[0] = self::$_cZeroZero - $es * (self::$_cZeroTwo + $es * (self::$_cZeroFour + $es * (self::$_cZeroSix + $es * self::$_cZeroEight)));
        $en[1] = es * (self::$_cTwoTwo - $es * (self::$_cZeroFour + $es * (self::$_cZeroSix + $es * self::$_cZeroEight)));
        $t     = $es * $es;
        $en[2] = $t * (self::$_cFourFour - $es * (self::$_cFourSix + $es * self::$_cFourEight));
        $t *= $es;
        $en[3] = $t * (self::$_cSixSix - $es * self::$_cSixEight);
        $en[4] = $t * $es * self::$_cEightEight;
        return $en;
    }

    public static function pjMlfn($phi, $sphi, $cphi, $en)
    {
        $cphi *= $sphi;
        $sphi *= $sphi;
        return($en[0] * $phi - $cphi * ($en[1] + $sphi * ($en[2] + $sphi * ($en[3] + $sphi * $en[4]))));
    }

    public static function pjInvMlfn($arg, $es, $en)
    {
        $k   = 1. / (1. - $es);
        $phi = $arg;
        for ($i   = self::$maxIter; $i; --$i) { /* rarely goes over 2 iterations */
            $s = sin($phi);
            $t = 1. - $es * $s * $s;
            //$t = self::$pj_mlfn($phi, $s, cos($phi), $en) - $arg;
            //$phi -= $t * ($t * sqrt($t)) * $k;
            $t = (self::$pjMlfn($phi, $s, cos($phi), $en) - $arg) * ($t * sqrt($t)) * $k;
            $phi -= $t;
            if (abs($t) < self::$epsln) return $phi;
        }

        ProjFourphp::reportError("cass:pjInvMlfn: Convergence error");

        return $phi;
    }

}