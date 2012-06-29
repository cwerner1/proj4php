<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by Proj4php from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
class Proj4php_ProjAeqd
{

    public function init()
    {
        $this->sinPOneTwo = sin($this->latZero);
        $this->cosPOneTwo = cos($this->latZero);
    }

    /**
     * 
     * @param type $p
     * @return type 
     */
    public function forward($p)
    {

        #$lon = $p->x;
        #$lat = $p->y;
        #$ksp;

        $sinphi = sin($p->y);
        $cosphi = cos($p->y);
        $dlon = Proj4php_Common::adjustLon(lon - $this->longZero);
        $coslon = cos($dlon);
        $g = $this->sinPOneTwo * $sinphi + $this->cosPOneTwo * $cosphi * $coslon;
        if (abs(abs($g) - 1.0) < Proj4php_Common::$epsln) {
            $ksp = 1.0;
            if ($g < 0.0) {
                Proj4php::reportError("aeqd:Fwd:PointError");
                return;
            }
        } else {
            $z = acos($g);
            $ksp = $z / sin($z);
        }
        $p->x = $this->xZero + $this->a * $ksp * $cosphi * sin($dlon);
        $p->y = $this->yZero + $this->a * $ksp * ($this->cosPOneTwo * $sinphi - $this->sinPOneTwo * $cosphi * $coslon);

        return $p;
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function inverse($p)
    {

        $p->x -= $this->xZero;
        $p->y -= $this->yZero;

        $rh = sqrt($p->x * $p->x + $p->y * $p->y);
        if ($rh > (2.0 * Proj4php_Common::$halfPi * $this->a)) {
            Proj4php::reportError("aeqdInvDataError");
            return;
        }
        $z = $rh / $this->a;

        $sinz = sin($z);
        $cosz = cos($z);

        $lon = $this->longZero;
        #$lat;
        if (abs($rh) <= Proj4php_Common::$epsln) {
            $lat = $this->latZero;
        } else {
            $lat = Proj4php_Common::asinz($cosz * $this->sinPOneTwo + ($p->y * $sinz * $this->cosPOneTwo) / $rh);
            $con = abs($this->latZero) - Proj4php_Common::$halfPi;
            if (abs($con) <= Proj4php_Common::$epsln) {
                if ($this->latZero >= 0.0) {
                    $lon = Proj4php_Common::adjustLon($this->longZero + atan2($p->x, -$p->y));
                } else {
                    $lon = Proj4php_Common::adjustLon($this->longZero - atan2(-$p->x, $p->y));
                }
            } else {
                $con = $cosz - $this->sinPOneTwo * sin($lat);
                if ((abs($con) < Proj4php_Common::$epsln) && (abs($p->x) < Proj4php_Common::$epsln)) {
                    //no-op, just keep the lon value as is
                } else {
                    #$temp = atan2( ($p->x * $sinz * $this->cosPOneTwo ), ($con * $rh ) ); // $temp is unused !?!
                    $lon = Proj4php_Common::adjustLon($this->longZero + atan2(($p->x * $sinz * $this->cosPOneTwo), ($con * $rh)));
                }
            }
        }

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

Proj4php::$proj['aeqd'] = new Proj4php_ProjAeqd();