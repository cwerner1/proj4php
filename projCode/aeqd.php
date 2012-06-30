<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
class ProjFourphp_ProjAeqd
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
        $dlon = ProjFourphp_Common::adjustLon(lon - $this->longZero);
        $coslon = cos($dlon);
        $g = $this->sinPOneTwo * $sinphi + $this->cosPOneTwo * $cosphi * $coslon;
        if (abs(abs($g) - 1.0) < ProjFourphp_Common::$epsln) {
            $ksp = 1.0;
            if ($g < 0.0) {
                ProjFourphp::reportError("aeqd:Fwd:PointError");
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
        if ($rh > (2.0 * ProjFourphp_Common::$halfPi * $this->a)) {
            ProjFourphp::reportError("aeqdInvDataError");
            return;
        }
        $z = $rh / $this->a;

        $sinz = sin($z);
        $cosz = cos($z);

        $lon = $this->longZero;
        #$lat;
        if (abs($rh) <= ProjFourphp_Common::$epsln) {
            $lat = $this->latZero;
        } else {
            $lat = ProjFourphp_Common::asinz($cosz * $this->sinPOneTwo + ($p->y * $sinz * $this->cosPOneTwo) / $rh);
            $con = abs($this->latZero) - ProjFourphp_Common::$halfPi;
            if (abs($con) <= ProjFourphp_Common::$epsln) {
                if ($this->latZero >= 0.0) {
                    $lon = ProjFourphp_Common::adjustLon($this->longZero + atan2($p->x, -$p->y));
                } else {
                    $lon = ProjFourphp_Common::adjustLon($this->longZero - atan2(-$p->x, $p->y));
                }
            } else {
                $con = $cosz - $this->sinPOneTwo * sin($lat);
                if ((abs($con) < ProjFourphp_Common::$epsln) && (abs($p->x) < ProjFourphp_Common::$epsln)) {
                    //no-op, just keep the lon value as is
                } else {
                    #$temp = atan2( ($p->x * $sinz * $this->cosPOneTwo ), ($con * $rh ) ); // $temp is unused !?!
                    $lon = ProjFourphp_Common::adjustLon($this->longZero + atan2(($p->x * $sinz * $this->cosPOneTwo), ($con * $rh)));
                }
            }
        }

        $p->x = $lon;
        $p->y = $lat;

        return $p;
    }

}

ProjFourphp::$proj['aeqd'] = new ProjFourphp_ProjAeqd();