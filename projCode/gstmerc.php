<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
class ProjFourphp_ProjGstmerc
{

    public function init()
    {

        // array of:  a, b, lonZero, latZero, kZero, xZero, yZero
        $temp   = $this->b / $this->a;
        $this->e = sqrt(1.0 - $temp * $temp);
        $this->lc = $this->longZero;
        $this->rs = sqrt(1.0 + $this->e * $this->e * pow(cos($this->latZero), 4.0) / (1.0 - $this->e * $this->e));
        $sinz   = sin($this->latZero);
        $pc     = asin($sinz / $this->rs);
        $sinzpc = sin($pc);
        $this->cp = ProjFourphp_Common::latiso(0.0, $pc, $sinzpc) -
            $this->rs *
            ProjFourphp_Common::latiso($this->e, $this->latZero, $sinz);
        $this->nTwo = $this->kZero * $this->a *
            sqrt(1.0 - $this->e * $this->e) /
            (1.0 - $this->e * $this->e * $sinz * $sinz);
        $this->xs = $this->xZero;
        $this->ys = $this->yZero - $this->nTwo * $pc;

        if (!$this->title) $this->title = "Gauss Schreiber transverse mercator";
    }

    // forward equations--mapping lat,long to x,y
    // -----------------------------------------------------------------
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $l      = $this->rs * ($lon - $this->lc);
        $ls     = $this->cp + ($this->rs * ProjFourphp_Common::latiso($this->e, $lat, sin($lat)));
        $latOne = asin(sin($l) / ProjFourphp_Common:: cosh($ls));
        $lsOne  = ProjFourphp_Common::atiso(0.0, $latOne, sin($latOne));
        $p->x = $this->xs + ($this->nTwo * $lsOne);
        $p->y = $this->ys + ($this->nTwo * atan(ProjFourphp_Common::sinh($ls) / cos($l)));
        return $p;
    }

    // inverse equations--mapping x,y to lat/long
    // -----------------------------------------------------------------
    public function inverse($p)
    {

        $x = $p->x;
        $y = $p->y;

        $l      = atan(ProjFourphp_Common::sinh(($x - $this->xs) / $this->nTwo) / cos(($y - $this->ys) / $this->nTwo));
        $latOne = asin(sin(($y - $this->ys) / $this->nTwo) / ProjFourphp_Common::cosh(($x - $this->xs) / $this->nTwo));
        $lc     = ProjFourphp_Common::latiso(0.0, $latOne, sin($latOne));
        $p->x = $this->lc + $l / $this->rs;
        $p->y = ProjFourphp_Common::invlatiso($this->e, ($lc - $this->cp) / $this->rs);
        return $p;
    }

}

ProjFourphp::$proj['gstmerc'] = new ProjFourphp_ProjGstmerc();
