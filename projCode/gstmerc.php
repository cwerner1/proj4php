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
        $temp = $this->b / $this->a;
        $this->e = sqrt(1.0 - $temp * $temp);
        $this->lc = $this->longZero;
        $this->rs = sqrt(1.0 + $this->e * $this->e * pow(cos($this->latZero), 4.0) / (1.0 - $this->e * $this->e));
        $sinz = sin($this->latZero);
        $pc = asin($sinz / $this->rs);
        $sinzpc = sin($pc);
        $this->cp = ProjFourphp::$common->latiso(0.0, $pc, $sinzpc) - $this->rs * ProjFourphp::$common->latiso($this->e, $this->latZero, $sinz);
        $this->nTwo = $this->kZero * $this->a * sqrt(1.0 - $this->e * $this->e) / (1.0 - $this->e * $this->e * $sinz * $sinz);
        $this->xs = $this->xZero;
        $this->ys = $this->yZero - $this->nTwo * $pc;

        if (!$this->title)
            $this->title = "Gauss Schreiber transverse mercator";
    }

    // forward equations--mapping lat,long to x,y
    // -----------------------------------------------------------------
    public function forward($p)
    {

        $lon = $p->x;
        $lat = $p->y;

        $L = $this->rs * ($lon - $this->lc);
        $Ls = $this->cp + ($this->rs * ProjFourphp::$common->latiso($this->e, $lat, sin($lat)));
        $latOne = asin(sin($L) / ProjFourphp::$common . cosh($Ls));
        $lsOne = ProjFourphp::$common . latiso(0.0, $latOne, sin($latOne));
        $p->x = $this->xs + ($this->nTwo * $lsOne);
        $p->y = $this->ys + ($this->nTwo * atan(ProjFourphp::$common->sinh($Ls) / cos($L)));
        return $p;
    }

    // inverse equations--mapping x,y to lat/long
    // -----------------------------------------------------------------
    public function inverse($p)
    {

        $x = $p->x;
        $y = $p->y;

        $L = atan(ProjFourphp::$common . sinh(($x - $this->xs) / $this->nTwo) / cos(($y - $this->ys) / $this->nTwo));
        $latOne = asin(sin(($y - $this->ys) / $this->nTwo) / ProjFourphp::$common . cosh(($x - $this->xs) / $this->nTwo));
        $LC = ProjFourphp::$common . latiso(0.0, $latOne, sin($latOne));
        $p->x = $this->lc + $L / $this->rs;
        $p->y = ProjFourphp::$common . invlatiso($this->e, ($LC - $this->cp) / $this->rs);
        return $p;
    }

}

ProjFourphp::$proj['gstmerc'] = new ProjFourphp_ProjGstmerc();
