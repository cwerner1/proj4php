<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourphp from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodma$p->com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
/* * *************************************************************************
 * ****
  NAME                       SWISS OBLIQUE MERCATOR

  PURPOSE:	Swiss projection.
  WARNING:  X and Y are inverted (weird) in the swiss coordinate system. Not
  here, since we want X to be horizontal and Y vertical.

  ALGORITHM REFERENCES
  1. "Formules et constantes pour le Calcul pour la
  projection cylindrique conforme à axe oblique et pour la transformation entre
  des systèmes de référence".
  http://www.swisstopo.admin.ch/internet/swisstopo/
 * fr/home/topics/survey/sys/refsys/switzerland.parsysrelatedOne.31216.
 * downloadList.77004.DownloadFile.tmp/swissprojectionfr.pdf

 * ************************************************************************* */

class ProjFourphp_ProjSomerc
{

    /**
     * 
     */
    public function init()
    {
        $phyZero       = $this->latZero;
        $this->lambdaZero = $this->longZero;
        $sinPhyZero    = sin($phyZero);
        $semiMajorAxis = $this->a;
        $invF          = $this->rf;
        $flattening    = 1 / $invF;
        $eTwo          = 2 * $flattening - pow($flattening, 2);
        $e             = $this->e = sqrt($eTwo);
        $this->r = $this->kZero *
            $semiMajorAxis * sqrt(1 - $eTwo) / (1 - $eTwo * pow($sinPhyZero, 2.0));
        $this->alpha = sqrt(1 + $eTwo / (1 - $eTwo) * pow(cos($phyZero), 4.0));
        $this->bZero = asin($sinPhyZero / $this->alpha);
        $this->k = log(tan(projFourphp_Common::$pi / 4.0 + $this->bZero / 2.0))
            - $this->alpha
            * log(tan(projFourphp_Common::$pi / 4.0 + $phyZero / 2.0))
            + $this->alpha
            * $e / 2
            * log((1 + $e * $sinPhyZero) / (1 - $e * $sinPhyZero));
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function forward($p)
    {
        $saOne = log(tan(projFourphp_Common::$pi / 4.0 - $p->y / 2.0));
        $saTwo = $this->e / 2.0
            * log((1 + $this->e * sin($p->y)) / (1 - $this->e * sin($p->y)));
        $s     = -$this->alpha * ($saOne + $saTwo) + $this->k;

        // spheric latitude
        $b = 2.0 * (atan(exp($s)) - projFourphp_Common::$pi / 4.0);

        // spheric longitude
        $i = $this->alpha * ($p->x - $this->lambdaZero);

        // psoeudo equatorial rotation
        $rotI = atan(sin($i) / (sin($this->bZero) * tan($b) + cos($this->bZero) * cos($i)));

        $rotB = asin(cos($this->bZero) * sin($b) - sin($this->bZero) * cos($b) * cos($i));

        $p->y = $this->r / 2.0
            * log((1 + sin($rotB)) / (1 - sin($rotB)))
            + $this->yZero;

        $p->x = $this->r * $rotI + $this->xZero;

        return $p;
    }

    /**
     *
     * @param type $p
     * @return type 
     */
    public function inverse($p)
    {

        $y = $p->x - $this->xZero;
        $x = $p->y - $this->yZero;

        $rotI = $y / $this->r;
        $rotB = 2 * (atan(exp($x / $this->r)) - projFourphp_Common::$pi / 4.0);

        $b = asin(cos($this->bZero) * sin($rotB) + sin($this->bZero) * cos($rotB) * cos($rotI));
        $i = atan(sin($rotI) / (cos($this->bZero) * cos($rotI) - sin($this->bZero) * tan($rotB)));

        $lambda = $this->lambdaZero + $i / $this->alpha;

        $s         = 0.0;
        $phy       = $b;
        $prevPhy   = -1000.0;
        $iteration = 0;
        while (abs($phy - $prevPhy) > 0.0000001) {
            if (++$iteration > 20) {
                ProjFourphp::reportError("omercFwdInfinity");
                return;
            }
            //S = log(tan(projFourphp_Common::$pi / 4.0 + phy / 2.0));
            $s       = 1.0
                / $this->alpha
                * (log(tan(projFourphp_Common::$pi / 4.0 + $b / 2.0)) - $this->k)
                + $this->e
                * log(tan(projFourphp_Common::$pi / 4.0 + asin($this->e * sin($phy)) / 2.0));
            $prevPhy = $phy;
            $phy     = 2.0 * atan(exp($s)) - projFourphp_Common::$pi / 2.0;
        }

        $p->x = $lambda;
        $p->y = $phy;

        return $p;
    }

}

ProjFourphp::$proj['somerc'] = new ProjFourphp_ProjSomerc();