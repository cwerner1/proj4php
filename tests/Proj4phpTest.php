<?php

require_once 'Proj4php/proj4php.php';

class Proj4phpTest extends PHPUnit_Framework_TestCase
{

    public function testProjFour()
    {


        $proj4 = new Proj4php();
        $projL93 = new Proj4php_Proj('EPSG:2154', $proj4);
        $projWGS84 = new Proj4php_Proj('EPSG:4326', $proj4);
        $projLI = new Proj4php_Proj('EPSG:27571', $proj4);
        $projLSud = new Proj4php_Proj('EPSG:27563', $proj4);
        $projL72 = new Proj4php_Proj('EPSG:31370', $proj4);


        $pointSrc = new proj4php_Point('652709.401', '6859290.946');
        $this->assertEquals('652709.401 6859290.946', $pointSrc->toShortString());

        $pointDest = $proj4->transform($projL93, $projWGS84, $pointSrc);
        $this->assertEquals('2.3557811127971 48.831938054369', $pointDest->toShortString());

        $pointDest = $proj4->transform($projWGS84, $projL72, $pointSrc);
        $this->assertEquals('2354.4969810662 -51359.251012595', $pointDest->toShortString());

        $pointDest = $proj4->transform($projL72, $projWGS84, $pointSrc);
        $this->assertEquals('2.3557810993491 48.831938051733', $pointDest->toShortString());

        $pointDest = $proj4->transform($projWGS84, $projLSud, $pointSrc);
        $this->assertEquals('601419.93647681 726554.08663424', $pointDest->toShortString());

        $pointDest = $proj4->transform($projLSud, $projWGS84, $pointSrc);
        $this->assertEquals('2.3557810993491 48.831938051718', $pointDest->toShortString());

        $pointDest = $proj4->transform($projWGS84, $projLI, $pointSrc);
        
        $this->assertEquals('601415.06988072 1125718.0309796', $pointDest->toShortString());
        $pointDest = $proj4->transform($projLI, $projL93, $pointSrc);

        $this->assertEquals('652709.40001126 6859290.9458141', $pointDest->toShortString());
    }

}
