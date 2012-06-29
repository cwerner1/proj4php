<?php

require_once 'Proj4php/proj4php.php';

/**
 * @group Proj4 
 */
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
        $this->assertEquals('2266.955994725 -51381.924173856', $pointDest->toShortString());

        $pointDest = $proj4->transform($projL72, $projWGS84, $pointSrc);
        $this->assertEquals('2.3557811127971 48.831938054355', $pointDest->toShortString());

        $pointDest = $proj4->transform($projWGS84, $projLSud, $pointSrc);
        $this->assertEquals('601366.63728352 726894.08105472', $pointDest->toShortString());

        $pointDest = $proj4->transform($projLSud, $projWGS84, $pointSrc);
        $this->assertEquals('2.3557811127971 48.83193805434', $pointDest->toShortString());

        $pointDest = $proj4->transform($projWGS84, $projLI, $pointSrc);

        $this->assertEquals('601361.94628612 1126056.8581154', $pointDest->toShortString());
        
        $pointDest = $proj4->transform($projLI, $projL93, $pointSrc);

        $this->assertEquals('652709.4010008 6859290.9460975', $pointDest->toShortString());
    }

}
