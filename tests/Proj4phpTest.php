<?php

require_once 'Proj4php/proj4php.php';

/**
 * @group Proj4 
 */
class Proj4phpTest extends PHPUnit_Framework_TestCase
{

    public function testProjFour()
    {


        $projFour = new Proj4php();
        $projLNintyThree = new Proj4php_Proj('EPSG:2154', $projFour);
        $projWgsEightyFour = new Proj4php_Proj('EPSG:4326', $projFour);
        $projLI = new Proj4php_Proj('EPSG:27571', $projFour);
        $projLSud = new Proj4php_Proj('EPSG:27563', $projFour);
        $projLSeventyTwo = new Proj4php_Proj('EPSG:31370', $projFour);


        $pointSrc = new proj4php_Point('652709.401', '6859290.946');
        $this->assertEquals('652709.401 6859290.946', $pointSrc->toShortString());

        $pointDest = $projFour->transform($projLNintyThree, $projWgsEightyFour, $pointSrc);
        $this->assertEquals('2.3557811127971 48.831938054369', $pointDest->toShortString());

        $pointDest = $projFour->transform($projWgsEightyFour, $projLSeventyTwo, $pointSrc);
        $this->assertEquals('2266.955994725 -51381.924173856', $pointDest->toShortString());

        $pointDest = $projFour->transform($projLSeventyTwo, $projWgsEightyFour, $pointSrc);
        $this->assertEquals('2.3557811127971 48.831938054355', $pointDest->toShortString());

        $pointDest = $projFour->transform($projWgsEightyFour, $projLSud, $pointSrc);
        $this->assertEquals('601366.63728352 726894.08105472', $pointDest->toShortString());

        $pointDest = $projFour->transform($projLSud, $projWgsEightyFour, $pointSrc);
        $this->assertEquals('2.3557811127971 48.83193805434', $pointDest->toShortString());

        $pointDest = $projFour->transform($projWgsEightyFour, $projLI, $pointSrc);

        $this->assertEquals('601361.94628612 1126056.8581154', $pointDest->toShortString());
        
        $pointDest = $projFour->transform($projLI, $projLNintyThree, $pointSrc);

        $this->assertEquals('652709.4010008 6859290.9460975', $pointDest->toShortString());
    }

}
