<?php

/**
 * @group ProjFour 
 */
require_once realpath(dirname(__FILE__)) . '/../projFourphp.php';

class ProjFourphpTest
    extends PHPUnit_Framework_TestCase
{

    public function testProjFour()
    {


        $projFour          = new ProjFourphp();
        $projLNintyThree   = new ProjFourphp_Proj('EPSG:2154', $projFour);
        $projWgsEightyFour = new ProjFourphp_Proj('EPSG:4326', $projFour);
        $projLI            = new ProjFourphp_Proj('EPSG:27571', $projFour);
        $projLSud          = new ProjFourphp_Proj('EPSG:27563', $projFour);
        $projLSeventyTwo   = new ProjFourphp_Proj('EPSG:31370', $projFour);


        $pointSrc = new projFourphp_Point('652709.401', '6859290.946');
        $this->assertEquals('652709.401 6859290.946', $pointSrc->toShortString());

        $pointDest =
            $projFour->transform($projLNintyThree, $projWgsEightyFour, $pointSrc);
        $this->assertEquals('2.3557811127971 48.831938054369', $pointDest->toShortString());

        $pointDest =
            $projFour->transform($projWgsEightyFour, $projLSeventyTwo, $pointSrc);
        $this->assertEquals('2354.4969810662 -51359.251012595', $pointDest->toShortString());

        $pointDest =
            $projFour->transform($projLSeventyTwo, $projWgsEightyFour, $pointSrc);
        $this->assertEquals('2.3557810993491 48.831938051733', $pointDest->toShortString());

        $pointDest =
            $projFour->transform($projWgsEightyFour, $projLSud, $pointSrc);
        $this->assertEquals('601419.93647681 726554.08663424', $pointDest->toShortString());

        $pointDest =
            $projFour->transform($projLSud, $projWgsEightyFour, $pointSrc);
        $this->assertEquals('2.3557810993491 48.831938051718', $pointDest->toShortString());

        $pointDest =
            $projFour->transform($projWgsEightyFour, $projLI, $pointSrc);

        $this->assertEquals('601415.06988072 1125718.0309796', $pointDest->toShortString());

        $pointDest =
            $projFour->transform($projLI, $projLNintyThree, $pointSrc);

        $this->assertEquals('652709.40001126 6859290.9458141', $pointDest->toShortString());
    }

}
