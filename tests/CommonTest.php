<?php



class ProjFourphp_CommonTest extends PHPUnit_Framework_TestCase
{

    public function testsign()
    {
        $this->assertEquals(-1, ProjFourphp_Common::sign(-111));
        $this->assertEquals(-1, ProjFourphp_Common::sign(-111.2));
        $this->assertEquals(1, ProjFourphp_Common::sign(1));
        $this->assertEquals(1, ProjFourphp_Common::sign(200));
    }

}
