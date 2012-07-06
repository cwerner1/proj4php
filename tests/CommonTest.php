<?php

/**
 * @group Common 
 */
class ProjFourphp_CommonTest
    extends PHPUnit_Framework_TestCase
{

    public function testSign()
    {
        $this->assertEquals(-1, ProjFourphp_Common::sign(-111));
        $this->assertEquals(-1, ProjFourphp_Common::sign(-111.2));
        $this->assertEquals(1, ProjFourphp_Common::sign(1));
        $this->assertEquals(1, ProjFourphp_Common::sign(200));
    }

    public function testMsfnz()
    {
        $ret = ProjFourphp_Common::msfnz(0.12, 0.30, 0.40);
        $this->assertEquals("0.40025945221481", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::msfnz(0.2, 0.23, 0.10);
        $this->assertEquals("0.10010596820122", substr(strval($ret), 0, 16));
    }

    public function testTsfnz()
    {


        $ret = ProjFourphp_Common::tsfnz(0.12, 0.30, 0.40);
        $this->assertEquals("0.74167840619598", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::tsfnz(0.4, 0.10, 0.80);

        $this->assertEquals("1.0330253798791", substr(strval($ret), 0, 16));
    }

    public function testAsinz()
    {


        $ret = ProjFourphp_Common::asinz(10);
        $this->assertEquals("1.5707963267949", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::asinz(-100);
        $this->assertEquals("-1.5707963267949", substr(strval($ret), 0, 16));


        $ret = ProjFourphp_Common::asinz(-240);
        $this->assertEquals("-1.5707963267949", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::asinz(-370);
        $this->assertEquals("-1.5707963267949", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::asinz(310);
        $this->assertEquals("1.5707963267949", substr(strval($ret), 0, 16));
    }

    public function testeZerofn()
    {
        $ret = ProjFourphp_Common::eZerofn(0.35363122);
        $this->assertEquals("0.90486650238871", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::eZerofn(0.31245122);
        $this->assertEquals("0.91671521990135", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::eZerofn(0.1257483412);
        $this->assertEquals("0.96778286074154", substr(strval($ret), 0, 16));
    }

    public function testeOnefn()
    {
        $ret = ProjFourphp_Common::eOnefn(0.112341);
        $this->assertEquals("0.0433733525487", substr(strval($ret), 0, 15));

        $ret = ProjFourphp_Common::eOnefn(0.12141321122);
        $this->assertEquals("0.0469905908072", substr(strval($ret), 0, 15));

        $ret = ProjFourphp_Common::eOnefn(0.12544522);
        $this->assertEquals("0.04860400576082", substr(strval($ret), 0, 16));
    }

    public function testeTwofn()
    {
        $ret = ProjFourphp_Common::eTwofn(0.22253223);
        $this->assertEquals("0.00338587145", substr(strval($ret), 0, 13));

        $ret = ProjFourphp_Common::eTwofn(0.1212);
        $this->assertEquals("0.00093894785718", substr(strval($ret), 0, 16));

        $ret = ProjFourphp_Common::eTwofn(0.1422);
        $this->assertEquals("0.00131117534683", substr(strval($ret), 0, 16));
    }

}

