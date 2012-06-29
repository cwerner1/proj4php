<?php

require_once 'Proj4php/proj4php.php';

class Proj4php_CommonTest extends PHPUnit_Framework_TestCase
{

    public function testsign()
    {
        $this->assertEquals(-1, Proj4php_Common::sign(-111));
        $this->assertEquals(-1, Proj4php_Common::sign(-111.2));
        $this->assertEquals(1, Proj4php_Common::sign(1));
        $this->assertEquals(1, Proj4php_Common::sign(200));
    }

}
