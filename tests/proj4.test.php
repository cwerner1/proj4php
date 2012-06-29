<?php
/*
require_once 'ProjFourphp/projFourphp.php';

$projFour = new ProjFourphp();
$projL93 = new ProjFourphp_Proj('EPSG:2154', $projFour);
$projWGS84 = new ProjFourphp_Proj('EPSG:4326', $projFour);
$projLI = new ProjFourphp_Proj('EPSG:27571', $projFour);
$projLSud = new ProjFourphp_Proj('EPSG:27563', $projFour);
$projL72 = new ProjFourphp_Proj('EPSG:31370', $projFour);

// GPS
// latitude        longitude
// 48,831938       2,355781
// 48°49'54.977''  2°21'20.812''
//
// L93
// 652709.401   6859290.946
//
// LI
// 601413.709   1125717.730

$pointSrc = new projFourphp_Point('652709.401', '6859290.946');
echo "Source : " . $pointSrc->toShortString() . " in L93 <br>";
$pointDest = $projFour->transform($projL93, $projWGS84, $pointSrc);
echo "Conversion : " . $pointDest->toShortString() . " in WGS84<br><br>";

$pointSrc = $pointDest;
echo "Source : " . $pointSrc->toShortString() . " in WGS84<br>";
$pointDest = $projFour->transform($projWGS84, $projL72, $pointSrc);
echo "Conversion : " . $pointDest->toShortString() . " in Lambert 1972<br><br>";

$pointSrc = $pointDest;
echo "Source : " . $pointSrc->toShortString() . " in Lambert 1972<br>";
$pointDest = $projFour->transform($projL72, $projWGS84, $pointSrc);
echo "Conversion : " . $pointDest->toShortString() . " in WGS84<br><br>";

$pointSrc = $pointDest;
echo "Source : " . $pointSrc->toShortString() . " in WGS84<br>";
$pointDest = $projFour->transform($projWGS84, $projLSud, $pointSrc);
echo "Conversion : " . $pointDest->toShortString() . " in Lambert Sud<br><br>";

$pointSrc = $pointDest;
echo "Source : " . $pointSrc->toShortString() . " in Lambert Sud<br>";
$pointDest = $projFour->transform($projLSud, $projWGS84, $pointSrc);
echo "Conversion : " . $pointDest->toShortString() . " in WGS84<br><br>";

$pointSrc = $pointDest;
echo "Source : " . $pointSrc->toShortString() . " in WGS84<br>";
$pointDest = $projFour->transform($projWGS84, $projLI, $pointSrc);
echo "Conversion : " . $pointDest->toShortString() . " in LI <br><br>";

$pointSrc = $pointDest;
echo "Source : " . $pointSrc->toShortString() . " in LI<br>";
$pointDest = $projFour->transform($projLI, $projL93, $pointSrc);
echo "Conversion : " . $pointDest->toShortString() . " in L93<br><br>";
 * 
 * 
 * 
 */
