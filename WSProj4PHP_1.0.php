<?php

include_once("ProjFourphp/projFourphp.php");
return;

$error = false;

/**
 * Geometry-Points 
 */
if (isset($_GET['GEOM'])) {
    list($x, $y) = explode(' ', $_GET['GEOM']);
} else {
    if (isset($_GET['x'])) {
        $x = $_GET['x'];
    }
    else
        $error = true;

    if (isset($_GET['y'])) {
        $y = $_GET['y'];
    }
    else
        $error = true;
}

/**
 * Source-CRS 
 */
if (isset($_GET['SOURCECRS'])) {
    $srcProjection = str_replace('::', ':', $_GET['SOURCECRS']);
} else if (isset($_GET['projectionxy'])) {
    $srcProjection = $_GET['projectionxy'];
    $srcProjection = str_replace('::', ':', $srcProjection);
}
else
    $srcProjection = 'EPSG:2154';

/**
 * Target-CRS 
 */
if (isset($_GET['TARGETCRS'])) {
    $tgtProjection = str_replace('::', ':', $_GET['TARGETCRS']);
} else if (isset($_GET['projection'])) {
    $tgtProjection = $_GET['projection'];
    $tgtProjection = str_replace('::', ':', $tgtProjection);
}
else
    $tgtProjection = 'EPSG:4326';

/**
 * Format
 */
if (isset($_GET['format'])) {
    $format = $_GET['format'];
    if (!($format == 'xml' || $format == 'json'))
        $error = true;
}
else
    $format = 'xml';


$projFour = new ProjFourphp();
$projsource = new ProjFourphp_Proj($srcProjection, $projFour);
$projdest = new ProjFourphp_Proj($tgtProjection, $projFour);

// check the projections
if (ProjFourphp::$defs[$srcProjection] == ProjFourphp::$defs['WGS84'] && $srcProjection != 'EPSG:4326')
    $error = true;
if (ProjFourphp::$defs[$tgtProjection] == ProjFourphp::$defs['WGS84'] && $tgtProjection != 'EPSG:4326')
    $error = true;

if ($error === true) {
    if ($format == 'json') {
        echo "{\"status\":\"error\", \"erreur\": {\"code\": 2, \"message\": \"Wrong parameters.\"} }";
        exit;
    } else {
        echo "<reponse>";
        echo "  <erreur>";
        echo "    <code>2</code>";
        echo "    <message>Wrong parameters</message>";
        echo "  </erreur>";
        echo "</reponse>";
    }
}

$pointSrc = new projFourphp_Point($x, $y);
$pointDest = $projFour->transform($projsource, $projdest, $pointSrc);

$tgtProjection = str_replace(':', '::', $tgtProjection);

if ($format == 'json') {
    echo "{\"status\" :\"success\", \"point\" : {\"x\":" . $pointDest->x . ", \"y\":" . $pointDest->y . ",\"projection\" :\"" . $tgtProjection . "\"}}";
} else {
    header("Content-Type:text/xml");
    echo "<reponse>";
    echo "<point>";
    echo "<x>" . $pointDest->x . "</x>";
    echo "<y>" . $pointDest->y . "</y>";
    echo "<projection>" . $tgtProjection . "</projection>";
    echo "</point>";
    echo "</reponse>";
}