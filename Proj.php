<?php

/**
 * Author : Julien Moquet
 * 
 * Inspired by ProjFourjs from Mike Adair madairATdmsolutions.ca
 *                      and Richard Greenwood rich@greenwoodmap.com 
 * License: LGPL as per: http://www.gnu.org/copyleft/lesser.html 
 */
class ProjFourphp_Proj
{

    /**
     * Property: readyToUse
     * Flag to indicate if initialization is complete for $this Proj object
     */
    public $readyToUse = false;

    /**
     * Property: title
     * The title to describe the projection
     */
    public $title = null;

    /**
     * Property: projName
     * The projection class for $this projection, e.g. lcc (lambert conformal conic,
     * or merc for mercator).  These are exactly equivalent to their ProjFour 
     * counterparts.
     */
    public $projName = null;

    /**
     * Property: projection
     * The projection object for $this projection. */
    public $projection = null;

    /**
     * Property: units
     * The units of the projection.  Values include 'm' and 'degrees'
     */
    public $units = null;

    /**
     * Property: datum
     * The datum specified for the projection
     */
    public $datum = null;

    /**
     * Property: xZero
     * The x coordinate origin
     */
    public $xZero = 0;

    /**
     * Property: yZero
     * The y coordinate origin
     */
    public $yZero = 0;

    /**
     * Property: localCS
     * Flag to indicate if the projection is a local one in which no transforms
     * are required.
     */
    public $localCS = false;

    /**
     *
     * @var type 
     */
    protected $_wktRE = '/^(\w+)\[(.*)\]$/';

    /**
     * Constructor: initialize
     * Constructor for ProjFourphp::Proj objects
     *
     * Parameters:
     * $srsCode - a code for map projection definition parameters.  These are usually
     * (but not always) EPSG codes.
     */
    public function __construct($srsCode)
    {

        $this->srsCodeInput = $srsCode;

        //check to see if $this is a WKT string
        if ((strpos($srsCode, 'GEOGCS') !== false) ||
                (strpos($srsCode, 'GEOCCS') !== false) ||
                (strpos($srsCode, 'PROJCS') !== false) ||
                (strpos($srsCode, 'LOCAL_CS') !== false)) {
            $this->parseWKT($srsCode);
            $this->deriveConstants();
            $this->loadProjCode($this->projName);
            return;
        }

        // DGR 2008-08-03 : support urn and url
        if (strpos($srsCode, 'urn:') === 0) {
            //urn:ORIGINATOR:def:crs:CODESPACE:VERSION:ID
            $urn = explode(':', $srsCode);
            if (($urn[1] == 'ogc' || $urn[1] == 'x-ogc') &&
                    ($urn[2] == 'def') &&
                    ($urn[3] == 'crs')) {
                $srsCode = $urn[4] . ':' . $urn[strlen($urn) - 1];
            }
        } else if (strpos($srsCode, 'http://') === 0) {
            //url#ID
            $url = explode('#', $srsCode);
            if (preg_match("/epsg.org/", $url[0])) {
                // http://www.epsg.org/#
                $srsCode = 'EPSG:' . $url[1];
            } else if (preg_match("/RIG.xml/", $url[0])) {
                //http://librairies.ign.fr/geoportail/resources/RIG.xml#
                //http://interop.ign.fr/registers/ign/RIG.xml#
                $srsCode = 'IGNF:' . $url[1];
            }
        }
        $this->srsCode = strtoupper($srsCode);
        if (strpos($this->srsCode, "EPSG") === 0) {
            $this->srsCode = $this->srsCode;
            $this->srsAuth = 'epsg';
            $this->srsProjNumber = substr($this->srsCode, 5);
            // DGR 2007-11-20 : authority IGNF
        } else if (strpos($this->srsCode, "IGNF") === 0) {
            $this->srsCode = $this->srsCode;
            $this->srsAuth = 'IGNF';
            $this->srsProjNumber = substr($this->srsCode, 5);
            // DGR 2008-06-19 : pseudo-authority CRS for WMS
        } else if (strpos($this->srsCode, "CRS") === 0) {
            $this->srsCode = $this->srsCode;
            $this->srsAuth = 'CRS';
            $this->srsProjNumber = substr($this->srsCode, 4);
        } else {
            $this->srsAuth = '';
            $this->srsProjNumber = $this->srsCode;
        }
        $this->loadProjDefinition();
    }

    /**
     * Function: loadProjDefinition
     *    Loads the coordinate system initialization string if required.
     *    Note that dynamic loading happens asynchronously so an application must 
     *    wait for the readyToUse property is set to true.
     *    To prevent dynamic loading, include the defs through a script tag in
     *    your application.
     *
     */
    public function loadProjDefinition()
    {

        //check in memory
        if (array_key_exists($this->srsCode, ProjFourphp::$defs)) {
            $this->defsLoaded();
            return;
        }
        //else check for def on the server
        $filename = dirname(__FILE__) . '/defs/' . strtoupper($this->srsAuth) . $this->srsProjNumber . '.php';

        try
        {
            ProjFourphp::loadScript($filename);
            $this->defsLoaded(); // succes
        } catch (Exception $e)
        {
            $this->loadFromService(); // fail
        }
    }

    /**
     * Function: loadFromService
     *    Creates the REST URL for loading the definition from a web service and 
     *    loads it.
     *
     *
     * DO IT AGAIN. : SHOULD PHP CODE BE GET BY WEBSERVICES ?
     */
    public function loadFromService()
    {

        //else load from web service
        $url = ProjFourphp::$defsLookupService . '/' . $this->srsAuth . '/' . $this->srsProjNumber . '/projFour/';
        try
        {
            ProjFourphp::$defs[strtoupper($this->srsAuth) . ":" . $this->srsProjNumber] = ProjFourphp::loadScript($url);
        } catch (Exception $e)
        {
            $this->defsFailed();
        }
    }

    /**
     * Function: defsLoaded
     * Continues the Proj object initilization once the def file is loaded
     *
     */
    public function defsLoaded()
    {

        $this->parseDefs();
        $this->loadProjCode($this->projName);
    }

    /**
     * Function: checkDefsLoaded
     *    $this is the loadCheck method to see if the def object exists
     *
     */
    public function checkDefsLoaded()
    {
        return isset(ProjFourphp::$defs[$this->srsCode]) && !empty(ProjFourphp::$defs[$this->srsCode]);
    }

    /**
     * Function: defsFailed
     *    Report an error in loading the defs file, but continue on using WGS84
     *
     */
    public function defsFailed()
    {
        ProjFourphp::reportError('failed to load projection definition for: ' . $this->srsCode);
        ProjFourphp::$defs[$this->srsCode] = ProjFourphp::$defs['WGS84'];  //set it to something so it can at least continue
        $this->defsLoaded();
    }

    /**
     * Function: loadProjCode
     *    Loads projection class code dynamically if required.
     *     Projection code may be included either through a script tag or in
     *     a built version of projFourphp
     *
     * An exception occurs if the projection is not found.
     */
    public function loadProjCode($projName)
    {

        if (array_key_exists($projName, ProjFourphp::$proj)) {
            $this->initTransforms();
            return;
        }
        //the filename for the projection code
        $filename = dirname(__FILE__) . '/projCode/' . $projName . '.php';

        try
        {
            ProjFourphp::loadScript($filename);
            $this->loadProjCodeSuccess($projName);
        } catch (Exception $e)
        {
            $this->loadProjCodeFailure($projName);
        }
    }

    /**
     * Function: loadProjCodeSuccess
     *    Loads any proj dependencies or continue on to final initialization.
     *
     */
    public function loadProjCodeSuccess($projName)
    {

        if (isset(ProjFourphp::$proj[$projName]->dependsOn) && !empty(ProjFourphp::$proj[$projName]->dependsOn)) {
            $this->loadProjCode(ProjFourphp::$proj[$projName]->dependsOn);
        } else {
            $this->initTransforms();
        }
    }

    /**
     * Function: defsFailed
     *    Report an error in loading the proj file.  Initialization of the Proj
     *    object has failed and the readyToUse flag will never be set.
     *
     */
    public function loadProjCodeFailure($projName)
    {
        ProjFourphp::reportError("failed to find projection file for: " . $projName);
        //TBD initialize with identity transforms so proj will still work?
    }

    /**
     * Function: checkCodeLoaded
     *    $this is the loadCheck method to see if the projection code is loaded
     *
     */
    public function checkCodeLoaded($projName)
    {

        return isset(ProjFourphp::$proj[$projName]) && !empty(ProjFourphp::$proj[$projName]);
    }

    /**
     * Function: initTransforms
     *    Finalize the initialization of the Proj object
     *
     */
    public function initTransforms()
    {

        $this->projection = clone(ProjFourphp::$proj[$this->projName]);
        ProjFourphp::extend($this->projection, $this);
        $this->init();

        // initiate depending class
        if (false !== ($dependsOn = isset($this->projection->dependsOn) && !empty($this->projection->dependsOn) ? $this->projection->dependsOn : false)) {
            ProjFourphp::extend(ProjFourphp::$proj[$dependsOn], $this->projection) &&
                    ProjFourphp::$proj[$dependsOn]->init() &&
                    ProjFourphp::extend($this->projection, ProjFourphp::$proj[$dependsOn]);
        }

        $this->readyToUse = true;
    }

    /**
     * 
     */
    public function init()
    {
        $this->projection->init();
    }

    /**
     *
     * @param type $pt
     * @return type 
     */
    public function forward($pt)
    {
        return $this->projection->forward($pt);
    }

    /**
     *
     * @param type $pt
     * @return type 
     */
    public function inverse($pt)
    {
        return $this->projection->inverse($pt);
    }

    /**
     * Function: parseWKT
     * Parses a WKT string to get initialization parameters
     *
     */
    public function parseWKT($wkt)
    {

        if (false === ($match = preg_match($this->_wktRE, $wkt, $wktMatch)))
            return;

        $wktObject = $wktMatch[1];
        $wktContent = $wktMatch[2];
        $wktTemp = explode(",", $wktContent);

        $wktName = (strtoupper($wktObject) == "TOWGS84") ? "TOWGS84" : array_shift($wktTemp);
        $wktName = preg_replace('/^\"/', "", $wktName);
        $wktName = preg_replace('/\"$/', "", $wktName);

        /*
          $wktContent = implode(",",$wktTemp);
          $wktArray = explode("],",$wktContent);
          for ($i=0; i<sizeof($wktArray)-1; ++$i) {
          $wktArray[$i] .= "]";
          }
         */

        $wktArray = array();
        $bkCount = 0;
        $obj = "";

        foreach ($wktTemp as $token) {

            $bkCount = substr_count($token, "[") - substr_count($token, "]");

            // ???
            $obj .= $token;
            if ($bkCount === 0) {
                array_push($wktArray, $obj);
                $obj = "";
            } else {
                $obj .= ",";
            }
        }

        //do something based on the type of the wktObject being parsed
        //add in variations in the spelling as required
        switch ($wktObject) {
            case 'LOCAL_CS':
                $this->projName = 'identity';
                $this->localCS = true;
                $this->srsCode = $wktName;
                break;
            case 'GEOGCS':
                $this->projName = 'longlat';
                $this->geocsCode = $wktName;
                if (!$this->srsCode)
                    $this->srsCode = $wktName;
                break;
            case 'PROJCS':
                $$this->srsCode = $wktName;
                break;
            case 'GEOCCS':
                break;
            case 'PROJECTION':
                $this->projName = ProjFourphp::$wktProjections[$wktName];
                break;
            case 'DATUM':
                $this->datumName = $wktName;
                break;
            case 'LOCAL_DATUM':
                $this->datumCode = 'none';
                break;
            case 'SPHEROID':
                $this->ellps = $wktName;
                $this->a = floatval(array_shift($wktArray));
                $this->rf = floatval(array_shift($wktArray));
                break;
            case 'PRIMEM':
                $this->from_greenwich = floatval(array_shift($wktArray)); //to radians?
                break;
            case 'UNIT':
                $this->units = $wktName;
                $this->unitsPerMeter = floatval(array_shift($wktArray));
                break;
            case 'PARAMETER':
                $name = strtolower($wktName);
                $value = floatval(array_shift($wktArray));
                //there may be many variations on the wktName values, add in case
                //statements as required
                switch ($name) {
                    case 'false_easting':
                        $this->xZero = $value;
                        break;
                    case 'false_northing':
                        $this->yZero = $value;
                        break;
                    case 'scale_factor':
                        $this->kZero = $value;
                        break;
                    case 'central_meridian':
                        $this->longZero = $value * ProjFourphp_Common::$dToR;
                        break;
                    case 'latitude_of_origin':
                        $this->latZero = $value * ProjFourphp_Common::$dToR;
                        break;
                    case 'more_here':
                        break;
                    default:
                        break;
                }
                break;
            case 'TOWGS84':
                $this->datum_params = $wktArray;
                break;
            //DGR 2010-11-12: AXIS
            case 'AXIS':
                $name = strtolower($wktName);
                $value = array_shift($wktArray);
                switch ($value) {
                    case 'EAST' : $value = 'e';
                        break;
                    case 'WEST' : $value = 'w';
                        break;
                    case 'NORTH': $value = 'n';
                        break;
                    case 'SOUTH': $value = 's';
                        break;
                    case 'UP' : $value = 'u';
                        break;
                    case 'DOWN' : $value = 'd';
                        break;
                    case 'OTHER':
                    default : $value = ' ';
                        break; //FIXME
                }
                if (!$this->axis) {
                    $this->axis = "enu";
                }
                switch ($name) {
                    case 'X': $this->axis = $value . substr($this->axis, 1, 2);
                        break;
                    case 'Y': $this->axis = substr($this->axis, 0, 1) . $value . substr($this->axis, 2, 1);
                        break;
                    case 'Z': $this->axis = substr($this->axis, 0, 2) . $value;
                        break;
                    default : break;
                }
            case 'MORE_HERE':
                break;
            default:
                break;
        }

        foreach ($wktArray as $wktArrayContent)
            $this->parseWKT($wktArrayContent);
    }

    /**
     * Function: parseDefs
     * Parses the PROJ.4 initialization string and sets the associated properties.
     *
     */
    public function parseDefs()
    {

        $this->defData = ProjFourphp::$defs[$this->srsCode];
        #$paramName;
        #$paramVal;
        if (!$this->defData) {
            return;
        }
        $paramArray = explode("+", $this->defData);
        for ($prop = 0; $prop < sizeof($paramArray); $prop++) {
            if (strlen($paramArray[$prop]) == 0)
                continue;
            $property = explode("=", $paramArray[$prop]);
            $paramName = strtolower($property[0]);
            if (sizeof($property) >= 2) {
                $paramVal = $property[1];
            }

            switch (trim($paramName)) {  // trim out spaces
                case "": break;   // throw away nameless parameter
                case "title": $this->title = $paramVal;
                    break;
                case "proj": $this->projName = trim($paramVal);
                    break;
                case "units": $this->units = trim($paramVal);
                    break;
                case "datum": $this->datumCode = trim($paramVal);
                    break;
                case "nadgrids": $this->nagrids = trim($paramVal);
                    break;
                case "ellps": $this->ellps = trim($paramVal);
                    break;
                case "a": $this->a = floatval($paramVal);
                    break;  // semi-major radius
                case "b": $this->b = floatval($paramVal);
                    break;  // semi-minor radius
                // DGR 2007-11-20
                case "rf": $this->rf = floatval(paramVal);
                    break; // inverse flattening rf= a/(a-b)
                case "lat_0": $this->latZero = $paramVal * ProjFourphp_Common::$dToR;
                    break;        // phiZero, central latitude
                case "lat_1": $this->latOne = $paramVal * ProjFourphp_Common::$dToR;
                    break;        //standard parallel 1
                case "lat_2": $this->latTwo = $paramVal * ProjFourphp_Common::$dToR;
                    break;        //standard parallel 2
                case "latTs": $this->latTs = $paramVal * ProjFourphp_Common::$dToR;
                    break;      // used in merc and eqc
                case "lon_0": $this->longZero = $paramVal * ProjFourphp_Common::$dToR;
                    break;       // lamZero, central longitude
                case "alpha": $this->alpha = floatval($paramVal) * ProjFourphp_Common::$dToR;
                    break;  //for somerc projection
                case "lonc": $this->longc = paramVal * ProjFourphp_Common::$dToR;
                    break;       //for somerc projection
                case "x_0": $this->xZero = floatval($paramVal);
                    break;  // false easting
                case "y_0": $this->yZero = floatval($paramVal);
                    break;  // false northing
                case "k_0": $this->kZero = floatval($paramVal);
                    break;  // projection scale factor
                case "k": $this->kZero = floatval($paramVal);
                    break;  // both forms returned
                case "r_a": $this->R_A = true;
                    break;                 // sphere--area of ellipsoid
                case "zone": $this->zone = intval($paramVal, 10);
                    break;  // UTM Zone
                case "south": $this->utmSouth = true;
                    break;  // UTM north/south
                case "towgs84": $this->datum_params = explode(",", $paramVal);
                    break;
                case "to_meter": $this->to_meter = floatval($paramVal);
                    break; // cartesian scaling
                case "from_greenwich": $this->from_greenwich = $paramVal * ProjFourphp_Common::$dToR;
                    break;
                // DGR 2008-07-09 : if pm is not a well-known prime meridian take
                // the value instead of 0.0, then convert to radians
                case "pm": $paramVal = trim($paramVal);
                    $this->from_greenwich = ProjFourphp::$primeMeridian[$paramVal] ? ProjFourphp::$primeMeridian[$paramVal] : floatval($paramVal);
                    $this->from_greenwich *= ProjFourphp_Common::$dToR;
                    break;
                // DGR 2010-11-12: axis
                case "axis": $paramVal = trim($paramVal);
                    $legalAxis = "ewnsud";
                    if (strlen(paramVal) == 3 &&
                            strpos($legalAxis, substr($paramVal, 0, 1)) !== false &&
                            strpos($legalAxis, substr($paramVal, 1, 1)) !== false &&
                            strpos($legalAxis, substr($paramVal, 2, 1)) !== false) {
                        $this->axis = $paramVal;
                    } //FIXME: be silent ?
                    break;
                case "no_defs": break;
                default: //alert("Unrecognized parameter: " . paramName);
            } // switch()
        } // for paramArray
        $this->deriveConstants();
    }

    /**
     * Function: deriveConstants
     * Sets several derived constant values and initialization of datum and ellipse parameters.
     *
     */
    public function deriveConstants()
    {

        if (isset($this->nagrids) && $this->nagrids == '@null')
            $this->datumCode = 'none';

        if (isset($this->datumCode) && $this->datumCode != 'none') {

            $datumDef = ProjFourphp::$datum[$this->datumCode];

            if (is_array($datumDef)) {
                $this->datum_params = array_key_exists('towgs84', $datumDef) ? explode(',', $datumDef['towgs84']) : null;
                $this->ellps = $datumDef['ellipse'];
                $this->datumName = array_key_exists('datumName', $datumDef) ? $datumDef['datumName'] : $this->datumCode;
            }
        }
        if (!isset($this->a)) {    // do we have an ellipsoid?
            if (!isset($this->ellps) || strlen($this->ellps) == 0 || !array_key_exists($this->ellps, ProjFourphp::$ellipsoid))
                $ellipse = ProjFourphp::$ellipsoid['WGS84'];
            else {
                $ellipse = ProjFourphp::$ellipsoid[$this->ellps];
            }

            ProjFourphp::extend($this, $ellipse);
        }

        if (isset($this->rf) && !isset($this->b))
            $this->b = (1.0 - 1.0 / $this->rf) * $this->a;

        if ((isset($this->rf) && $this->rf === 0) || abs($this->a - $this->b) < ProjFourphp_Common::$epsln) {
            $this->sphere = true;
            $this->b = $this->a;
        }
        $this->aTwo = $this->a * $this->a;          // used in geocentric
        $this->b2 = $this->b * $this->b;          // used in geocentric
        $this->es = ($this->aTwo - $this->b2) / $this->aTwo;  // e ^ 2
        $this->e = sqrt($this->es);        // eccentricity
        if (isset($this->R_A)) {
            $this->a *= 1. - $this->es * ( ProjFourphp_Common::$sixth + $this->es * ( ProjFourphp_Common::$ra4 + $this->es * ProjFourphp_Common::$ra6));
            $this->aTwo = $this->a * $this->a;
            $this->b2 = $this->b * $this->b;
            $this->es = 0.0;
        }
        $this->ep2 = ($this->aTwo - $this->b2) / $this->b2; // used in geocentric
        if (!isset($this->kZero))
            $this->kZero = 1.0;    //default value


            
//DGR 2010-11-12: axis
        if (!isset($this->axis)) {
            $this->axis = "enu";
        }

        $this->datum = new ProjFourphp_Datum($this);
    }

}