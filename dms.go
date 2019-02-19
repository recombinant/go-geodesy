package geodesy

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"
)

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Geodesy representation conversion functions                        (c) Chris Veness 2002-2017  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong.html                                                    */
/* www.movable-type.co.uk/scripts/geodesy/docs/module-dms.html                                    */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

func round(num float64) int {
	return int(num + math.Copysign(0.5, num))
}

// https://stackoverflow.com/questions/18390266
func ToFixed(num float64, precision uint) float64 {
	output := math.Pow(10, float64(precision))
	return float64(round(num*output)) / output
}

type DmsFormat int

const (
	FmtDMS DmsFormat = iota
	FmtDM
	FmtD
)

/**
 * Latitude/longitude points may be represented as decimal degrees, or subdivided into sexagesimal
 * minutes and seconds.
 */

/**
 * Functions for parsing and representing degrees / minutes / seconds.
 */

// note Unicode Degree = U+00B0. Prime = U+2032, Double prime = U+2033

// Compiled regular expressions for ParseDMS
var cre1 = regexp.MustCompile(`(?i)^\p{Zs}*[+−\-]?(?P<deg>\d+(?:\.(?:\d+)?)?)[°º]?\p{Zs}?[NSEW]?\p{Zs}*$`)
var cre2 = regexp.MustCompile(`(?i)^\p{Zs}*[+−\-]?(?P<deg>\d+)(?:°\p{Zs}?|\p{Zs})(?P<min>\d+(\.(\d+)?)?)['’′]?\p{Zs}?[NSEW]?\p{Zs}*$`)
var cre3 = regexp.MustCompile(`(?i)^\p{Zs}*[+−\-]?(?P<deg>\d+)(?:[°º]\p{Zs}*|\p{Zs})(?P<min>\d+)(?:['’′]\p{Zs}?|\p{Zs})(?P<sec>\d+(?:\.(?:\d+)?)?)[”"″]?\p{Zs}?[NSEW]?\p{Zs}*$`)

// U+2212 Minus Sign
var creSign = regexp.MustCompile(`^\p{Zs}*[\-−]|[WS]\p{Zs}*$`)

// Parses string representing degrees/minutes/seconds into numeric degrees.
//
// This is very flexible on formats, allowing signed decimal degrees, or deg-min-sec optionally
// suffixed by compass direction (NSEW). A variety of separators are accepted (eg 3° 37′ 09″W).
// Seconds and minutes may be omitted.
//
// @param   {string|number} dmsStr - Degrees or deg/min/sec in variety of formats.
// @returns {number} Degrees as decimal number.
//
// @example
//     var lat = ParseDMS('51° 28′ 40.12″ N')
//     var lon = ParseDMS('000° 00′ 05.31″ W')
//     var p1 = LatLon{lat, lon} // 51.4778°N, 000.0015°W
func ParseDMS(dmsStr string) float64 {
	degrees := math.NaN()

	{
		match1 := cre1.FindStringSubmatch(dmsStr)
		if len(match1) != 0 {
			results := map[string]string{}
			for i, name := range match1 {
				results[cre1.SubexpNames()[i]] = name
			}

			degrees, _ = strconv.ParseFloat(results["deg"], 64)

			if creSign.MatchString(dmsStr) {
				degrees *= -1.0
			}
			return degrees
		}
	}

	{
		match2 := cre2.FindStringSubmatch(dmsStr)
		if len(match2) != 0 {
			results := map[string]string{}
			for i, name := range match2 {
				results[cre2.SubexpNames()[i]] = name
			}

			d, _ := strconv.Atoi(results["deg"])
			m, _ := strconv.ParseFloat(results["min"], 64)

			degrees = float64(d)/1.0 + m/60.0

			if creSign.MatchString(dmsStr) {
				degrees *= -1.0
			}
			return degrees
		}
	}

	{
		match3 := cre3.FindStringSubmatch(dmsStr)
		if len(match3) != 0 {
			results := map[string]string{}
			for i, name := range match3 {
				results[cre3.SubexpNames()[i]] = name
			}

			d, _ := strconv.Atoi(results["deg"])
			m, _ := strconv.Atoi(results["min"])
			s, _ := strconv.ParseFloat(results["sec"], 64)

			degrees = float64(d)/1.0 + float64(m)/60.0 + s/3600.0

			if creSign.MatchString(dmsStr) {
				degrees *= -1.0
			}
			return degrees
		}
	}
	return degrees
}

// Separator character to be used to separate degrees, minutes, seconds, and cardinal directions.
//
// Set to '\u202f' (narrow no-break space) for improved formatting.
//
// @example
//   var p = new LatLon(51.2, 0.33);  // 51°12′00.0″N, 000°19′48.0″E
//   Separator = '\u202f';        // narrow no-break space
//   var pʹ = new LatLon(51.2, 0.33); // 51° 12′ 00.0″ N, 000° 19′ 48.0″ E
var Separator = ""

func ToDMS1(deg float64) *string {
	return ToDMS2(deg, FmtDMS)
}

// Converts decimal degrees to deg/min/sec format
//  - degree, prime, double-prime symbols are added, but sign is discarded, though no compass
//    direction is added.
//
// @param   {number} deg - Degrees to be formatted as specified.
// @param   {string} [format=dms] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
// @param   {number} [dp=0|2|4] - Number of decimal places to use – default 0 for dms, 2 for dm, 4 for d.
// @returns {string} Degrees formatted as deg/min/secs according to specified format.
func ToDMS2(deg float64, format DmsFormat) *string {
	var dp uint
	switch format {
	case FmtD:
		dp = 4
	case FmtDM:
		dp = 2
	case FmtDMS:
	default:
		dp = 0
	}
	return ToDMS3(deg, format, dp)
}

func ToDMS3(deg float64, format DmsFormat, dp uint) *string {
	if math.IsNaN(deg) {
		return nil
	}

	deg = math.Abs(deg) // (unsigned result ready for appending compass dir'n)

	var s string
	switch format {
	case FmtD:
		degrees := ToFixed(deg, dp) // does appropriate rounding
		degrees, decimalDegrees := math.Modf(degrees)
		if dp == 0 {
			s = fmt.Sprintf("%03d°", int(degrees))
		} else {
			decimalDegrees = ToFixed(decimalDegrees, dp)
			format := fmt.Sprintf("%%.0%df", dp)
			decimalDegreesString := fmt.Sprintf(format, decimalDegrees)[1:]
			s = fmt.Sprintf("%03d%s°", int(degrees), decimalDegreesString)
		}
	case FmtDM:
		degrees := math.Floor(deg)
		minutes := ToFixed(math.Mod(deg*60, 60), dp)
		if minutes == 60.00 {
			minutes = 0.0
			degrees += 1.0
		}
		minutes, decimalMinutes := math.Modf(minutes)
		if dp == 0 {
			s = fmt.Sprintf("%03d°%s%02d′", int(degrees), Separator, int(minutes))
		} else {
			format := fmt.Sprintf("%%.0%df", dp)
			decimalMinutesString := fmt.Sprintf(format, decimalMinutes)[1:]
			s = fmt.Sprintf("%03d°%s%02d%s′", int(degrees), Separator, int(minutes), decimalMinutesString)
		}
	case FmtDMS:
		degrees := math.Floor(deg)
		minutes := math.Floor(math.Mod(deg*60.0, 60))
		seconds := ToFixed(math.Mod(deg*3600, 60), dp)
		if seconds == 60.0 {
			seconds = 0.0
			minutes += 1.0
		}
		if minutes == 60.00 {
			minutes = 0.0
			degrees += 1.0
		}
		seconds, decimalSeconds := math.Modf(seconds)
		if dp == 0 {
			s = fmt.Sprintf("%03d°%s%02d′%s%02d″", int(degrees), Separator, int(minutes), Separator, int(seconds))
		} else {
			format := fmt.Sprintf("%%.0%df", dp)
			decimalSecondsString := fmt.Sprintf(format, decimalSeconds)[1:]
			s = fmt.Sprintf("%03d°%s%02d′%s%02d%s″", int(degrees), Separator, int(minutes), Separator, int(seconds), decimalSecondsString)
		}
	default:
		panic("write some code")
	}
	return &s
}

// Converts numeric degrees to deg/min/sec latitude (2-digit degrees, suffixed with N/S).
//
// @param   {number} deg - Degrees to be formatted as specified.
// @param   {string} [format=dms] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
// @param   {number} [dp=0|2|4] - Number of decimal places to use – default 0 for dms, 2 for dm, 4 for d.
// @returns {string} Degrees formatted as deg/min/secs according to specified format.
func ToLat3(deg float64, format DmsFormat, dp uint) string {
	return _toLat(deg, ToDMS3(deg, format, dp))
}

func ToLat2(deg float64, format DmsFormat) string {
	return _toLat(deg, ToDMS2(deg, format))
}

func _toLat(deg float64, result *string) string {
	if result == nil {
		return "-"
	}

	var c string // compass direction
	if deg < 0.0 {
		c = "S"
	} else {
		c = "N"
	}
	// knock off initial '0' for lat!
	return (*result)[1:] + Separator + c
}

// Convert numeric degrees to deg/min/sec longitude (3-digit degrees, suffixed with E/W)
//
//	@param   {number} deg - Degrees to be formatted as specified.
//	@param   {string} [format=dms] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
//	@param   {number} [dp=0|2|4] - Number of decimal places to use – default 0 for dms, 2 for dm, 4 for d.
//	@returns {string} Degrees formatted as deg/min/secs according to specified format.
func ToLon3(deg float64, format DmsFormat, dp uint) string {
	return _toLon(deg, ToDMS3(deg, format, dp))
}

func ToLon2(deg float64, format DmsFormat) string {
	return _toLon(deg, ToDMS2(deg, format))
}

func _toLon(deg float64, result *string) string {
	if result == nil {
		return "-"
	}

	var c string // compass direction
	if deg < 0.0 {
		c = "W"
	} else {
		c = "E"
	}

	return *result + Separator + c
}

//
// Converts numeric degrees to deg/min/sec as a bearing (0°..360°)
//
// @param   {number} deg - Degrees to be formatted as specified.
// @param   {string} [format=dms] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
// @param   {number} [dp=0|2|4] - Number of decimal places to use – default 0 for dms, 2 for dm, 4 for d.
// @returns {string} Degrees formatted as deg/min/secs according to specified format.
//
func ToBrng(deg float64, format DmsFormat, dp uint) string {
	// normalise -ve values to 180°..360°
	// normalise to range 0..360°
	deg = math.Mod(math.Mod(deg, 360.0)+360.0, 360.0)
	bearing := ToDMS3(deg, format, dp)
	if bearing == nil {
		return "-"
	}

	if strings.HasPrefix(*bearing, "360") {
		// just in case rounding took us up to 360°!
		return "000" + strings.TrimPrefix(*bearing, "360")
	}

	return *bearing
}

type compassPrecision int

const (
	CardinalPrecision compassPrecision = iota + 1
	InterCardinalPrecision
	SecondaryInterCardinalPrecision
)

func CompassPoint1(bearing float64) string {
	return CompassPoint2(bearing, SecondaryInterCardinalPrecision)
}

// Returns compass point (to given precision) for supplied bearing.
//
// @param   {number} bearing - Bearing in degrees from north.
// @param   {number} [precision=3] - Precision (1:cardinal / 2:intercardinal / 3:secondary-intercardinal).
// @returns {string} Compass point for supplied bearing.
//
// @example
//   var point = compassPoint(24)    // point = 'NNE'
//   var point = compassPoint(24, 1) // point = 'N'
func CompassPoint2(bearing float64, precision compassPrecision) string {
	// note precision could be extended to 4 for quarter-winds (eg NbNW), but I think they are little used

	// normalise to range 0..360°
	bearing = math.Mod(math.Mod(bearing, 360.0)+360.0, 360.0)

	var cardinals = [...]string{
		"N", "NNE", "NE", "ENE",
		"E", "ESE", "SE", "SSE",
		"S", "SSW", "SW", "WSW",
		"W", "WNW", "NW", "NNW",
		"N"} // the extra "N" eliminates a modulo n for idx below.

	// number of compass points at req’d precision
	var n float64
	switch precision {
	case CardinalPrecision:
		n = 4
	case InterCardinalPrecision:
		n = 8
	case SecondaryInterCardinalPrecision:
		n = 16
	}

	// There are 16 cardinals.
	idx := int(math.Round(bearing*n/360.0) * 16.0 / n)

	return cardinals[idx]
}
