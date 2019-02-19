package geodesy

import (
	"fmt"
	"math"
)

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Ordnance Survey Grid Reference functions                           (c) Chris Veness 2005-2018  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong-gridref.html                                            */
/* www.movable-type.co.uk/scripts/geodesy/docs/module-osgridref.html                              */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

//	RadToDeg for conversion from degrees to radians
//	DegToRad for conversion from radians to degrees
const (
	RadToDeg = 180 / math.Pi
	DegToRad = math.Pi / 180
)

/**
 * Convert OS grid references to/from OSGB latitude/longitude points.
 *
 * Formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is inferior
 * to Krüger as used by e.g. Karney 2011.
 *
 * www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf.
 *
 * Note OSGB grid references cover Great Britain only; Ireland and the Channel Islands have their
 * own references.
 *
 * Note that these formulae are based on ellipsoidal calculations, and according to the OS are
 * accurate to about 4–5 metres – for greater accuracy, a geoid-based transformation (OSTN15) must
 * be used.
 *
 * @module   osgridref
 * @requires latlon-ellipsoidal
 */

/*
 * Converted 2015 to work with WGS84 by default, OSGB36 as option;
 * www.ordnancesurvey.co.uk/blog/2014/12/confirmation-on-changes-to-latitude-and-longitude
 */

// OsGridRef for Ordnance Survey (OS) grid references.
//	easting - Easting in metres from OS false origin.
//	northing - Northing in metres from OS false origin.
//
//	@example
//		var grid = new OsGridRef(651409, 313177);
type OsGridRef struct {
	Easting, Northing float64
}

func NewOsGridRef(easting float64, northing float64) *OsGridRef {
	return &OsGridRef{Easting: easting, Northing: northing}
}

// LatLonToOsGrid converts latitude/longitude to Ordnance Survey grid reference easting/northing coordinate.
//
// Note formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is
// inferior to Krüger as used by e.g. Karney 2011.
//
//	@param   {LatLon}    point - latitude/longitude.
//	@returns {OsGridRef} OS Grid Reference easting/northing.
//
//	@example
//		var p = new LatLon(52.65798, 1.71605);
//		var grid = OsGridRef.latLonToOsGrid(p); // grid.toString(): TG 51409 13177
//		// for conversion of (historical) OSGB36 latitude/longitude point:
//		var p = new LatLon(52.65757, 1.71791, OSGB36);
func (point *LatLon) LatLonToOsGrid() OsGridRef {
	// if necessary convert to OSGB36 first
	if point.Datum != OSGB36 {
		point = point.convertDatum(OSGB36)
	}

	φ := point.Lat * DegToRad
	λ := point.Lon * DegToRad

	a, b := 6377563.396, 6356256.909     // Airy 1830 major & minor semi-axes
	F0 := 0.9996012717                   // NatGrid scale factor on central meridian
	φ0, λ0 := 49*DegToRad, (-2)*DegToRad // NatGrid true origin is 49°N,2°W
	N0, E0 := -100000.0, 400000.0        // northing & easting of true origin, metres
	e2 := 1 - (b*b)/(a*a)                // eccentricity squared
	n := (a - b) / (a + b)
	n2 := n * n
	n3 := n * n * n // n, n², n³

	cosφ, sinφ := math.Cos(φ), math.Sin(φ)
	ν := a * F0 / math.Sqrt(1-e2*sinφ*sinφ)                // nu = transverse radius of curvature
	ρ := a * F0 * (1 - e2) / math.Pow(1-e2*sinφ*sinφ, 1.5) // rho = meridional radius of curvature
	η2 := ν/ρ - 1                                          // eta = ?

	Ma := (1 + n + (5.0/4.0)*n2 + (5.0/4.0)*n3) * (φ - φ0)
	Mb := (3*n + 3*n*n + (21.0/8.0)*n3) * math.Sin(φ-φ0) * math.Cos(φ+φ0)
	Mc := ((15.0/8.0)*n2 + (15.0/8.0)*n3) * math.Sin(2*(φ-φ0)) * math.Cos(2*(φ+φ0))
	Md := (35.0 / 24.0) * n3 * math.Sin(3*(φ-φ0)) * math.Cos(3*(φ+φ0))
	M := b * F0 * (Ma - Mb + Mc - Md) // meridional arc

	cos3φ := cosφ * cosφ * cosφ
	cos5φ := cos3φ * cosφ * cosφ
	tan2φ := math.Tan(φ) * math.Tan(φ)
	tan4φ := tan2φ * tan2φ

	I := M + N0
	II := (ν / 2) * sinφ * cosφ
	III := (ν / 24) * sinφ * cos3φ * (5 - tan2φ + 9*η2)
	IIIA := (ν / 720) * sinφ * cos5φ * (61 - 58*tan2φ + tan4φ)
	IV := ν * cosφ
	V := (ν / 6) * cos3φ * (ν/ρ - tan2φ)
	VI := (ν / 120) * cos5φ * (5 - 18*tan2φ + tan4φ + 14*η2 - 58*tan2φ*η2)

	Δλ := λ - λ0
	Δλ2 := Δλ * Δλ
	Δλ3 := Δλ2 * Δλ
	Δλ4 := Δλ3 * Δλ
	Δλ5 := Δλ4 * Δλ
	Δλ6 := Δλ5 * Δλ

	N := I + II*Δλ2 + III*Δλ4 + IIIA*Δλ6
	E := E0 + IV*Δλ + V*Δλ3 + VI*Δλ5

	// TODO:
	//N = Number(N.toFixed(3)) // round to mm precision
	//E = Number(E.toFixed(3))

	return OsGridRef{E, N} // gets truncated to SW corner of 1m grid square
}

// Converts Ordnance Survey grid reference easting/northing coordinate to latitude/longitude
// (SW corner of grid square).
//
// Note formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is
// inferior to Krüger as used by e.g. Karney 2011.
//
//	@param   {OsGridRef}    gridref - Grid ref E/N to be converted to lat/long (SW corner of grid square).
//	@param   {LatLon.datum} [datum=WGS84] - Datum to convert grid reference into.
//	@returns {LatLon}       Latitude/longitude of supplied grid reference.
//
//	@example
//		var gridref = new OsGridRef(651409.903, 313177.270);
//		var pWgs84 = OsGridRef.osGridToLatLon(gridref);            // 52°39′28.723″N, 001°42′57.787″E
//		// to obtain (historical) OSGB36 latitude/longitude point:
//		var pOsgb = OsGridRef.osGridToLatLon(gridref, OSGB36);     // 52°39′27.253″N, 001°43′04.518″E
func (gridRef *OsGridRef) OsGridToLatLon(datum *Datum) *LatLon {
	if datum == nil {
		datum = WGS84
	}
	E := gridRef.Easting
	N := gridRef.Northing

	a, b := 6377563.396, 6356256.909       // Airy 1830 major & minor semi-axes
	F0 := 0.9996012717                     // NatGrid scale factor on central meridian
	φ0, λ0 := 49.0*DegToRad, -2.0*DegToRad // NatGrid true origin is 49°N,2°W
	N0, E0 := -100000.0, 400000.0          // northing & easting of true origin, metres
	e2 := 1.0 - (b*b)/(a*a)                // eccentricity squared
	n := (a - b) / (a + b)
	n2 := n * n
	n3 := n * n * n // n, n², n³

	φ := φ0
	var M float64 = 0
	for {
		φ = (N-N0-M)/(a*F0) + φ

		Ma := (1.0 + n + (5.0/4.0)*n2 + (5.0/4.0)*n3) * (φ - φ0)
		Mb := (3.0*n + 3.0*n*n + (21.0/8.0)*n3) * math.Sin(φ-φ0) * math.Cos(φ+φ0)
		Mc := ((15.0/8.0)*n2 + (15.0/8.0)*n3) * math.Sin(2.0*(φ-φ0)) * math.Cos(2.0*(φ+φ0))
		Md := (35.0 / 24.0) * n3 * math.Sin(3.0*(φ-φ0)) * math.Cos(3.0*(φ+φ0))
		M = b * F0 * (Ma - Mb + Mc - Md) // meridional arc

		if N-N0-M < 0.00001 {
			break // ie until < 0.01mm
		}
	}

	cosφ, sinφ := math.Cos(φ), math.Sin(φ)
	ν := a * F0 / math.Sqrt(1-e2*sinφ*sinφ)                    // nu = transverse radius of curvature
	ρ := a * F0 * (1.0 - e2) / math.Pow(1.0-e2*sinφ*sinφ, 1.5) // rho = meridional radius of curvature
	η2 := ν/ρ - 1.0                                            // eta = ?

	tanφ := math.Tan(φ)
	tan2φ := tanφ * tanφ
	tan4φ := tan2φ * tan2φ
	tan6φ := tan4φ * tan2φ
	secφ := 1.0 / cosφ
	ν3 := ν * ν * ν
	ν5 := ν3 * ν * ν
	ν7 := ν5 * ν * ν
	VII := tanφ / (2.0 * ρ * ν)
	VIII := tanφ / (24.0 * ρ * ν3) * (5.0 + 3.0*tan2φ + η2 - 9.0*tan2φ*η2)
	IX := tanφ / (720.0 * ρ * ν5) * (61.0 + 90.0*tan2φ + 45.0*tan4φ)
	X := secφ / ν
	XI := secφ / (6.0 * ν3) * (ν/ρ + 2.0*tan2φ)
	XII := secφ / (120.0 * ν5) * (5.0 + 28.0*tan2φ + 24.0*tan4φ)
	XIIA := secφ / (5040.0 * ν7) * (61.0 + 662.0*tan2φ + 1320.0*tan4φ + 720.0*tan6φ)

	dE := E - E0
	dE2 := dE * dE
	dE3 := dE2 * dE
	dE4 := dE2 * dE2
	dE5 := dE3 * dE2
	dE6 := dE4 * dE2
	dE7 := dE5 * dE2
	φ = φ - VII*dE2 + VIII*dE4 - IX*dE6
	λ := λ0 + X*dE - XI*dE3 + XII*dE5 - XIIA*dE7

	point := &LatLon{Lat: φ * RadToDeg, Lon: λ * RadToDeg, Datum: OSGB36}

	if datum != OSGB36 {
		point = point.convertDatum(datum)
	}
	return point
}

func toNumeric(value float64) string {
	intPart, floatPart := math.Modf(value)
	result := fmt.Sprintf("%06d", int(intPart))
	if floatPart > 0 {
		result += fmt.Sprintf("%.3f", floatPart)[1:]
	}
	return result
}

/**
 * Converts ‘gridRef’ numeric grid reference to standard OS grid reference.
 *
 * @param   {number} [digits=10] - Precision of returned grid reference (10 digits = metres);
 *   digits=0 will return grid reference in numeric format.
 * @returns {string} This grid reference in standard format.
 *
 * @example
 *   var ref = new OsGridRef(651409, 313177).toString(); // TG 51409 13177
 */
func (gridRef *OsGridRef) ToString(digits uint) string {
	e := gridRef.Easting
	n := gridRef.Northing

	if math.IsNaN(e) || math.IsNaN(n) {
		panic("Invalid grid reference")
	}

	// use digits = 0 to return numeric format (in metres, allowing for decimals & for northing > 1e6)
	if digits == 0 {
		return fmt.Sprintf("%s,%s", toNumeric(e), toNumeric(n))
	} else {
		panic("code not written")
	}

	//// get the 100km-grid indices
	//var e100k = math.Floor(e/100000), n100k = math.Floor(n/100000);
	//
	//if (e100k<0 || e100k>6 || n100k<0 || n100k>12) return '';
	//
	//// translate those into numeric equivalents of the grid letters
	//var l1 = (19-n100k) - (19-n100k)%5 + math.Floor((e100k+10)/5);
	//var l2 = (19-n100k)*5%25 + e100k%5;
	//
	//// compensate for skipped 'I' and build grid letter-pairs
	//if (l1 > 7) l1++;
	//if (l2 > 7) l2++;
	//var letterPair = String.fromCharCode(l1+'A'.charCodeAt(0), l2+'A'.charCodeAt(0));
	//
	//// strip 100km-grid indices from easting & northing, and reduce precision
	//e = math.Floor((e%100000)/math.Pow(10, 5-digits/2));
	//n = math.Floor((n%100000)/math.Pow(10, 5-digits/2));
	//
	//// pad eastings & northings with leading zeros (just in case, allow up to 16-digit (mm) refs)
	//e = ('00000000'+e).slice(-digits/2);
	//n = ('00000000'+n).slice(-digits/2);
	//
	//return letterPair + ' ' + e + ' ' + n;
}
