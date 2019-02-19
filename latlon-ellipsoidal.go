package geodesy

import "math"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Geodesy tools for an ellipsoidal earth model                       (c) Chris Veness 2005-2016  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong-convert-coords.html                                     */
/* www.movable-type.co.uk/scripts/geodesy/docs/module-latlon-ellipsoidal.html                     */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

/**
 * Library of geodesy functions for operations on an ellipsoidal earth model.
 *
 * Includes ellipsoid parameters and datums for different coordinate systems, and methods for
 * converting between them and to cartesian coordinates.
 *
 * q.v. Ordnance Survey ‘A guide to coordinate systems in Great Britain’ Section 6
 * www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf.
 *
 */

// LatLon (polar) point with latitude & longitude values, on a specified datum.
//
//	Lat - Geodetic latitude in degrees.
//	Lon - Longitude in degrees.
//	Datum=WGS84 - Datum this point is defined within.
//
//	@example
//		var p1 = LatLon{51.4778, -0.0016, WGS84}
//
// Latitude was accurate many years before accurate timepieces
// were developed to enable the accurate measurement of longitude.
// Hence Lat/Lon. (Somewhat counter-intuitive for people used to x/y axes).
//
type LatLon struct {
	Lat, Lon float64
	Datum    *Datum
}

func NewLatLon(lat float64, lon float64) *LatLon {
	return &LatLon{Lat: lat, Lon: lon, Datum: WGS84}
}

// Ellipsoid parameters; major axis (a), minor axis (b), and flattening (f) for each ellipsoid.
type Ellipsoid struct{ a, b, f float64 }

var (
	EllipsoidWGS84         = &Ellipsoid{a: 6378137, b: 6356752.314245, f: 1 / 298.257223563}
	EllipsoidAiry1830      = &Ellipsoid{a: 6377563.396, b: 6356256.909, f: 1 / 299.3249646}
	EllipsoidAiryModified  = &Ellipsoid{a: 6377340.189, b: 6356034.448, f: 1 / 299.3249646}
	EllipsoidBessel1841    = &Ellipsoid{a: 6377397.155, b: 6356078.962818, f: 1 / 299.1528128}
	EllipsoidClarke1866    = &Ellipsoid{a: 6378206.4, b: 6356583.8, f: 1 / 294.978698214}
	EllipsoidClarke1880IGN = &Ellipsoid{a: 6378249.2, b: 6356515.0, f: 1 / 293.466021294}
	EllipsoidGRS80         = &Ellipsoid{a: 6378137, b: 6356752.314140, f: 1 / 298.257222101}
	EllipsoidIntl1924      = &Ellipsoid{a: 6378388, b: 6356911.946, f: 1 / 297} // aka Hayford
	EllipsoidWGS72         = &Ellipsoid{a: 6378135, b: 6356750.5, f: 1 / 298.26}
)

type Transform struct {
	tx, ty, tz, s, rx, ry, rz float64
}

func (t *Transform) inverse() *Transform {
	return &Transform{-t.tx, -t.ty, -t.tz, -t.s, -t.rx, -t.ry, -t.rz}
}

type Datum struct {
	ellipsoid *Ellipsoid
	transform *Transform
}

// Datums; with associated ellipsoid, and Helmert transform parameters to convert from WGS 84 into
// given datum.
//
// Note that precision of various datums will vary, and WGS-84 (original) is not defined to be
// accurate to better than ±1 metre. No transformation should be assumed to be accurate to better
// than a meter; for many datums somewhat less.
var (
	// transforms: t in metres, s in ppm, r in arcseconds
	ED50       = &Datum{ellipsoid: EllipsoidIntl1924, transform: &Transform{tx: 89.5, ty: 93.8, tz: 123.1, s: -1.2, rx: 0.0, ry: 0.0, rz: 0.156}}                     // epsg.io/1311
	Irl1975    = &Datum{ellipsoid: EllipsoidAiryModified, transform: &Transform{tx: -482.530, ty: 130.596, tz: -564.557, s: -8.150, rx: 1.042, ry: 0.214, rz: 0.631}} // epsg.io/1954
	NAD27      = &Datum{ellipsoid: EllipsoidClarke1866, transform: &Transform{tx: 8, ty: -160, tz: -176, s: 0, rx: 0, ry: 0, rz: 0}}
	NAD83      = &Datum{ellipsoid: EllipsoidGRS80, transform: &Transform{tx: 0.9956, ty: -1.9103, tz: -0.5215, s: -0.00062, rx: 0.025915, ry: 0.009426, rz: 0.011599}}
	NTF        = &Datum{ellipsoid: EllipsoidClarke1880IGN, transform: &Transform{tx: 168, ty: 60, tz: -320, s: 0, rx: 0, ry: 0, rz: 0}}
	OSGB36     = &Datum{ellipsoid: EllipsoidAiry1830, transform: &Transform{tx: -446.448, ty: 125.157, tz: -542.060, s: 20.4894, rx: -0.1502, ry: -0.2470, rz: -0.8421}} // epsg.io/1314
	Potsdam    = &Datum{ellipsoid: EllipsoidBessel1841, transform: &Transform{tx: -582, ty: -105, tz: -414, s: -8.3, rx: 1.04, ry: 0.35, rz: -3.08}}
	TokyoJapan = &Datum{ellipsoid: EllipsoidBessel1841, transform: &Transform{tx: 148, ty: -507, tz: -685, s: 0, rx: 0, ry: 0, rz: 0}}
	WGS72      = &Datum{ellipsoid: EllipsoidWGS72, transform: &Transform{tx: 0, ty: 0, tz: -4.5, s: -0.22, rx: 0, ry: 0, rz: 0.554}}
	WGS84      = &Datum{ellipsoid: EllipsoidWGS84, transform: &Transform{tx: 0.0, ty: 0.0, tz: 0.0, s: 0.0, rx: 0.0, ry: 0.0, rz: 0.0}}
)

/* sources:
 * - ED50:       www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#pon-4
 * - Irl1975:    www.osi.ie/wp-content/uploads/2015/05/transformations_booklet.pdf
 * - NAD27:      en.wikipedia.org/wiki/Helmert_transformation
 * - NAD83:      www.uvm.edu/giv/resources/WGS84_NAD83.pdf [strictly, WGS84(G1150) -> NAD83(CORS96) @ epoch 1997.0]
 *               (note NAD83(1986) ≡ WGS84(Original); confluence.qps.nl/pages/viewpage.action?pageId=29855173)
 * - NTF:        Nouvelle Triangulation Francaise geodesie.ign.fr/contenu/fichiers/Changement_systeme_geodesique.pdf
 * - OSGB36:     www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf
 * - Potsdam:    kartoweb.itc.nl/geometrics/Coordinate%20transformations/coordtrans.html
 * - TokyoJapan: www.geocachingtoolbox.com?page=datumEllipsoidDetails
 * - WGS72:      www.icao.int/safety/pbn/documentation/eurocontrol/eurocontrol wgs 84 implementation manual.pdf
 *
 * more transform parameters are available from earth-info.nga.mil/GandG/coordsys/datums/NATO_DT.pdf,
 * www.fieldenmaps.info/cconv/web/cconv_params.js
 */
/* note:
 * - ETRS89 reference frames are coincident with WGS-84 at epoch 1989.0 (ie null transform) at the one metre level.
 */

// Converts ‘point’ lat/lon coordinate to new coordinate system.
//
//	@param   {Datum} toDatum - Datum this coordinate is to be converted to.
//	@returns {LatLon} This point converted to new datum.
//
//	@example
//		var pWGS84 = new LatLon(51.4778, -0.0016, WGS84);
//		var pOSGB = pWGS84.convertDatum(OSGB36); // 51.4773°N, 000.0000°E
func (point *LatLon) convertDatum(toDatum *Datum) *LatLon {
	if point.Datum == toDatum {
		return point
	}

	var transform *Transform = nil

	if point.Datum == WGS84 {
		// converting from WGS 84
		transform = toDatum.transform
	}
	if toDatum == WGS84 {
		// converting to WGS 84; use inverse transform (don't overwrite original!)
		transform = point.Datum.transform.inverse()
	}
	if transform == nil {
		// neither point.Datum nor toDatum are WGS84: convert ‘point’ to WGS84 first
		point = point.convertDatum(WGS84)
		transform = toDatum.transform
	}

	oldCartesian := point.toCartesian()                    // convert polar to cartesian...
	newCartesian := oldCartesian.applyTransform(transform) // ...apply transform...
	newLatLon := newCartesian.toLatLonE(toDatum)           // ...and convert cartesian to polar

	return &newLatLon
}

// Converts ‘point’ point from (geodetic) latitude/longitude coordinates to
// (geocentric) cartesian (x/y/z) coordinates.
//
//	@returns {Vector3d} Vector pointing to lat/lon point, with x, y, z in metres from earth centre.
func (point *LatLon) toCartesian() vector3d {
	φ, λ := point.Lat*DegToRad, point.Lon*DegToRad
	h := 0.0 // height above ellipsoid - not currently used
	a, f := point.Datum.ellipsoid.a, point.Datum.ellipsoid.f

	sinφ, cosφ := math.Sin(φ), math.Cos(φ)
	sinλ, cosλ := math.Sin(λ), math.Cos(λ)

	eSq := 2*f - f*f                    // 1st eccentricity squared ≡ (a²-b²)/a²
	ν := a / math.Sqrt(1-eSq*sinφ*sinφ) // radius of curvature in prime vertical

	x := (ν + h) * cosφ * cosλ
	y := (ν + h) * cosφ * sinλ
	z := (ν*(1-eSq) + h) * sinφ

	return vector3d{x, y, z}
}

// Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
// degrees+minutes+seconds.
//
//	@param   {string} [format=dms] - Format point as 'd', 'dm', 'dms'.
//	@param   {number} [dp=0|2|4] - Number of decimal places to use - default 0 for dms, 2 for dm, 4 for d.
//	@returns {string} Comma-separated latitude/longitude.
func (point *LatLon) ToString(format DmsFormat, dp uint) string {
	return ToLat3(point.Lat, format, dp) + ", " + ToLon3(point.Lon, format, dp)
}

// Converts ‘v’ (geocentric) cartesian (x/y/z) point to
// (ellipsoidal geodetic) latitude/longitude coordinates on specified datum.
//
// Uses Bowring’s (1985) formulation for μm precision in concise form.
//
//	@param datum - Datum to use when converting point.
func (v vector3d) toLatLonE(datum *Datum) LatLon {
	x, y, z := v.x, v.y, v.z
	a, b, f := datum.ellipsoid.a, datum.ellipsoid.b, datum.ellipsoid.f

	e2 := 2*f - f*f           // 1st eccentricity squared ≡ (a²-b²)/a²
	ε2 := e2 / (1 - e2)       // 2nd eccentricity squared ≡ (a²-b²)/b²
	p := math.Sqrt(x*x + y*y) // distance from minor axis
	R := math.Sqrt(p*p + z*z) // polar radius

	// parametric latitude (Bowring eqn 17, replacing tanβ = z·a / p·b)
	tanβ := (b * z) / (a * p) * (1 + ε2*b/R)
	sinβ := tanβ / math.Sqrt(1+tanβ*tanβ)
	cosβ := sinβ / tanβ

	// geodetic latitude (Bowring eqn 18: tanφ = z+ε²bsin³β / p−e²cos³β)
	var φ float64
	if math.IsNaN(cosβ) {
		φ = 0.0
	} else {
		φ = math.Atan2(z+ε2*b*sinβ*sinβ*sinβ, p-e2*a*cosβ*cosβ*cosβ)
	}

	// longitude
	λ := math.Atan2(y, x)

	// height above ellipsoid (Bowring eqn 7) [not currently used]
	// sinφ := math.Sin(φ)
	// cosφ := math.Cos(φ)
	// ν := a / math.Sqrt(1-e2*sinφ*sinφ) // length of the normal terminated by the minor axis
	// h := p*cosφ + z*sinφ - (a * a / ν)

	return LatLon{Lat: φ * RadToDeg, Lon: λ * RadToDeg, Datum: datum}
}

// Applies Helmert transform to ‘v’ point using transform parameters t.
//
//	@param   t - Transform to apply to this point.
//	@returns {Vector3} Transformed point.
func (v vector3d) applyTransform(t *Transform) vector3d {
	x1 := v.x
	y1 := v.y
	z1 := v.z

	// transform parameters
	tx := t.tx                     // x-shift
	ty := t.ty                     // y-shift
	tz := t.tz                     // z-shift
	s1 := t.s/1e6 + 1              // scale: normalise parts-per-million to (s+1)
	rx := (t.rx / 3600) * DegToRad // x-rotation: normalise arcseconds to radians
	ry := (t.ry / 3600) * DegToRad // y-rotation: normalise arcseconds to radians
	rz := (t.rz / 3600) * DegToRad // z-rotation: normalise arcseconds to radians

	// apply transform
	x2 := tx + x1*s1 - y1*rz + z1*ry
	y2 := ty + x1*rz + y1*s1 - z1*rx
	z2 := tz - x1*ry + y1*rx + z1*s1

	return vector3d{x2, y2, z2}
}
