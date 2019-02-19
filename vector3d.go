package geodesy

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Vector handling functions                                          (c) Chris Veness 2011-2016  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/geodesy/docs/module-vector3d.html                               */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

// Library of 3-d vector manipulation routines.
//
// In a geodesy context, these vectors may be used to represent:
// - n-vector representing a normal to point on Earth's surface
// - earth-centered, earth fixed vector (≡ Gade’s ‘p-vector’)
// - great circle normal to vector (on spherical earth model)
// - motion vector on Earth's surface
// - etc
//
// Functions return vectors as return results, so that operations can be chained.
//	@example var v = v1.cross(v2).dot(v3) // ≡ v1×v2⋅v3

// Vector3d is represents a 3-d vector.
//
// The vector may be normalised, or use x/y/z values for eg height relative to the sphere or
// ellipsoid, distance from earth centre, etc.
//
//	@param x - X component of vector.
//	@param y - Y component of vector.
//	@param z - Z component of vector.
type vector3d struct {
	x, y, z float64
}
