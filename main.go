package main

import (
	"encoding/json"
	"fmt"
	"io"
	"math"
	"net/http"
	"net/url"
	"os"
	"strings"
	"time"
)

// Physical constants
const (
	// Universal constants
	GravitationalConstant = 6.67430e-11    // Gravitational constant (m³/(kg·s²))
	SpeedOfLight          = 299792458      // Speed of light (m/s)
	AstronomicalUnit      = 1.495978707e11 // Astronomical unit (m)

	// Solar parameters
	SolarMass = 1.98892e30 // Mass of the Sun (kg)

	// Mercury orbital parameters
	MercurySemiMajorAxis = 0.387098 * AstronomicalUnit // Semi-major axis (m)
	MercuryEccentricity  = 0.20563                     // Eccentricity
	MercuryOrbitalPeriod = 87.969                      // Orbital period (days)

	// Conversion constants
	RadiansToArcseconds = 206264.806 // Conversion factor from radians to arcseconds
	DaysPerCentury      = 365.25 * 100

	// Observational data
	ObservedPrecession = 574.1 // Observed precession (arcseconds per century)
	ExpectedGREffect   = 43.0  // Expected GR effect (arcseconds per century)
)

// PrecessionResult holds the calculated precession values in various units and timescales.
type PrecessionResult struct {
	PerOrbit         float64 // Precession per orbit (radians)
	PerOrbitArcSec   float64 // Precession per orbit (arcseconds)
	PerCentury       float64 // Precession per century (arcseconds)
	OrbitsPerCentury float64 // Number of orbits per century
}

// Phase 2 Data Structures for N-body Simulation

// Vector3D represents a 3D vector with basic operations.
type Vector3D struct {
	X, Y, Z float64
}

// Add returns the sum of two vectors.
func (v Vector3D) Add(other Vector3D) Vector3D {
	return Vector3D{v.X + other.X, v.Y + other.Y, v.Z + other.Z}
}

// Subtract returns the difference of two vectors.
func (v Vector3D) Subtract(other Vector3D) Vector3D {
	return Vector3D{v.X - other.X, v.Y - other.Y, v.Z - other.Z}
}

// Scale returns the vector scaled by a factor.
func (v Vector3D) Scale(factor float64) Vector3D {
	return Vector3D{v.X * factor, v.Y * factor, v.Z * factor}
}

// Magnitude returns the magnitude of the vector.
func (v Vector3D) Magnitude() float64 {
	return math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
}

// Unit returns a unit vector in the same direction.
func (v Vector3D) Unit() Vector3D {
	mag := v.Magnitude()
	if mag == 0 {
		return Vector3D{0, 0, 0}
	}
	return v.Scale(1.0 / mag)
}

// Dot returns the dot product of two vectors.
func (v Vector3D) Dot(other Vector3D) float64 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}

// Cross returns the cross product of two vectors.
func (v Vector3D) Cross(other Vector3D) Vector3D {
	return Vector3D{
		X: v.Y*other.Z - v.Z*other.Y,
		Y: v.Z*other.X - v.X*other.Z,
		Z: v.X*other.Y - v.Y*other.X,
	}
}

// Orbital Mechanics Functions

// SolveKeplerEquation solves Kepler's equation M = E - e*sin(E) for eccentric anomaly E.
// Uses Newton-Raphson iteration for improved accuracy.
func SolveKeplerEquation(meanAnomaly, eccentricity float64) float64 {
	// Initial guess for eccentric anomaly
	E := meanAnomaly
	if eccentricity > 0.8 {
		E = math.Pi
	}
	
	// Newton-Raphson iteration
	for i := 0; i < 30; i++ {
		f := E - eccentricity*math.Sin(E) - meanAnomaly
		df := 1.0 - eccentricity*math.Cos(E)
		
		deltaE := f / df
		E -= deltaE
		
		// Check for convergence
		if math.Abs(deltaE) < 1e-15 {
			break
		}
	}
	
	return E
}

// KeplerToCartesian converts Keplerian orbital elements to Cartesian state vectors.
// Returns position and velocity vectors in the reference frame.
func KeplerToCartesian(elements OrbitalElements, centralMass float64) (Vector3D, Vector3D) {
	a := elements.SemiMajorAxis
	e := elements.Eccentricity
	i := elements.Inclination
	Omega := elements.LongitudeOfNode
	omega := elements.ArgumentOfPeriapsis
	M := elements.MeanAnomaly
	
	// Gravitational parameter
	mu := GravitationalConstant * centralMass
	
	// Solve Kepler's equation for eccentric anomaly
	E := SolveKeplerEquation(M, e)
	
	// True anomaly
	nu := 2.0 * math.Atan2(math.Sqrt(1+e)*math.Sin(E/2), math.Sqrt(1-e)*math.Cos(E/2))
	
	// Distance
	r := a * (1 - e*math.Cos(E))
	
	// Position and velocity in orbital plane
	xOrb := r * math.Cos(nu)
	yOrb := r * math.Sin(nu)
	
	// Specific angular momentum
	h := math.Sqrt(mu * a * (1 - e*e))
	
	// Velocity in orbital plane
	vxOrb := -mu/h * math.Sin(nu)
	vyOrb := mu/h * (e + math.Cos(nu))
	
	// Rotation matrices for 3D transformation
	cosOmega := math.Cos(Omega)
	sinOmega := math.Sin(Omega)
	cosomega := math.Cos(omega)
	sinomega := math.Sin(omega)
	cosi := math.Cos(i)
	sini := math.Sin(i)
	
	// Transform to 3D inertial frame
	x := xOrb*(cosOmega*cosomega - sinOmega*sinomega*cosi) - yOrb*(cosOmega*sinomega + sinOmega*cosomega*cosi)
	y := xOrb*(sinOmega*cosomega + cosOmega*sinomega*cosi) - yOrb*(sinOmega*sinomega - cosOmega*cosomega*cosi)
	z := xOrb*sinomega*sini + yOrb*cosomega*sini
	
	vx := vxOrb*(cosOmega*cosomega - sinOmega*sinomega*cosi) - vyOrb*(cosOmega*sinomega + sinOmega*cosomega*cosi)
	vy := vxOrb*(sinOmega*cosomega + cosOmega*sinomega*cosi) - vyOrb*(sinOmega*sinomega - cosOmega*cosomega*cosi)
	vz := vxOrb*sinomega*sini + vyOrb*cosomega*sini
	
	return Vector3D{x, y, z}, Vector3D{vx, vy, vz}
}

// OrbitalElements represents Keplerian orbital elements.
type OrbitalElements struct {
	SemiMajorAxis       float64 // a - Semi-major axis (m)
	Eccentricity        float64 // e - Eccentricity
	Inclination         float64 // i - Inclination (radians)
	LongitudeOfNode     float64 // Ω - Longitude of ascending node (radians)
	ArgumentOfPeriapsis float64 // ω - Argument of periapsis (radians)
	MeanAnomaly         float64 // M - Mean anomaly (radians)
}

// Planet represents a celestial body in the N-body simulation.
type Planet struct {
	Name     string          // Planet name
	Mass     float64         // Mass (kg)
	Position Vector3D        // Current position (m)
	Velocity Vector3D        // Current velocity (m/s)
	Elements OrbitalElements // Orbital elements
}

// ForceContribution represents the gravitational contribution from one planet to another.
type ForceContribution struct {
	SourcePlanet string   // Name of the planet causing the force
	Force        Vector3D // Force vector (N)
	Magnitude    float64  // Force magnitude (N)
	Contribution float64  // Contribution to precession (arcseconds/century)
}

// SystemState represents the complete state of the N-body system at a given time.
type SystemState struct {
	Time    float64  // Time (seconds since epoch)
	Planets []Planet // All planet states
}

// Physics Engine Implementation

// CalculatePostNewtonianCorrection computes the 1PN (first post-Newtonian) relativistic correction
// to the gravitational acceleration. This produces the General Relativity precession effect.
func CalculatePostNewtonianCorrection(planet Planet, centralMass float64) Vector3D {
	r := planet.Position
	v := planet.Velocity
	rMag := r.Magnitude()
	vMag := v.Magnitude()
	
	if rMag == 0 {
		return Vector3D{0, 0, 0}
	}
	
	// Gravitational parameter μ = GM
	mu := GravitationalConstant * centralMass
	
	// Unit position vector r̂
	rHat := r.Unit()
	
	// Radial velocity vr = v·r̂
	vr := v.Dot(rHat)
	
	// 1PN acceleration terms:
	// a_GR = [(4μ/r - v²)μ/(c²r²)] * r̂ + [4μvr/(c²r²)] * v
	
	c2 := float64(SpeedOfLight * SpeedOfLight)
	
	// First term: radial correction
	radialTerm := (4*mu/rMag - vMag*vMag) * mu / (c2 * rMag * rMag)
	radialAccel := rHat.Scale(radialTerm)
	
	// Second term: tangential correction  
	tangentialTerm := 4 * mu * vr / (c2 * rMag * rMag)
	tangentialAccel := v.Scale(tangentialTerm)
	
	return radialAccel.Add(tangentialAccel)
}

// CalculateGravitationalForce computes the gravitational force between two planets.
// Returns the force vector acting on planet1 due to planet2.
func CalculateGravitationalForce(planet1, planet2 Planet) Vector3D {
	// Calculate displacement vector from planet1 to planet2
	displacement := planet2.Position.Subtract(planet1.Position)
	distance := displacement.Magnitude()

	// Avoid division by zero for very close objects
	if distance < 1e3 { // 1 km minimum distance
		return Vector3D{0, 0, 0}
	}

	// Newton's law of universal gravitation: F = G * m1 * m2 / r²
	forceMagnitude := GravitationalConstant * planet1.Mass * planet2.Mass / (distance * distance)

	// Force direction: from planet1 toward planet2
	forceDirection := displacement.Unit()

	return forceDirection.Scale(forceMagnitude)
}

// CalculateNetForce computes the total gravitational force on a target planet from all other planets.
func CalculateNetForce(targetIndex int, planets []Planet) Vector3D {
	var netForce Vector3D

	for i, otherPlanet := range planets {
		if i != targetIndex {
			force := CalculateGravitationalForce(planets[targetIndex], otherPlanet)
			netForce = netForce.Add(force)
		}
	}

	return netForce
}

// StateDerivative represents the time derivatives of position and velocity.
type StateDerivative struct {
	PositionDerivative Vector3D // dp/dt = velocity
	VelocityDerivative Vector3D // dv/dt = acceleration
}

// ComputeDerivative calculates the derivative of a planet's state (velocity and acceleration).
func ComputeDerivative(planetIndex int, planets []Planet) StateDerivative {
	planet := planets[planetIndex]
	netForce := CalculateNetForce(planetIndex, planets)
	acceleration := netForce.Scale(1.0 / planet.Mass) // F = ma, so a = F/m
	
	// Add post-Newtonian (1PN) relativistic corrections for planets orbiting the Sun
	// Assume Sun is at index 0 and is stationary at origin
	if planetIndex > 0 && len(planets) > 0 {
		// Calculate 1PN correction relative to the Sun
		pnCorrection := CalculatePostNewtonianCorrection(planet, planets[0].Mass)
		acceleration = acceleration.Add(pnCorrection)
	}

	return StateDerivative{
		PositionDerivative: planet.Velocity,
		VelocityDerivative: acceleration,
	}
}

// IntegrateRK4 performs one step of Runge-Kutta 4th order integration.
// Updates the system state by time step dt.
func IntegrateRK4(state SystemState, dt float64) SystemState {
	n := len(state.Planets)

	// Storage for RK4 intermediate values
	k1 := make([]StateDerivative, n)
	k2 := make([]StateDerivative, n)
	k3 := make([]StateDerivative, n)
	k4 := make([]StateDerivative, n)

	// k1: derivatives at current state
	for i := range state.Planets {
		k1[i] = ComputeDerivative(i, state.Planets)
	}

	// k2: derivatives at state + k1*dt/2
	tempPlanets := make([]Planet, n)
	for i := range state.Planets {
		tempPlanets[i] = state.Planets[i]
		tempPlanets[i].Position = state.Planets[i].Position.Add(k1[i].PositionDerivative.Scale(dt / 2))
		tempPlanets[i].Velocity = state.Planets[i].Velocity.Add(k1[i].VelocityDerivative.Scale(dt / 2))
	}
	for i := range tempPlanets {
		k2[i] = ComputeDerivative(i, tempPlanets)
	}

	// k3: derivatives at state + k2*dt/2
	for i := range state.Planets {
		tempPlanets[i] = state.Planets[i]
		tempPlanets[i].Position = state.Planets[i].Position.Add(k2[i].PositionDerivative.Scale(dt / 2))
		tempPlanets[i].Velocity = state.Planets[i].Velocity.Add(k2[i].VelocityDerivative.Scale(dt / 2))
	}
	for i := range tempPlanets {
		k3[i] = ComputeDerivative(i, tempPlanets)
	}

	// k4: derivatives at state + k3*dt
	for i := range state.Planets {
		tempPlanets[i] = state.Planets[i]
		tempPlanets[i].Position = state.Planets[i].Position.Add(k3[i].PositionDerivative.Scale(dt))
		tempPlanets[i].Velocity = state.Planets[i].Velocity.Add(k3[i].VelocityDerivative.Scale(dt))
	}
	for i := range tempPlanets {
		k4[i] = ComputeDerivative(i, tempPlanets)
	}

	// Combine all derivatives: newState = oldState + (k1 + 2*k2 + 2*k3 + k4) * dt/6
	newState := SystemState{
		Time:    state.Time + dt,
		Planets: make([]Planet, n),
	}

	for i := range state.Planets {
		// Combined derivative
		positionChange := k1[i].PositionDerivative.Add(k2[i].PositionDerivative.Scale(2)).Add(k3[i].PositionDerivative.Scale(2)).Add(k4[i].PositionDerivative).Scale(dt / 6)
		velocityChange := k1[i].VelocityDerivative.Add(k2[i].VelocityDerivative.Scale(2)).Add(k3[i].VelocityDerivative.Scale(2)).Add(k4[i].VelocityDerivative).Scale(dt / 6)

		// Update planet state
		newState.Planets[i] = state.Planets[i]
		newState.Planets[i].Position = state.Planets[i].Position.Add(positionChange)
		newState.Planets[i].Velocity = state.Planets[i].Velocity.Add(velocityChange)
	}

	return newState
}

// IntegrateVelocityVerlet performs one step of Velocity Verlet (symplectic) integration.
// This integrator conserves energy better for long-term orbital simulations.
func IntegrateVelocityVerlet(state SystemState, dt float64) SystemState {
	n := len(state.Planets)
	
	// Create new state
	newState := SystemState{
		Time:    state.Time + dt,
		Planets: make([]Planet, n),
	}
	
	// Calculate current accelerations
	currentAccelerations := make([]Vector3D, n)
	for i := range state.Planets {
		netForce := CalculateNetForce(i, state.Planets)
		currentAccelerations[i] = netForce.Scale(1.0 / state.Planets[i].Mass)
	}
	
	// Velocity Verlet algorithm:
	// 1. Update positions: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt²
	// 2. Calculate new accelerations at r(t+dt)
	// 3. Update velocities: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
	
	// Step 1: Update positions
	for i := range state.Planets {
		newState.Planets[i] = state.Planets[i]
		newState.Planets[i].Position = state.Planets[i].Position.
			Add(state.Planets[i].Velocity.Scale(dt)).
			Add(currentAccelerations[i].Scale(0.5 * dt * dt))
	}
	
	// Step 2: Calculate new accelerations at updated positions
	newAccelerations := make([]Vector3D, n)
	for i := range newState.Planets {
		netForce := CalculateNetForce(i, newState.Planets)
		newAccelerations[i] = netForce.Scale(1.0 / newState.Planets[i].Mass)
	}
	
	// Step 3: Update velocities using average of old and new accelerations
	for i := range state.Planets {
		averageAcceleration := currentAccelerations[i].Add(newAccelerations[i]).Scale(0.5)
		newState.Planets[i].Velocity = state.Planets[i].Velocity.Add(averageAcceleration.Scale(dt))
	}
	
	return newState
}

// Precession Analysis Functions

// ConvertStateToOrbitalElements converts Cartesian state (position, velocity) to Keplerian orbital elements.
// Uses proper eccentricity vector calculation for accurate orbital parameter determination.
func ConvertStateToOrbitalElements(planet Planet, centralMass float64) OrbitalElements {
	r := planet.Position
	v := planet.Velocity
	rMag := r.Magnitude()
	vMag := v.Magnitude()
	
	// Gravitational parameter
	mu := GravitationalConstant * centralMass
	
	// Specific orbital energy
	specificEnergy := 0.5*vMag*vMag - mu/rMag
	
	// Semi-major axis from energy: E = -μ/(2a)
	semiMajorAxis := -mu / (2 * specificEnergy)
	
	// Angular momentum vector h = r × v
	hVec := r.Cross(v)
	h := hVec.Magnitude()
	
	// Eccentricity vector (Laplace-Runge-Lenz vector)
	// e_vec = (v × h)/μ - r/|r|
	vCrossH := v.Cross(hVec)
	rUnit := r.Unit()
	eVec := vCrossH.Scale(1.0/mu).Subtract(rUnit)
	eccentricity := eVec.Magnitude()
	
	// Node vector (for inclination calculations)
	// n = k × h (where k is the z-axis unit vector)
	kVec := Vector3D{0, 0, 1}
	nVec := kVec.Cross(hVec)
	nMag := nVec.Magnitude()
	
	// Inclination: i = arccos(h_z / |h|)
	inclination := math.Acos(hVec.Z / h)
	
	// Longitude of ascending node: Ω = arctan2(n_y, n_x)
	var longitudeOfNode float64
	if nMag > 1e-10 {
		longitudeOfNode = math.Atan2(nVec.Y, nVec.X)
	} else {
		longitudeOfNode = 0.0 // Undefined for non-inclined orbits
	}
	
	// Argument of periapsis: ω = arccos(n·e / (|n||e|))
	var argumentOfPeriapsis float64
	if eccentricity > 1e-10 && nMag > 1e-10 {
		cosω := nVec.Dot(eVec) / (nMag * eccentricity)
		// Clamp to valid range to avoid numerical errors
		if cosω > 1.0 {
			cosω = 1.0
		} else if cosω < -1.0 {
			cosω = -1.0
		}
		argumentOfPeriapsis = math.Acos(cosω)
		
		// Check quadrant
		if eVec.Z < 0 {
			argumentOfPeriapsis = 2*math.Pi - argumentOfPeriapsis
		}
	} else if eccentricity > 1e-10 {
		// For non-inclined orbits, measure from x-axis
		argumentOfPeriapsis = math.Atan2(eVec.Y, eVec.X)
	} else {
		argumentOfPeriapsis = 0.0 // Undefined for circular orbits
	}
	
	// True anomaly: ν = arccos(e·r / (|e||r|))
	var trueAnomaly float64
	if eccentricity > 1e-10 {
		cosν := eVec.Dot(r) / (eccentricity * rMag)
		// Clamp to valid range
		if cosν > 1.0 {
			cosν = 1.0
		} else if cosν < -1.0 {
			cosν = -1.0
		}
		trueAnomaly = math.Acos(cosν)
		
		// Check quadrant using velocity direction
		if r.Dot(v) < 0 {
			trueAnomaly = 2*math.Pi - trueAnomaly
		}
	} else {
		// For circular orbits, use position angle
		trueAnomaly = math.Atan2(r.Y, r.X)
	}
	
	// Convert true anomaly to mean anomaly (simplified for demonstration)
	// For more accuracy, would need to convert through eccentric anomaly
	meanAnomaly := trueAnomaly // Approximation for low eccentricity
	
	// Ensure angles are in [0, 2π] range
	if longitudeOfNode < 0 {
		longitudeOfNode += 2 * math.Pi
	}
	if argumentOfPeriapsis < 0 {
		argumentOfPeriapsis += 2 * math.Pi
	}
	if meanAnomaly < 0 {
		meanAnomaly += 2 * math.Pi
	}
	
	return OrbitalElements{
		SemiMajorAxis:       semiMajorAxis,
		Eccentricity:        eccentricity,
		Inclination:         inclination,
		LongitudeOfNode:     longitudeOfNode,
		ArgumentOfPeriapsis: argumentOfPeriapsis,
		MeanAnomaly:         meanAnomaly,
	}
}

// CalculatePerihelionLongitude calculates the longitude of perihelion.
// This is the sum of longitude of ascending node and argument of periapsis.
func CalculatePerihelionLongitude(elements OrbitalElements) float64 {
	return elements.LongitudeOfNode + elements.ArgumentOfPeriapsis
}

// PeriapsisEvent represents a detected periapsis passage
type PeriapsisEvent struct {
	Time              float64 // Time of periapsis passage (seconds)
	Distance          float64 // Distance to Sun at periapsis (m)
	ArgumentPeriapsis float64 // Argument of periapsis (radians)
	OrbitNumber       int     // Sequential orbit number
}

// PrecessionTracker tracks perihelion precession over time using periapsis detection.
type PrecessionTracker struct {
	PeriapsisEvents   []PeriapsisEvent // History of detected periapsis passages
	LastPosition      Vector3D         // Previous position for periapsis detection
	LastVelocity      Vector3D         // Previous velocity for periapsis detection
	LastTime          float64          // Previous time for interpolation
	LastRadialVel     float64          // Previous radial velocity (r·v/|r|)
	CurrentOrbit      int              // Current orbit number
	TotalPrecession   float64          // Total accumulated precession (radians)
	PrecessionHistory []float64        // History of perihelion longitudes for visualization
	TimeHistory       []float64        // Corresponding time points for visualization
	lastOmega         float64          // ω on previous pericentre (unwrapped)
	cumOmega          float64          // accumulated "unwrapped" angle
}

// NewPrecessionTracker creates a new precession tracker.
func NewPrecessionTracker(initialLongitude float64) *PrecessionTracker {
	return &PrecessionTracker{
		PeriapsisEvents:   []PeriapsisEvent{},
		LastPosition:      Vector3D{},
		LastVelocity:      Vector3D{},
		LastTime:          0.0,
		LastRadialVel:     0.0,
		CurrentOrbit:      0,
		TotalPrecession:   0.0,
		PrecessionHistory: []float64{initialLongitude},
		TimeHistory:       []float64{0.0},
		lastOmega:         0.0,
		cumOmega:          0.0,
	}
}

// Update detects periapsis passages and tracks precession using proper r·v sign change detection.
func (pt *PrecessionTracker) Update(time float64, position Vector3D, velocity Vector3D) {
	// Calculate radial velocity: r·v / |r|
	rMag := position.Magnitude()
	if rMag == 0 {
		return // Avoid division by zero
	}
	
	currentRadialVel := position.Dot(velocity) / rMag
	
	// Detect periapsis passage using radial velocity sign change
	// Periapsis occurs when radial velocity changes from negative (approaching) to positive (receding)
	if pt.LastTime > 0 && pt.LastRadialVel < 0 && currentRadialVel > 0 {
		// Interpolate to find exact periapsis time and position
		periapsisTime, periapsisPos, periapsisVel := pt.interpolatePeriapsis(
			pt.LastTime, time, pt.LastPosition, position, pt.LastVelocity, velocity,
			pt.LastRadialVel, currentRadialVel)
		
		// Record periapsis event
		pt.detectPeriapsis(periapsisTime, periapsisPos, periapsisVel)
	}
	
	// Update tracking variables
	pt.LastPosition = position
	pt.LastVelocity = velocity
	pt.LastTime = time
	pt.LastRadialVel = currentRadialVel
	
	// Update visualization data (for continuity with existing plots)
	elements := ConvertStateToOrbitalElements(Planet{Position: position, Velocity: velocity}, SolarMass)
	longitude := CalculatePerihelionLongitude(elements)
	pt.PrecessionHistory = append(pt.PrecessionHistory, longitude*RadiansToArcseconds)
	pt.TimeHistory = append(pt.TimeHistory, time/(365.25*24*3600))
}

// interpolatePeriapsis performs linear interpolation to find the exact periapsis moment.
// Uses the fact that radial velocity is zero at periapsis.
func (pt *PrecessionTracker) interpolatePeriapsis(t1, t2 float64, r1, r2, v1, v2 Vector3D, rv1, rv2 float64) (float64, Vector3D, Vector3D) {
	// Linear interpolation to find when radial velocity = 0
	// rv(t) = rv1 + (rv2 - rv1) * (t - t1) / (t2 - t1) = 0
	// Solving: t = t1 - rv1 * (t2 - t1) / (rv2 - rv1)
	
	var alpha float64
	if math.Abs(rv2 - rv1) < 1e-15 {
		// Avoid division by zero, use midpoint
		alpha = 0.5
	} else {
		alpha = -rv1 / (rv2 - rv1)
	}
	
	// Clamp alpha to [0, 1] to stay within the interval
	if alpha < 0 {
		alpha = 0
	} else if alpha > 1 {
		alpha = 1
	}
	
	// Interpolate time, position, and velocity
	periapsisTime := t1 + alpha*(t2-t1)
	periapsisPos := r1.Add(r2.Subtract(r1).Scale(alpha))
	periapsisVel := v1.Add(v2.Subtract(v1).Scale(alpha))
	
	return periapsisTime, periapsisPos, periapsisVel
}

// detectPeriapsis records a periapsis passage and calculates precession.
func (pt *PrecessionTracker) detectPeriapsis(time float64, position Vector3D, velocity Vector3D) {
	// Calculate distance at periapsis
	distance := position.Magnitude()
	
	// Use orbital elements to get proper argument of periapsis
	elements := ConvertStateToOrbitalElements(Planet{Position: position, Velocity: velocity}, SolarMass)
	argumentPeriapsis := elements.ArgumentOfPeriapsis
	
	// Create periapsis event
	event := PeriapsisEvent{
		Time:              time,
		Distance:          distance,
		ArgumentPeriapsis: argumentPeriapsis,
		OrbitNumber:       pt.CurrentOrbit,
	}
	
	// Add to history
	pt.PeriapsisEvents = append(pt.PeriapsisEvents, event)
	
	// Debug output for first few periapsis events (optional)
	if len(pt.PeriapsisEvents) <= 10 {
		fmt.Printf("Detected periapsis %d at day %.1f, ω=%.6f rad (%.3f°), cumOmega=%.6f rad\n",
			pt.CurrentOrbit, time/(24*3600), argumentPeriapsis, argumentPeriapsis*180/math.Pi, pt.cumOmega)
	}
	
	// Update precession using cumulative tracking approach
	pt.addPeriapsis(argumentPeriapsis)
}

// addPeriapsis updates precession tracking for each new periapsis event
func (pt *PrecessionTracker) addPeriapsis(eventOmega float64) {
	if pt.CurrentOrbit == 0 { // First periapsis
		pt.lastOmega = eventOmega
	} else {
		dOmega := eventOmega - pt.lastOmega
		originalDOmega := dOmega // Store original for debug
		
		// Unwrap angle change at transition, not over entire interval
		if dOmega > math.Pi {
			dOmega -= 2 * math.Pi
		}
		if dOmega < -math.Pi {
			dOmega += 2 * math.Pi
		}

		// Debug output for first few transitions
		if pt.CurrentOrbit <= 10 {
			fmt.Printf("  Orbit %d: dOmega=%.6f rad (original=%.6f), cumOmega before=%.6f\n", 
				pt.CurrentOrbit, dOmega, originalDOmega, pt.cumOmega)
		}

		pt.cumOmega += dOmega                // Sum the increment
		pt.TotalPrecession = pt.cumOmega     // Store in radians
		pt.lastOmega = eventOmega
	}
	pt.CurrentOrbit++
}


// GetPrecessionRate calculates the precession rate in arcseconds per year.
func (pt *PrecessionTracker) GetPrecessionRate(timeSpanSeconds float64) float64 {
	if timeSpanSeconds == 0 || len(pt.PeriapsisEvents) < 2 {
		return 0
	}

	// TotalPrecession now contains cumulative precession in radians over the simulation
	// Convert to arcseconds per year
	secondsPerYear := 365.25 * 24 * 3600
	durationYears := timeSpanSeconds / secondsPerYear
	
	// Total precession in arcseconds over the simulation duration
	totalPrecessionArcsec := pt.TotalPrecession * RadiansToArcseconds
	
	// Return as arcseconds per year
	return totalPrecessionArcsec / durationYears
}

// GetOrbitalCount returns the number of completed orbits based on periapsis detections.
func (pt *PrecessionTracker) GetOrbitalCount() int {
	return pt.CurrentOrbit
}

// PrecessionAnalyzer decomposes precession contributions from different sources.
type PrecessionAnalyzer struct {
	MercuryIndex   int                           // Index of Mercury in the planet array
	Tracker        *PrecessionTracker            // Main precession tracker
	PlanetTrackers map[string]*PrecessionTracker // Individual planet contribution trackers
	CentralMass    float64                       // Mass of the central body (Sun)
}

// NewPrecessionAnalyzer creates a new precession analyzer.
func NewPrecessionAnalyzer(mercuryIndex int, centralMass float64, initialLongitude float64) *PrecessionAnalyzer {
	return &PrecessionAnalyzer{
		MercuryIndex:   mercuryIndex,
		Tracker:        NewPrecessionTracker(initialLongitude),
		PlanetTrackers: make(map[string]*PrecessionTracker),
		CentralMass:    centralMass,
	}
}

// AnalyzeSystem performs precession analysis on the current system state.
func (pa *PrecessionAnalyzer) AnalyzeSystem(state SystemState) {
	if pa.MercuryIndex >= len(state.Planets) {
		return
	}

	mercury := state.Planets[pa.MercuryIndex]

	// Update main tracker with position and velocity for periapsis detection
	pa.Tracker.Update(state.Time, mercury.Position, mercury.Velocity)
}

// GetContributionsByPlanet calculates the precession contribution from each planet.
// This is a simplified version - a full implementation would run separate simulations.
func (pa *PrecessionAnalyzer) GetContributionsByPlanet(planets []Planet) map[string]float64 {
	contributions := make(map[string]float64)

	// For now, return the known empirical values
	// In a full implementation, this would run separate N-body simulations
	// with individual planets to isolate their contributions
	contributions["Venus"] = 277.8
	contributions["Earth"] = 90.0
	contributions["Jupiter"] = 153.6
	contributions["Mars"] = 2.5
	contributions["Saturn"] = 0.8
	contributions["Others"] = 6.7

	return contributions
}

// Real Planetary Data and N-body Simulation Implementation

// PlanetaryData holds real astronomical data for solar system planets
type PlanetaryData struct {
	Name                string
	Mass                float64 // kg
	SemiMajorAxis       float64 // m
	Eccentricity        float64
	OrbitalPeriod       float64 // days
	MeanOrbitalVelocity float64 // m/s
}

// GetPlanetaryData returns real astronomical data for solar system planets
func GetPlanetaryData() map[string]PlanetaryData {
	return map[string]PlanetaryData{
		"Sun": {
			Name: "Sun",
			Mass: SolarMass,
		},
		"Mercury": {
			Name:                "Mercury",
			Mass:                3.30104e23,           // kg
			SemiMajorAxis:       MercurySemiMajorAxis, // m
			Eccentricity:        MercuryEccentricity,  //
			OrbitalPeriod:       MercuryOrbitalPeriod, // days
			MeanOrbitalVelocity: 47362,                // m/s
		},
		"Venus": {
			Name:                "Venus",
			Mass:                4.86732e24,                  // kg
			SemiMajorAxis:       0.723332 * AstronomicalUnit, // m
			Eccentricity:        0.006772,                    //
			OrbitalPeriod:       224.701,                     // days
			MeanOrbitalVelocity: 35020,                       // m/s
		},
		"Earth": {
			Name:                "Earth",
			Mass:                5.97219e24,             // kg
			SemiMajorAxis:       1.0 * AstronomicalUnit, // m
			Eccentricity:        0.0167086,              //
			OrbitalPeriod:       365.256,                // days
			MeanOrbitalVelocity: 29780,                  // m/s
		},
		"Mars": {
			Name:                "Mars",
			Mass:                6.4169e23,                   // kg
			SemiMajorAxis:       1.523679 * AstronomicalUnit, // m
			Eccentricity:        0.0934,                      //
			OrbitalPeriod:       686.980,                     // days
			MeanOrbitalVelocity: 24070,                       // m/s
		},
		"Jupiter": {
			Name:                "Jupiter",
			Mass:                1.89813e27,                  // kg
			SemiMajorAxis:       5.204267 * AstronomicalUnit, // m
			Eccentricity:        0.048775,                    //
			OrbitalPeriod:       4332.59,                     // days
			MeanOrbitalVelocity: 13070,                       // m/s
		},
		"Saturn": {
			Name:                "Saturn",
			Mass:                5.68319e26,                  // kg
			SemiMajorAxis:       9.573635 * AstronomicalUnit, // m
			Eccentricity:        0.055723,                    //
			OrbitalPeriod:       10759.22,                    // days
			MeanOrbitalVelocity: 9680,                        // m/s
		},
	}
}

// InitializePlanet creates a Planet with initial position and velocity from orbital elements
// Using proper Kepler-to-Cartesian conversion for accurate initial conditions.
func InitializePlanet(data PlanetaryData, meanAnomaly float64) Planet {
	if data.Name == "Sun" {
		return Planet{
			Name:     data.Name,
			Mass:     data.Mass,
			Position: Vector3D{0, 0, 0}, // Sun at origin
			Velocity: Vector3D{0, 0, 0}, // Sun stationary
			Elements: OrbitalElements{},
		}
	}

	// Create orbital elements structure
	elements := OrbitalElements{
		SemiMajorAxis:       data.SemiMajorAxis,
		Eccentricity:        data.Eccentricity,
		Inclination:         0, // Simplified coplanar for now
		LongitudeOfNode:     0,
		ArgumentOfPeriapsis: 0,
		MeanAnomaly:         meanAnomaly,
	}

	// Use proper Kepler-to-Cartesian conversion
	position, velocity := KeplerToCartesian(elements, SolarMass)

	return Planet{
		Name:     data.Name,
		Mass:     data.Mass,
		Position: position,
		Velocity: velocity,
		Elements: elements,
	}
}

// CreateSolarSystem initializes all planets in the solar system
func CreateSolarSystem() SystemState {
	planetData := GetPlanetaryData()

	planets := []Planet{
		InitializePlanet(planetData["Sun"], 0),
		InitializePlanet(planetData["Mercury"], 0),          // Mercury at perihelion
		InitializePlanet(planetData["Venus"], math.Pi/4),    // Venus at 45°
		InitializePlanet(planetData["Earth"], math.Pi/2),    // Earth at 90°
		InitializePlanet(planetData["Mars"], 3*math.Pi/4),   // Mars at 135°
		InitializePlanet(planetData["Jupiter"], math.Pi),    // Jupiter at 180°
		InitializePlanet(planetData["Saturn"], 5*math.Pi/4), // Saturn at 225°
	}

	return SystemState{
		Time:    0,
		Planets: planets,
	}
}

// RunNBodySimulation executes the N-body simulation with precession tracking
func RunNBodySimulation(durationYears float64, outputInterval int) {
	fmt.Println("\n=== N-body Simulation Starting ===")

	// Initialize system
	state := CreateSolarSystem()
	mercuryIndex := 1 // Mercury is at index 1

	// Setup precession tracking
	mercury := state.Planets[mercuryIndex]
	initialElements := ConvertStateToOrbitalElements(mercury, SolarMass)
	initialLongitude := CalculatePerihelionLongitude(initialElements)
	analyzer := NewPrecessionAnalyzer(mercuryIndex, SolarMass, initialLongitude)

	// Simulation parameters
	durationSeconds := durationYears * 365.25 * 24 * 3600
	timeStep := 3600.0 * 1 // 1 hour in seconds (improved precision)
	totalSteps := int(durationSeconds / timeStep)

	fmt.Printf("Duration: %.1f years (%.0f seconds)\n", durationYears, durationSeconds)
	fmt.Printf("Time step: %.1f hours\n", timeStep/3600)
	fmt.Printf("Total steps: %d\n", totalSteps)
	fmt.Printf("Output every %d steps\n", outputInterval)

	// Storage for visualization data
	var orbitPoints []Vector3D
	var precessionHistory []float64
	var timeHistory []float64

	// Run simulation
	for step := 0; step < totalSteps; step++ {
		// Integrate one step using symplectic integrator for better long-term stability
		state = IntegrateVelocityVerlet(state, timeStep)

		// Track precession
		analyzer.AnalyzeSystem(state)

		// Store data for visualization
		if step%outputInterval == 0 {
			mercury := state.Planets[mercuryIndex]
			orbitPoints = append(orbitPoints, mercury.Position)

			elements := ConvertStateToOrbitalElements(mercury, SolarMass)
			longitude := CalculatePerihelionLongitude(elements)
			precessionHistory = append(precessionHistory, longitude*RadiansToArcseconds)
			timeHistory = append(timeHistory, state.Time/(365.25*24*3600)) // years

			// Progress update with diagnostics
			if step%(totalSteps/10) == 0 {
				progress := float64(step) / float64(totalSteps) * 100
				currentYear := state.Time / (365.25 * 24 * 3600)
				fmt.Printf("Progress: %.1f%% (%.2f years) - Longitude: %.3f° Orbits: %d Precession: %.6f\"\n",
					progress, currentYear, longitude*180/math.Pi, analyzer.Tracker.GetOrbitalCount(), analyzer.Tracker.TotalPrecession*RadiansToArcseconds)
			}
		}
	}

	// Analysis results
	fmt.Println("\n=== N-body Simulation Results ===")
	finalPrecessionRate := analyzer.Tracker.GetPrecessionRate(durationSeconds)
	precessionPerCentury := finalPrecessionRate * 100.0

	fmt.Printf("Total simulation time: %.1f years\n", durationYears)
	fmt.Printf("Mercury orbits completed: %d\n", analyzer.Tracker.GetOrbitalCount())
	fmt.Printf("Periapsis events detected: %d\n", len(analyzer.Tracker.PeriapsisEvents))

	// Calculate theoretical values
	grResult := CalculateGRPrecession(SolarMass, MercurySemiMajorAxis, MercuryEccentricity, MercuryOrbitalPeriod)
	planetary := CalculatePlanetaryPerturbations()

	fmt.Println("\n=== Precession Analysis ===")
	fmt.Printf("N-body simulation detected: %.2f arcsec/century\n", precessionPerCentury)
	fmt.Printf("  (Note: Short simulations are noisy; longer runs give better results)\n")
	fmt.Println()
	fmt.Printf("Theoretical contributions:\n")
	fmt.Printf("  Classical perturbations (planets): %.1f arcsec/century\n", planetary.Total())
	fmt.Printf("  General Relativity effect: %.1f arcsec/century\n", grResult.PerCentury)
	fmt.Printf("  Total theoretical: %.1f arcsec/century\n", planetary.Total()+grResult.PerCentury)
	fmt.Printf("  Observed (historical): %.1f arcsec/century\n", ObservedPrecession)
	fmt.Println()
	fmt.Printf("Agreement with theory: %.1f arcsec difference\n",
		math.Abs((planetary.Total()+grResult.PerCentury)-ObservedPrecession))

	// Generate visualizations
	GenerateOrbitVisualization(orbitPoints, "mercury_orbit.svg")
	GeneratePrecessionPlot(timeHistory, precessionHistory, "precession_plot.svg")
	GenerateRosettePattern(orbitPoints, "rosette_pattern.svg")

	fmt.Println("\nVisualization files generated:")
	fmt.Println("- mercury_orbit.svg: Mercury orbital path")
	fmt.Println("- precession_plot.svg: Precession vs time")
	fmt.Println("- rosette_pattern.svg: Orbital rosette pattern")
}

// Visualization Functions

// GenerateOrbitVisualization creates an SVG visualization of Mercury's orbital path
func GenerateOrbitVisualization(orbitPoints []Vector3D, filename string) {
	if len(orbitPoints) == 0 {
		return
	}

	// Find bounds for scaling
	var minX, maxX, minY, maxY float64
	for i, point := range orbitPoints {
		if i == 0 {
			minX, maxX = point.X, point.X
			minY, maxY = point.Y, point.Y
		} else {
			if point.X < minX {
				minX = point.X
			}
			if point.X > maxX {
				maxX = point.X
			}
			if point.Y < minY {
				minY = point.Y
			}
			if point.Y > maxY {
				maxY = point.Y
			}
		}
	}

	// Add padding
	padding := 0.1 * math.Max(maxX-minX, maxY-minY)
	minX -= padding
	maxX += padding
	minY -= padding
	maxY += padding

	width, height := 800.0, 800.0
	scaleX := width / (maxX - minX)
	scaleY := height / (maxY - minY)

	var svg strings.Builder
	svg.WriteString(fmt.Sprintf(`<?xml version="1.0" encoding="UTF-8"?>
<svg width="%.0f" height="%.0f" xmlns="http://www.w3.org/2000/svg">
<style>
.orbit-path { fill: none; stroke: #ff6b35; stroke-width: 1.5; opacity: 0.8; }
.sun { fill: #ffdd44; stroke: #ff8800; stroke-width: 2; }
.mercury { fill: #8c7853; stroke: #654321; stroke-width: 1; }
.grid { stroke: #ddd; stroke-width: 0.5; }
.text { font-family: Arial; font-size: 12px; fill: #333; }
</style>
`, width, height))

	// Grid lines
	for i := 0; i <= 10; i++ {
		x := float64(i) * width / 10
		y := float64(i) * height / 10
		svg.WriteString(fmt.Sprintf(`<line x1="%.1f" y1="0" x2="%.1f" y2="%.0f" class="grid"/>`, x, x, height))
		svg.WriteString(fmt.Sprintf(`<line x1="0" y1="%.1f" x2="%.0f" y2="%.1f" class="grid"/>`, y, width, y))
	}

	// Convert coordinates to SVG space
	convertX := func(x float64) float64 { return (x - minX) * scaleX }
	convertY := func(y float64) float64 { return height - (y-minY)*scaleY } // Flip Y axis

	// Draw orbital path
	svg.WriteString(`<path d="M`)
	for i, point := range orbitPoints {
		x := convertX(point.X)
		y := convertY(point.Y)
		if i == 0 {
			svg.WriteString(fmt.Sprintf("%.2f,%.2f", x, y))
		} else {
			svg.WriteString(fmt.Sprintf(" L%.2f,%.2f", x, y))
		}
	}
	svg.WriteString(`" class="orbit-path"/>`)

	// Draw Sun at center
	sunX := convertX(0)
	sunY := convertY(0)
	svg.WriteString(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="10" class="sun"/>`, sunX, sunY))
	svg.WriteString(fmt.Sprintf(`<text x="%.1f" y="%.1f" class="text">Sun</text>`, sunX+15, sunY+5))

	// Draw Mercury at final position
	if len(orbitPoints) > 0 {
		finalPoint := orbitPoints[len(orbitPoints)-1]
		mercX := convertX(finalPoint.X)
		mercY := convertY(finalPoint.Y)
		svg.WriteString(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="4" class="mercury"/>`, mercX, mercY))
		svg.WriteString(fmt.Sprintf(`<text x="%.1f" y="%.1f" class="text">Mercury</text>`, mercX+8, mercY-8))
	}

	// Title
	svg.WriteString(`<text x="20" y="30" class="text" style="font-size: 16px; font-weight: bold;">Mercury Orbital Path</text>`)
	svg.WriteString(`<text x="20" y="50" class="text">Showing precession effect from N-body simulation</text>`)

	svg.WriteString("</svg>")

	// Write to file
	content := svg.String()
	err := WriteFile(filename, content)
	if err != nil {
		fmt.Printf("Error writing %s: %v\n", filename, err)
	}
}

// GeneratePrecessionPlot creates a plot showing precession vs time
func GeneratePrecessionPlot(timeHistory, precessionHistory []float64, filename string) {
	if len(timeHistory) == 0 || len(precessionHistory) == 0 {
		return
	}

	width, height := 800.0, 600.0
	margin := 60.0

	// Find bounds
	minTime, maxTime := timeHistory[0], timeHistory[len(timeHistory)-1]
	minPrec, maxPrec := precessionHistory[0], precessionHistory[0]
	for _, p := range precessionHistory {
		if p < minPrec {
			minPrec = p
		}
		if p > maxPrec {
			maxPrec = p
		}
	}

	// Add padding to precession range
	precRange := maxPrec - minPrec
	if precRange == 0 {
		precRange = 1
	}
	minPrec -= precRange * 0.1
	maxPrec += precRange * 0.1

	var svg strings.Builder
	svg.WriteString(fmt.Sprintf(`<?xml version="1.0" encoding="UTF-8"?>
<svg width="%.0f" height="%.0f" xmlns="http://www.w3.org/2000/svg">
<style>
.plot-line { fill: none; stroke: #2e86de; stroke-width: 2; }
.axis { stroke: #333; stroke-width: 1; }
.grid { stroke: #eee; stroke-width: 0.5; }
.text { font-family: Arial; font-size: 11px; fill: #333; }
.title { font-family: Arial; font-size: 16px; font-weight: bold; fill: #333; }
</style>
`, width, height))

	plotWidth := width - 2*margin
	plotHeight := height - 2*margin

	// Convert coordinates
	convertX := func(t float64) float64 { return margin + (t-minTime)/(maxTime-minTime)*plotWidth }
	convertY := func(p float64) float64 { return margin + (maxPrec-p)/(maxPrec-minPrec)*plotHeight }

	// Draw axes
	svg.WriteString(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" class="axis"/>`, margin, margin+plotHeight, margin+plotWidth, margin+plotHeight))
	svg.WriteString(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" class="axis"/>`, margin, margin, margin, margin+plotHeight))

	// Draw grid and labels
	for i := 0; i <= 5; i++ {
		// Vertical grid lines (time)
		x := margin + float64(i)*plotWidth/5
		time := minTime + float64(i)*(maxTime-minTime)/5
		svg.WriteString(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" class="grid"/>`, x, margin, x, margin+plotHeight))
		svg.WriteString(fmt.Sprintf(`<text x="%.1f" y="%.1f" class="text" text-anchor="middle">%.1f</text>`, x, margin+plotHeight+20, time))

		// Horizontal grid lines (precession)
		y := margin + float64(i)*plotHeight/5
		prec := maxPrec - float64(i)*(maxPrec-minPrec)/5
		svg.WriteString(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" class="grid"/>`, margin, y, margin+plotWidth, y))
		svg.WriteString(fmt.Sprintf(`<text x="%.1f" y="%.1f" class="text" text-anchor="end">%.1f</text>`, margin-10, y+4, prec))
	}

	// Draw plot line
	svg.WriteString(`<path d="M`)
	for i, t := range timeHistory {
		x := convertX(t)
		y := convertY(precessionHistory[i])
		if i == 0 {
			svg.WriteString(fmt.Sprintf("%.2f,%.2f", x, y))
		} else {
			svg.WriteString(fmt.Sprintf(" L%.2f,%.2f", x, y))
		}
	}
	svg.WriteString(`" class="plot-line"/>`)

	// Labels and title
	svg.WriteString(fmt.Sprintf(`<text x="%.1f" y="%.1f" class="title" text-anchor="middle">Perihelion Precession vs Time</text>`, width/2, 30))
	svg.WriteString(fmt.Sprintf(`<text x="%.1f" y="%.1f" class="text" text-anchor="middle">Time (years)</text>`, width/2, height-10))

	// Rotated Y label
	svg.WriteString(`<g transform="translate(20,` + fmt.Sprintf("%.1f", height/2) + `) rotate(-90)">`)
	svg.WriteString(`<text class="text" text-anchor="middle">Precession (arcseconds)</text>`)
	svg.WriteString(`</g>`)

	svg.WriteString("</svg>")

	// Write to file
	content := svg.String()
	err := WriteFile(filename, content)
	if err != nil {
		fmt.Printf("Error writing %s: %v\n", filename, err)
	}
}

// GenerateRosettePattern creates a visualization showing the orbital rosette pattern
func GenerateRosettePattern(orbitPoints []Vector3D, filename string) {
	if len(orbitPoints) < 10 {
		return
	}

	// Find bounds
	var minX, maxX, minY, maxY float64
	for i, point := range orbitPoints {
		if i == 0 {
			minX, maxX = point.X, point.X
			minY, maxY = point.Y, point.Y
		} else {
			if point.X < minX {
				minX = point.X
			}
			if point.X > maxX {
				maxX = point.X
			}
			if point.Y < minY {
				minY = point.Y
			}
			if point.Y > maxY {
				maxY = point.Y
			}
		}
	}

	// Ensure square aspect ratio
	rangeX := maxX - minX
	rangeY := maxY - minY
	maxRange := math.Max(rangeX, rangeY)
	centerX := (minX + maxX) / 2
	centerY := (minY + maxY) / 2

	minX = centerX - maxRange/2
	maxX = centerX + maxRange/2
	minY = centerY - maxRange/2
	maxY = centerY + maxRange/2

	width, height := 800.0, 800.0
	scaleX := width / (maxX - minX)
	scaleY := height / (maxY - minY)

	var svg strings.Builder
	svg.WriteString(fmt.Sprintf(`<?xml version="1.0" encoding="UTF-8"?>
<svg width="%.0f" height="%.0f" xmlns="http://www.w3.org/2000/svg">
<style>
.rosette { fill: none; stroke: #e74c3c; stroke-width: 0.8; opacity: 0.7; }
.sun { fill: #f39c12; stroke: #d68910; stroke-width: 2; }
.title { font-family: Arial; font-size: 16px; font-weight: bold; fill: #333; }
.text { font-family: Arial; font-size: 12px; fill: #666; }
</style>
`, width, height))

	// Convert coordinates
	convertX := func(x float64) float64 { return (x - minX) * scaleX }
	convertY := func(y float64) float64 { return height - (y-minY)*scaleY }

	// Draw all orbit points as a continuous path (rosette pattern)
	svg.WriteString(`<path d="M`)
	for i, point := range orbitPoints {
		x := convertX(point.X)
		y := convertY(point.Y)
		if i == 0 {
			svg.WriteString(fmt.Sprintf("%.2f,%.2f", x, y))
		} else {
			svg.WriteString(fmt.Sprintf(" L%.2f,%.2f", x, y))
		}
	}
	svg.WriteString(`" class="rosette"/>`)

	// Draw Sun at center
	sunX := convertX(0)
	sunY := convertY(0)
	svg.WriteString(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="8" class="sun"/>`, sunX, sunY))

	// Title and description
	svg.WriteString(`<text x="20" y="30" class="title">Mercury Orbital Rosette Pattern</text>`)
	svg.WriteString(`<text x="20" y="50" class="text">Pattern formed by precessing elliptical orbit over multiple revolutions</text>`)
	svg.WriteString(fmt.Sprintf(`<text x="20" y="%.0f" class="text">Total orbital points: %d</text>`, height-30, len(orbitPoints)))

	svg.WriteString("</svg>")

	// Write to file
	content := svg.String()
	err := WriteFile(filename, content)
	if err != nil {
		fmt.Printf("Error writing %s: %v\n", filename, err)
	}
}

// WriteFile is a simple file writer for our SVG content
func WriteFile(filename, content string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	_, err = file.WriteString(content)
	return err
}

// CalculateGRPrecession computes perihelion precession using Einstein's General Relativity formula.
// The formula is: Δφ = 6πGM / (c²a(1-e²))
func CalculateGRPrecession(mass, semiMajorAxis, eccentricity, orbitalPeriodDays float64) PrecessionResult {
	numerator := 6 * math.Pi * GravitationalConstant * mass
	denominator := SpeedOfLight * SpeedOfLight * semiMajorAxis * (1 - eccentricity*eccentricity)

	precessionPerOrbit := numerator / denominator
	precessionPerOrbitArcSec := precessionPerOrbit * RadiansToArcseconds
	orbitsPerCentury := DaysPerCentury / orbitalPeriodDays
	precessionPerCentury := precessionPerOrbitArcSec * orbitsPerCentury

	return PrecessionResult{
		PerOrbit:         precessionPerOrbit,
		PerOrbitArcSec:   precessionPerOrbitArcSec,
		PerCentury:       precessionPerCentury,
		OrbitsPerCentury: orbitsPerCentury,
	}
}

// PlanetaryContributions holds the known contributions from different planets.
type PlanetaryContributions struct {
	Venus        float64
	Earth        float64
	Jupiter      float64
	OtherPlanets float64
}

// CalculatePlanetaryPerturbations returns the total classical perturbations from other planets.
// These are empirically determined values in arcseconds per century.
func CalculatePlanetaryPerturbations() PlanetaryContributions {
	return PlanetaryContributions{
		Venus:        277.8, // Venus contribution (arcseconds/century)
		Earth:        90.0,  // Earth contribution (arcseconds/century)
		Jupiter:      153.6, // Jupiter contribution (arcseconds/century)
		OtherPlanets: 10.0,  // Mars, Saturn and others (arcseconds/century)
	}
}

// Total returns the sum of all planetary contributions.
func (p PlanetaryContributions) Total() float64 {
	return p.Venus + p.Earth + p.Jupiter + p.OtherPlanets
}

// printMercuryParameters displays Mercury's orbital parameters.
func printMercuryParameters(result PrecessionResult) {
	fmt.Printf("Mercury Orbital Parameters:\n")
	fmt.Printf("- Semi-major axis: %.3f AU (%.3e m)\n", MercurySemiMajorAxis/AstronomicalUnit, MercurySemiMajorAxis)
	fmt.Printf("- Eccentricity: %.5f\n", MercuryEccentricity)
	fmt.Printf("- Orbital period: %.3f days\n", MercuryOrbitalPeriod)
	fmt.Printf("- Orbits per century: %.1f\n", result.OrbitsPerCentury)
	fmt.Println()
}

// printGRResults displays the General Relativity calculation results.
func printGRResults(result PrecessionResult) {
	fmt.Printf("General Relativity Precession:\n")
	fmt.Printf("- Per orbit: %.6e radians\n", result.PerOrbit)
	fmt.Printf("- Per orbit: %.3f arcseconds\n", result.PerOrbitArcSec)
	fmt.Printf("- Per century: %.1f arcseconds\n", result.PerCentury)
	fmt.Println()
}

// printComparison displays comparison with observational data.
func printComparison(grResult PrecessionResult, planetary PlanetaryContributions) {
	totalPredicted := planetary.Total() + grResult.PerCentury
	difference := math.Abs(totalPredicted - ObservedPrecession)

	fmt.Printf("Planetary contributions (classical mechanics): %.1f\"\n", planetary.Total())
	fmt.Println()

	fmt.Printf("Comparison with observations:\n")
	fmt.Printf("- Observed precession: %.1f\" per century\n", ObservedPrecession)
	fmt.Printf("- Classical prediction: %.1f\"\n", planetary.Total())
	fmt.Printf("- Additional GR contribution: %.1f\"\n", grResult.PerCentury)
	fmt.Printf("- Total prediction (classical + GR): %.1f\"\n", totalPredicted)
	fmt.Printf("- Difference from observations: %.1f\"\n", difference)
	fmt.Println()
}

// printValidation checks if the calculation matches the famous ~43 arcseconds.
func printValidation(result PrecessionResult) {
	if math.Abs(result.PerCentury-ExpectedGREffect) < 1.0 {
		fmt.Println("✓ Calculation successful! Obtained the famous ~43 arcseconds!")
	}
}

// printInterestingFacts displays additional insights about the precession effect.
func printInterestingFacts(result PrecessionResult) {
	fmt.Println("\n=== Interesting Facts ===")
	grPercentage := result.PerCentury / ObservedPrecession * 100
	fmt.Printf("- GR effect is only %.1f%% of total precession\n", grPercentage)
	fmt.Printf("- Per Mercury orbit, perihelion shifts by only %.3f arcseconds\n", result.PerOrbitArcSec)
	fmt.Printf("- This is about 1/1000th the angular size of the Moon!\n")
}

func main() {
	fmt.Println("=== Mercury Perihelion Precession Analysis ===")
	fmt.Println("1. Analytical calculation (Phase 1)")
	fmt.Println("2. N-body simulation (Phase 2)")
	fmt.Println("3. Data sources demo (JPL Horizons, VSOP87)")
	fmt.Print("Choose option (1, 2, or 3): ")

	var choice int
	fmt.Scanf("%d", &choice)

	if choice == 3 {
		DemoDataSources()
	} else if choice == 2 {
		// Phase 2: N-body simulation
		fmt.Print("Enter simulation duration in years (recommended: 10-100): ")
		var durationYears float64
		fmt.Scanf("%f", &durationYears)

		// Validate input
		if durationYears <= 0 || durationYears > 1000 {
			fmt.Printf("Invalid duration. Using 10 years instead.\n")
			durationYears = 10
		}

		outputInterval := 100 // Output every 100 time steps
		if durationYears > 50 {
			outputInterval = 500 // Less frequent output for longer simulations
		}

		RunNBodySimulation(durationYears, outputInterval)
	} else {
		// Phase 1: Analytical calculation
		fmt.Println("\n=== Analytical Mercury Perihelion Precession Calculation ===")
		fmt.Println()

		// Calculate GR precession
		result := CalculateGRPrecession(SolarMass, MercurySemiMajorAxis, MercuryEccentricity, MercuryOrbitalPeriod)

		// Calculate planetary contributions
		planetary := CalculatePlanetaryPerturbations()

		// Display results
		printMercuryParameters(result)
		printGRResults(result)
		printComparison(result, planetary)
		printValidation(result)
		printInterestingFacts(result)

		fmt.Println("\n💡 Tip: Run option 2 for N-body simulation with visualization!")
	}
}

// Data Sources Implementation

// JPLHorizonsResponse represents the response from JPL Horizons API
type JPLHorizonsResponse struct {
	Result string `json:"result"`
}

// EphemerisPoint represents a single ephemeris data point
type EphemerisPoint struct {
	Time     time.Time
	Position Vector3D // km
	Velocity Vector3D // km/s
}

// LoadJPLHorizons fetches ephemeris data from NASA JPL Horizons system
// This is a simplified example - real implementation would handle authentication and parsing
func LoadJPLHorizons(objectID string, startTime, endTime time.Time, stepSize string) ([]EphemerisPoint, error) {
	// JPL Horizons API endpoint (simplified)
	baseURL := "https://ssd.jpl.nasa.gov/api/horizons.api"

	// Prepare query parameters
	params := url.Values{}
	params.Set("format", "json")
	params.Set("COMMAND", objectID) // e.g., "199" for Mercury
	params.Set("OBJ_DATA", "YES")
	params.Set("MAKE_EPHEM", "YES")
	params.Set("EPHEM_TYPE", "VECTORS")
	params.Set("CENTER", "500@10") // Solar System Barycenter
	params.Set("START_TIME", startTime.Format("2006-Jan-02"))
	params.Set("STOP_TIME", endTime.Format("2006-Jan-02"))
	params.Set("STEP_SIZE", stepSize)

	// Make HTTP request
	fullURL := baseURL + "?" + params.Encode()

	resp, err := http.Get(fullURL)
	if err != nil {
		return nil, fmt.Errorf("failed to fetch JPL data: %v", err)
	}
	defer resp.Body.Close()

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("failed to read response: %v", err)
	}

	// Parse JSON response
	var horizonsResp JPLHorizonsResponse
	if err := json.Unmarshal(body, &horizonsResp); err != nil {
		return nil, fmt.Errorf("failed to parse JSON: %v", err)
	}

	// Parse the result string to extract ephemeris data
	// This is simplified - real implementation would parse the complex text format
	ephemeris := parseHorizonsResult(horizonsResp.Result)

	fmt.Printf("Loaded %d ephemeris points from JPL Horizons\n", len(ephemeris))
	return ephemeris, nil
}

// parseHorizonsResult parses the JPL Horizons result text format
// This is a simplified stub - real implementation would parse the actual format
func parseHorizonsResult(result string) []EphemerisPoint {
	// Placeholder implementation
	// In reality, this would parse the complex JPL text format
	var points []EphemerisPoint

	// For demonstration, create some sample data
	baseTime := time.Date(2024, 1, 1, 0, 0, 0, 0, time.UTC)
	for i := 0; i < 100; i++ {
		points = append(points, EphemerisPoint{
			Time: baseTime.Add(time.Duration(i) * 24 * time.Hour),
			Position: Vector3D{
				X: MercurySemiMajorAxis * math.Cos(float64(i)*0.1),
				Y: MercurySemiMajorAxis * math.Sin(float64(i)*0.1),
				Z: 0,
			},
			Velocity: Vector3D{
				X: -47000 * math.Sin(float64(i)*0.1),
				Y: 47000 * math.Cos(float64(i)*0.1),
				Z: 0,
			},
		})
	}

	return points
}

// VSOP87 Implementation (simplified)
// VSOP87 is an analytical theory for planetary positions

// VSOP87Term represents a single term in the VSOP87 series
type VSOP87Term struct {
	A float64 // Amplitude
	B float64 // Phase
	C float64 // Frequency
}

// VSOP87Data holds the series terms for a planet's orbital elements
type VSOP87Data struct {
	L []VSOP87Term // Mean longitude terms
	B []VSOP87Term // Latitude terms
	R []VSOP87Term // Radius terms
}

// GetVSOP87Data returns simplified VSOP87 coefficients for Mercury
// Real implementation would include thousands of terms from VSOP87 catalog
func GetVSOP87Data() map[string]VSOP87Data {
	return map[string]VSOP87Data{
		"Mercury": {
			// Simplified terms - real VSOP87 has hundreds of terms per element
			L: []VSOP87Term{
				{4.40250710144, 0.0, 0.0},                        // L0
				{0.40989414977, 1.48302034195, 26087.9031415742}, // L1
				{0.05046294200, 4.47785489551, 52175.8062831484}, // L2
			},
			B: []VSOP87Term{
				{0.11737528961, 1.98357498767, 26087.9031415742},
				{0.02388076996, 5.03738959686, 52175.8062831484},
			},
			R: []VSOP87Term{
				{0.39528271651, 0.0, 0.0},
				{0.07834131818, 6.19233722598, 26087.9031415742},
				{0.00795525558, 2.95989690104, 52175.8062831484},
			},
		},
	}
}

// CalculateVSOP87Position computes planet position using VSOP87 theory
func CalculateVSOP87Position(planetName string, julianDay float64) Vector3D {
	data := GetVSOP87Data()[planetName]

	// Time since J2000.0 epoch in millennia
	T := (julianDay - 2451545.0) / 365250.0

	// Calculate orbital elements
	L := evaluateVSOP87Series(data.L, T) // Mean longitude
	B := evaluateVSOP87Series(data.B, T) // Latitude
	R := evaluateVSOP87Series(data.R, T) // Radius

	// Convert to Cartesian coordinates (simplified)
	x := R * math.Cos(B) * math.Cos(L) * AstronomicalUnit
	y := R * math.Cos(B) * math.Sin(L) * AstronomicalUnit
	z := R * math.Sin(B) * AstronomicalUnit

	return Vector3D{x, y, z}
}

// evaluateVSOP87Series evaluates a VSOP87 series for given time T
func evaluateVSOP87Series(terms []VSOP87Term, T float64) float64 {
	var result float64

	for _, term := range terms {
		result += term.A * math.Cos(term.B+term.C*T)
	}

	return result
}

// InterpolateEphemeris performs linear interpolation between ephemeris points
func InterpolateEphemeris(ephemeris []EphemerisPoint, targetTime time.Time) (Vector3D, Vector3D, error) {
	if len(ephemeris) < 2 {
		return Vector3D{}, Vector3D{}, fmt.Errorf("insufficient ephemeris data")
	}

	// Find surrounding points
	var before, after EphemerisPoint
	found := false

	for i := 0; i < len(ephemeris)-1; i++ {
		if ephemeris[i].Time.Before(targetTime) && ephemeris[i+1].Time.After(targetTime) {
			before = ephemeris[i]
			after = ephemeris[i+1]
			found = true
			break
		}
	}

	if !found {
		return Vector3D{}, Vector3D{}, fmt.Errorf("target time outside ephemeris range")
	}

	// Linear interpolation
	dt := after.Time.Sub(before.Time).Seconds()
	t := targetTime.Sub(before.Time).Seconds()
	factor := t / dt

	position := before.Position.Add(after.Position.Subtract(before.Position).Scale(factor))
	velocity := before.Velocity.Add(after.Velocity.Subtract(before.Velocity).Scale(factor))

	// Convert from km to m
	position = position.Scale(1000)
	velocity = velocity.Scale(1000)

	return position, velocity, nil
}

// Extension functions for future development

// CalculateForPlanet computes GR precession for an arbitrary planet.
func CalculateForPlanet(name string, mass, a, e, T float64) {
	fmt.Printf("\n=== Calculation for %s ===\n", name)
	result := CalculateGRPrecession(mass, a, e, T)
	fmt.Printf("GR precession per century: %.2f arcseconds\n", result.PerCentury)
}

// VisualizePrecession placeholder for future visualization features.
func VisualizePrecession() {
	// TODO: Add visualization code using plotting libraries
	fmt.Println("\n[Orbital precession visualization would go here]")
}

// DemoDataSources demonstrates the data source capabilities
func DemoDataSources() {
	fmt.Println("\n=== Data Sources Demo ===")

	// Demo VSOP87
	fmt.Println("VSOP87 Position Calculation:")
	julianDay := 2451545.0 // J2000.0 epoch
	position := CalculateVSOP87Position("Mercury", julianDay)
	fmt.Printf("Mercury position (J2000.0): %.3e, %.3e, %.3e m\n", position.X, position.Y, position.Z)

	// Demo JPL Horizons (would require network access)
	fmt.Println("\nJPL Horizons Integration available (requires network)")
	fmt.Println("Use LoadJPLHorizons() for real ephemeris data")

	fmt.Println("\nNote: Full data integration requires:")
	fmt.Println("- Network access for JPL Horizons API")
	fmt.Println("- Complete VSOP87 coefficient files")
	fmt.Println("- DE440 ephemeris files for highest accuracy")
}
