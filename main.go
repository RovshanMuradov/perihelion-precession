package main

import (
	"fmt"
	"math"
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

// Precession Analysis Functions

// ConvertStateToOrbitalElements converts Cartesian state (position, velocity) to Keplerian orbital elements.
// This is a simplified version focusing on the key elements needed for precession analysis.
func ConvertStateToOrbitalElements(planet Planet, centralMass float64) OrbitalElements {
	r := planet.Position.Magnitude()
	v := planet.Velocity.Magnitude()

	// Specific orbital energy
	specificEnergy := 0.5*v*v - GravitationalConstant*centralMass/r

	// Semi-major axis from energy: E = -GM/(2a)
	semiMajorAxis := -GravitationalConstant * centralMass / (2 * specificEnergy)

	// Angular momentum vector
	angularMomentum := Vector3D{
		X: planet.Position.Y*planet.Velocity.Z - planet.Position.Z*planet.Velocity.Y,
		Y: planet.Position.Z*planet.Velocity.X - planet.Position.X*planet.Velocity.Z,
		Z: planet.Position.X*planet.Velocity.Y - planet.Position.Y*planet.Velocity.X,
	}
	h := angularMomentum.Magnitude()

	// Eccentricity from angular momentum and energy
	eccentricity := math.Sqrt(1 + 2*specificEnergy*h*h/(GravitationalConstant*GravitationalConstant*centralMass*centralMass))

	// For simplified 2D case, set other elements to reasonable defaults
	// In a full implementation, these would be calculated from the 3D vectors
	inclination := 0.0 // Assume coplanar for now
	longitudeOfNode := 0.0

	// True anomaly approximation (simplified)
	trueAnomaly := math.Atan2(planet.Position.Y, planet.Position.X)

	// Argument of periapsis (simplified)
	argumentOfPeriapsis := 0.0 // Would need eccentricity vector calculation

	// Mean anomaly from true anomaly (simplified)
	meanAnomaly := trueAnomaly // Approximate for low eccentricity

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

// PrecessionTracker tracks perihelion precession over time.
type PrecessionTracker struct {
	InitialLongitude  float64   // Initial perihelion longitude (radians)
	PreviousLongitude float64   // Previous perihelion longitude (radians)
	TotalPrecession   float64   // Total accumulated precession (radians)
	PrecessionHistory []float64 // History of perihelion longitudes
	TimeHistory       []float64 // Corresponding time points
	OrbitalCount      int       // Number of completed orbits
}

// NewPrecessionTracker creates a new precession tracker.
func NewPrecessionTracker(initialLongitude float64) *PrecessionTracker {
	return &PrecessionTracker{
		InitialLongitude:  initialLongitude,
		PreviousLongitude: initialLongitude,
		TotalPrecession:   0.0,
		PrecessionHistory: []float64{initialLongitude},
		TimeHistory:       []float64{0.0},
		OrbitalCount:      0,
	}
}

// Update adds a new measurement to the precession tracker.
func (pt *PrecessionTracker) Update(time float64, longitude float64) {
	// Handle angle wrapping (longitude goes from 0 to 2π)
	longitudeDiff := longitude - pt.PreviousLongitude

	// Correct for angle wrapping
	if longitudeDiff > math.Pi {
		longitudeDiff -= 2 * math.Pi
	} else if longitudeDiff < -math.Pi {
		longitudeDiff += 2 * math.Pi
	}

	pt.TotalPrecession += longitudeDiff
	pt.PreviousLongitude = longitude
	pt.PrecessionHistory = append(pt.PrecessionHistory, longitude)
	pt.TimeHistory = append(pt.TimeHistory, time)

	// Count completed orbits (rough approximation)
	totalLongitudeChange := math.Abs(pt.TotalPrecession)
	if totalLongitudeChange > float64(pt.OrbitalCount+1)*2*math.Pi {
		pt.OrbitalCount++
	}
}

// GetPrecessionRate calculates the precession rate in arcseconds per year.
func (pt *PrecessionTracker) GetPrecessionRate(timeSpanSeconds float64) float64 {
	if timeSpanSeconds == 0 {
		return 0
	}

	// Convert from radians/second to arcseconds/year
	radiansPerSecond := pt.TotalPrecession / timeSpanSeconds
	secondsPerYear := 365.25 * 24 * 3600
	return radiansPerSecond * secondsPerYear * RadiansToArcseconds
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

	// Calculate current orbital elements
	elements := ConvertStateToOrbitalElements(mercury, pa.CentralMass)
	longitude := CalculatePerihelionLongitude(elements)

	// Update main tracker
	pa.Tracker.Update(state.Time, longitude)
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
	fmt.Println("=== Mercury Perihelion Precession Calculation ===")
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
