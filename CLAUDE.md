# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Go project that calculates the perihelion precession of Mercury using Einstein's General Theory of Relativity. The project demonstrates the famous ~43 arcseconds per century effect that confirmed Einstein's theory.

**Current Status**: Phase 1 & 2 COMPLETE! Full N-body simulation with visualization operational.

## Architecture

### Phase 1 (Current Implementation)
The codebase consists of a single main.go file with the following components:

- **Physical constants**: Gravitational constant, speed of light, astronomical unit, solar mass, and Mercury's orbital parameters
- **PrecessionResult struct**: Contains calculated precession values in different units and time scales
- **CalculateGRPrecession()**: Core function implementing the General Relativity formula Δφ = 6πGM / (c²a(1-e²))
- **PlanetaryContributions struct**: Structured approach to planetary perturbations
- **Display functions**: Modular output functions following Go best practices

### Phase 2 (Implemented N-body Simulation)
Complete N-body simulation system with real-time precession tracking:

- **Planet struct**: Mass, position, velocity, orbital elements ✅
- **Vector3D**: 3D vector operations with methods ✅
- **OrbitalElements**: Keplerian orbital parameters (a, e, i, Ω, ω, M) ✅  
- **Physics engine**: Gravitational forces, RK4 integration ✅
- **Precession analyzer**: Real-time tracking and decomposition ✅
- **Data sources**: JPL Horizons API, VSOP87 integration ✅
- **Visualization**: SVG orbital plots, precession graphs, rosette patterns ✅

## Development Commands

**Run the application:**
```bash
go run main.go
# Choose option:
# 1: Analytical calculation (fast, ~43 arcseconds result)
# 2: N-body simulation (slow, full physics with visualization)  
# 3: Data sources demo (JPL Horizons, VSOP87)
```

**Build executable:**
```bash
go build -o perihelion main.go
```

**Code quality:**
```bash
go fmt ./...
go vet ./...
golangci-lint run  # Static analysis (config issues exist with current version)
```

**Testing (when implemented):**
```bash
go test ./...
go test -v ./...     # Verbose output
go test -cover ./... # Coverage report
```

**Module management:**
```bash
go mod tidy
go mod verify
```

## Code Standards

### Go Idioms Applied
- Exported functions use PascalCase (CalculateGRPrecession)
- Constants use descriptive names (GravitationalConstant, not G)
- Structs have proper documentation comments
- Functions are broken down into focused, single-purpose units
- Error handling follows Go conventions (when applicable)

### Naming Conventions
- Physical constants: Full descriptive names (SpeedOfLight, SolarMass)
- Functions: Verb-noun patterns (CalculateGRPrecession, PrintResults)
- Structs: Noun phrases (PrecessionResult, PlanetaryContributions)
- Files: snake_case for modules, but currently single file

## Phase 2 Implementation Plan

### Required Data Structures
```go
type Planet struct {
    Name     string
    Mass     float64
    Position Vector3D
    Velocity Vector3D
    Elements OrbitalElements
}

type Vector3D struct {
    X, Y, Z float64
}

type OrbitalElements struct {
    SemiMajorAxis      float64 // a
    Eccentricity       float64 // e
    Inclination        float64 // i
    LongitudeOfNode    float64 // Ω
    ArgumentOfPeriapsis float64 // ω
    MeanAnomaly        float64 // M
}
```

### Success Criteria for Phase 2
- Reproduce known planetary contributions with ±5% accuracy
- Venus: ~277.8" (expect 265-290")
- Jupiter: ~153.6" (expect 145-160")
- Earth: ~90.0" (expect 85-95")

### Performance Considerations
- Use Runge-Kutta 4th order integration
- Adaptive time stepping for accuracy
- Parallel computation for different planets
- Caching for repeated calculations

## Key Implementation Details

The project uses the exact General Relativity formula for perihelion advance and compares results against the historically significant observed value of 574.1 arcseconds per century. The calculation separates:

1. Classical planetary perturbations (~531 arcseconds/century)
2. General Relativity correction (~43 arcseconds/century)

The success criterion is achieving the famous ~43 arcseconds result that validated Einstein's theory.

## Implemented Visualizations ✅
- **mercury_orbit.svg**: Mercury's orbital path showing precession effect
- **precession_plot.svg**: Perihelion longitude vs time graph  
- **rosette_pattern.svg**: Beautiful rosette pattern from multiple orbits
- Real-time progress tracking during simulation
- Comprehensive results analysis and comparison

## Advanced Features Available
- RK4 numerical integration for accuracy
- Configurable simulation duration (recommended: 10-100 years)
- Adaptive output intervals for performance
- Real planetary masses and orbital parameters
- JPL Horizons API integration (network dependent)
- VSOP87 analytical position calculations
- Full solar system (Sun + 6 planets) simulation