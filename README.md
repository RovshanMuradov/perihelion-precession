# Mercury Perihelion Precession Calculator

A Go implementation that calculates and simulates the perihelion precession of Mercury using Einstein's General Theory of Relativity, reproducing one of the most famous confirmations of Einstein's revolutionary theory.

## 🌟 Scientific Background

### The Mercury Problem: A Crisis in Classical Physics

In the 19th century, astronomers discovered a puzzling anomaly in Mercury's orbit that classical Newtonian mechanics couldn't explain. Mercury's perihelion (the point in its orbit closest to the Sun) was precessing (rotating) at a rate of **574.1 arcseconds per century**. 

Classical gravitational theory, accounting for perturbations from other planets, could only explain **531 arcseconds per century**, leaving a mysterious **43 arcseconds per century** unexplained. This discrepancy was small but persistent, creating a crisis in our understanding of gravity.

### Einstein's Revolutionary Solution

In 1915, Albert Einstein published his General Theory of Relativity, which described gravity not as a force, but as the curvature of spacetime caused by mass and energy. When Einstein applied his new theory to Mercury's orbit, he derived the famous formula:

```
Δφ = 6πGM / (c²a(1-e²))
```

Where:
- `Δφ` = perihelion advance per orbit
- `G` = gravitational constant
- `M` = mass of the Sun
- `c` = speed of light
- `a` = semi-major axis of Mercury's orbit
- `e` = eccentricity of Mercury's orbit

This formula predicted exactly **42.98 arcseconds per century** - matching the unexplained remainder with stunning precision! This was the first major confirmation of General Relativity and marked the beginning of modern physics.

## 🎯 Project Purpose

This project demonstrates and calculates the Mercury perihelion precession effect using multiple approaches:

1. **Analytical Calculation**: Direct application of Einstein's GR formula
2. **N-body Simulation**: Full physics simulation with gravitational interactions
3. **Planetary Analysis**: Detailed breakdown of contributions from different planets
4. **Data Integration**: Real astronomical data from JPL Horizons and VSOP87

### Why This Matters

- **Historical Significance**: Reproduces one of the most important physics discoveries
- **Educational Value**: Demonstrates the transition from classical to modern physics
- **Computational Physics**: Shows how theoretical physics translates to numerical simulations
- **Astronomical Precision**: Uses real orbital parameters and gravitational data

## 🚀 Features

### Phase 1: Analytical Calculation
- **Einstein's Formula**: Direct implementation of the GR perihelion advance equation
- **Precise Results**: Calculates the famous ~43 arcseconds per century
- **Multiple Units**: Results in radians, arcseconds, per orbit, and per century
- **Comparison Data**: Shows observed vs. theoretical values

### Phase 2: N-body Simulation  
- **Full Physics Engine**: Gravitational forces between all planets
- **RK4 Integration**: 4th-order Runge-Kutta for numerical accuracy
- **Real-time Tracking**: Monitors precession during simulation
- **Orbital Elements**: Complete Keplerian orbital parameter tracking

### Phase 3: Data Sources
- **JPL Horizons API**: Real NASA/JPL astronomical data
- **VSOP87 Theory**: High-precision planetary position calculations
- **Historical Data**: Comparison with observed astronomical measurements

### Phase 4: Advanced Analysis
- **Parallel Processing**: Concurrent analysis of planetary influences
- **Adaptive Time Stepping**: Automatic accuracy optimization
- **Decomposition**: Separates GR effects from classical perturbations
- **Statistical Analysis**: Comprehensive results breakdown

## 📊 Visualization

The project generates beautiful SVG visualizations:

- **`mercury_orbit.svg`**: Mercury's orbital path showing precession
- **`precession_plot.svg`**: Perihelion longitude vs. time graph  
- **`rosette_pattern.svg`**: Artistic rosette pattern from orbital precession

## 🔬 Scientific Accuracy

### Expected Results
- **General Relativity Effect**: ~43 arcseconds/century
- **Total Observed Precession**: 574.1 arcseconds/century
- **Classical Perturbations**: ~531 arcseconds/century

### Planetary Contributions (per century)
- **Venus**: ~277.8 arcseconds
- **Jupiter**: ~153.6 arcseconds  
- **Earth**: ~90.0 arcseconds
- **Other planets**: ~9.6 arcseconds

## 🛠️ Installation & Usage

### Prerequisites
```bash
# Requires Go 1.19+
go version
```

### Running the Simulation
```bash
# Clone and run
git clone <repository>
cd perihelion-precession
go run main.go

# Or build and run
go build -o perihelion main.go
./perihelion
```

### Interactive Menu
```
=== Mercury Perihelion Precession Analysis ===
1. Analytical calculation (Phase 1)     # Fast, direct formula application
2. N-body simulation (Phase 2)          # Full physics simulation  
3. Data sources demo (JPL Horizons)     # Real astronomical data
4. Parallel planetary analysis          # Advanced multi-planet analysis
```

### Example Output
```
General Relativity Precession Results:
Precession per orbit: 5.0155e-07 radians (0.1034 arcseconds)
Mercury completes 415.2 orbits per century
Total GR precession: 42.98 arcseconds per century

Historical comparison:
- Observed total: 574.1 arcseconds/century
- Classical theory: 531.1 arcseconds/century
- Einstein's prediction: 42.98 arcseconds/century
- Difference: 0.02 arcseconds/century ✓
```

## 🧮 Physics Implementation

### Core Equation
```go
// Einstein's General Relativity formula
numerator := 6 * π * G * M_sun
denominator := c² * a * (1 - e²)
precession_per_orbit := numerator / denominator
```

### N-body Gravitational Forces
```go
// Newton's law of universal gravitation
F = G * m1 * m2 / r²
```

### Numerical Integration
- **Method**: 4th-order Runge-Kutta (RK4)
- **Adaptive Stepping**: Automatic error control
- **Conservation**: Energy and angular momentum monitoring

## 📚 References

- Einstein, A. (1915). "Explanation of the Perihelion Motion of Mercury from General Relativity Theory"
- Clemence, G. M. (1947). "The Relativity Effect in Planetary Motions"
- Will, C. M. (2014). "The Confrontation between General Relativity and Experiment"
- JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/
- VSOP87 Theory: Bureau des Longitudes, Paris

## 🎖️ Historical Context

This calculation represents one of science's greatest triumphs:

> *"Scarcely anyone who truly understands [General Relativity] can escape from its magic."* - Albert Einstein

The Mercury perihelion precession was:
- **1859**: First observed by Urbain Le Verrier
- **1915**: Explained by Einstein's General Relativity  
- **1919**: Confirmed alongside eclipse observations
- **Today**: Measured with extraordinary precision by space missions

This project lets you experience firsthand the moment when Einstein's theory proved that Newton's universe, while magnificent, was incomplete.

---

*"The most beautiful thing we can experience is the mysterious."* - Albert Einstein