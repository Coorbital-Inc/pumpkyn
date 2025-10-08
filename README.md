# 🎃Pumpkyn Toolbox

**Pumpkyn** is an open-source MATLAB toolbox for exploring and analyzing tulip-shaped orbits in the Earth–Moon Circular Restricted Three-Body Problem (CR3BP). It provides a hands-on environment for generating, propagating, and diagnosing periodic orbit families that are relevant to cislunar navigation, communications, and science missions.

---

### 🧩 Installation
1. Download `Pumpkyn.mltbx` or clone this repository.  
2. In MATLAB, double-click the `.mltbx` file or run:

```matlab
    matlab.addons.install('Pumpkyn.mltbx');
```

---

### Using MATLAB Live Scripts
To get started quickly, open the **interactive Live Script tutorial** included in the `/docs` folder:  

```matlab
    open('docs/intro.mlx')
```
---

### 🚀 Features
- Generate and analyze tulip- and pumpkin-shaped periodic orbits in the Earth–Moon CR3BP  
- Propagate and visualize trajectories in both rotating and inertial reference frames  
- Assess orbit stability, energy levels, and maintenance ΔV requirements  
- Explore zero-velocity surfaces, invariant manifolds, and low-energy transfer pathways  
- Perform orbit continuation and design cislunar constellations for communications or SDA missions  
- Evaluate surface access, coverage, and navigation geometry using built-in analysis tools  
- Includes fully documented Live Script tutorials for hands-on learning and research  

---
./toolbox/docs/
### 📚 Tutorials
The toolbox includes guided Live Script tutorials located in the `/docs` folder:
| Tutorial | Description |
|-----------|--------------|
| [**Extracting Tulip-Shaped Orbits and Properties**](./toolbox/docs/md/intro_tulip.md) | Generates tulip orbits in the Earth–Moon CR3BP. Retrieves full orbit families for a chosen petal count and computes properties such as Jacobi constant, stability index, perilune/apolune altitudes, and lunar occultation. |
| [**Creating a Pumpkin Orbit**](./toolbox/docs/md/intro_pumpkyn.md) | Generates, propagates, and visualizes a 13-petal pumpkin-shaped orbit. Computes stability index, Jacobi constant, and altitude extremes to assess orbit stability and suitability for cislunar operations. |
| [**Zero-Velocity Surface Analysis**](./toolbox/docs/md/zeroVelocity.md) | Visualizes zero-velocity surfaces for a chosen Jacobi constant and overlays a 9-petal tulip orbit (C = 3.1777). Illustrates how energy level determines accessibility through L1 and confinement near the Moon. |
| [**Station-Keeping Estimation Demonstration**](./toolbox/docs/md/stationKeepingEst.md) | Estimates ΔV requirements for maintaining a pumpkin-shaped orbit. Applies linear-stability theory to compute upper and lower ΔV bounds for various correction intervals, providing insight into maneuver frequency and annual maintenance costs. |
| [**Pseudo-Arclength Continuation (Tulips → Pumpkins)**](./toolbox/docs/md/tulipContinuation.md) | Extends a tulip orbit using pseudo-arclength continuation. Starting from a seed orbit, the routine traces new family members across increasing sidereal periods, visualizing geometric evolution from tulip to pumpkin. |
| [**Finding Unstable Manifolds**](./toolbox/docs/md/manifolds.md) | Computes and visualizes invariant manifolds departing from a periodic tulip orbit. Highlights dynamical pathways that enable low-energy transfers and natural escape/capture trajectories. |
| [**Designing a Tulip-Shaped Constellation**](./toolbox/docs/md/tulipContinuation.md) | Builds a time-phased constellation along a multi-petal tulip orbit in the CR3BP rotating frame, illustrating constellation design for persistent lunar coverage or SDA missions. |
| [**Minimum-Energy Transfers**](./toolbox/docs/md/minEnergy.md) | Estimates the minimum ΔV required for a tulip-orbit spacecraft to reach the lunar surface. Maps energy contours to reveal optimal access regions and low-energy transfer opportunities. |
| [**Coordinate Transformations: J2000 ↔ CR3BP**](./toolbox/docs/md/j2k2cr3bp.md) | Demonstrates conversion between Moon-centered inertial and barycentric rotating frames. Bridges classical two-body dynamics with multi-body CR3BP analysis for integrated mission design. |
| [**Lunar Surface Coverage Analysis**](./toolbox/docs/md/surfCoverage.md) | Simulates a tulip-orbit constellation and computes lunar surface coverage based on a minimum elevation angle. Plots maximum uncovered-gap times to evaluate communications and sensing performance. |
| [**Lunar South Pole GDOP Analysis**](./toolbox/docs/md/intro_gdop.md) | Computes geometric dilution of precision (GDOP) at the lunar south pole for a 7-petal, 6:5-resonant tulip constellation. Evaluates navigation geometry and PNT performance over time. |
---

## 🛰 References  

- Koblick, D. C., and McCarthy, B. (2026). *Station-Keeping Techniques for Sidereal Resonant Tulip-Shaped Three-Body Orbits.* AIAA SciTech Forum, Orlando, FL.  
- Koblick, D. C., and Kelly, P. (2025). *Novel Three-Body Tulip-Shaped Orbit Families for Lunar Missions.* *Journal of the Astronautical Sciences*, 72(4): 32.  
- Koblick, D. C. (2023). *Novel Tulip-Shaped Three-Body Orbits for Cislunar Space Domain Awareness Missions.* Proceedings of the Advanced Maui Optical and Space Surveillance (AMOS) Technologies Conference, Maui, HI.  

---

## Coorbital Inc.

Source code is actively maintained by Coorbital Inc, with ongoing updates to expand functionality and improve usability. As new mission needs and research opportunities emerge, we welcome collaboration and can provide custom features, tailored analyses, or specialized extensions upon request. For inquiries, suggestions, or feature requests, please reach out.

---

**Disclaimer:** Pumpkyn Toolbox is provided "as is," without warranty of any kind, express or implied. Coorbital Inc. does not guarantee the accuracy or reliability of any results produced and is not responsible for any damages, losses, or consequences arising from its use. See our [license](./LICENSE) file for additional information.
