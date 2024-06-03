# Monte-Carlo Electron Transport Simulation

This project simulates the transport of electrons in a cylindrical Geiger counter filled with Helium gas. The simulation models the complex dynamics of electron scattering under an electric field, incorporating the physics of electron-atom interactions and the effects of the Lorentz force.

## Simulation Details

### Overview

The Monte-Carlo simulation handles the probabilistic nature of quantum scattering and kinetic energy changes due to electron acceleration in an electric field. This challenge is addressed using the null-collision method, which effectively balances the probabilistic time steps for scattering events with the necessity of a constant, discrete simulation time step.

### Features

- **Dynamic Simulation**: Start the simulation with variable charges and voltage using the `run(charges, voltage)` function.
- **Geometry**: Simulates a cylindrical Geiger counter with precise dimensions for the anode and cathode, oriented along the z-axis.
- **Scattering Mechanisms**: Includes both elastic and inelastic scattering processes, with the latter leading to ionization and multiplication of electrons.
- **Collision Recording**: Records a subset of all collision events to manage memory and performance during simulation.

### Technical Specifications

- **Input Parameters**:
  - `charges`: List containing the initial positions of electrons to be transported.
  - `voltage`: Electric potential applied to the anode to simulate the electric field inside the Geiger counter.
- **Simulation Output**:
  - List of electron interaction locations in 3D for visualization.
  - Total transport times for each electron.

### Null-Collision Method

This method involves simulating both real and "null" collisions to maintain a consistent simulation time step. It ensures that the physical behavior of electrons under varying kinetic energies is accurately modeled without the need for continuous adjustments to the time step based on energy changes.

## Results and Output

The simulation generates detailed logs of electron paths and interactions within the Geiger counter. Visualization scripts are provided to help interpret the complex dynamics observed and to draw scatter plots of interaction locations.
