# Atom module

## Helium 
The `Helium` module provides atomic data and calculations specifically for Helium atoms.

### Methods
1. **`Initialisation`**
   - Default physical parameters:
     - Mass
     - Atomic wavelength
     - Atomic width
     - Atomic omega
     - Scattering length
   - Default magnetic parameters:
     - Lande g-factor
     - Zeeman state
   
2. **`get_alpha(self, wavelength=1550e-9, unit="SI")`**
   - Computes the polarizability.
   - Parameters:
     - `wavelength`: Laser wavelength (in meters)
     - `unit`: 'SI' (default) or 'au' (atomic units)
   - Returns:
     - Polarizability (in SI / au)
   
3. **`get_scattering_rate(self, intensity, wavelength=1550e-9)`**
   - Computes the scattering rate at given intensity & detuning.
   - Parameters:
     - `intensity`: Laser intensity (in W/m^2)
     - `wavelength`: Laser wavelength (in meters)
   - Returns:
     - Scattering rate (in s^-1)
   
4. **`convert_speed_to_lattice_momentum(self, v, wavelength=1064.0, theta=166)`**
   - Converts atom speed (mm/s) to momentum in lattice momentum unit klatt.
   - Parameters:
     - `v`: Speed of the atom in millimeter per second
     - `wavelength`: Wavelength of the lattice (in nanometers), default: 1064
     - `theta`: Angle in degrees between the two beams, default: 166
   - Returns:
     - Relative speed compared to lattice momentum
   
5. **`convert_lattice_momentum_to_speed(self, k, wavelength=1064.0, theta=166)`**
   - Converts atom speed (mm/s) to momentum in hbar unit (SI).
   - Parameters:
     - `k`: Momentum in lattice momentum unit klatt
     - `wavelength`: Wavelength of the lattice (in nanometers), default: 1064
     - `theta`: Angle in degrees between the two beams, default: 166
   - Returns:
     - Speed compared to lattice momentum
