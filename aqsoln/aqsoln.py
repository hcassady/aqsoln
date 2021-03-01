import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import fsolve
import sqlite3
import molmass


class Calculator:

  def __init__(self, salt, check_bounds=True):
    self.salt = salt
    self.solute_molar_mass = molmass.Formula(salt).mass
    self.solvent_molar_mass = molmass.Formula('H2O').mass
    self._load_data()
    self._initialize_interpolator()

  def _load_data(self):
    query = """
      --begin-sql
      SELECT
        weight_percent,  -- (g solute / g total) x 100%
        temperature,     -- (deg C)
        density,         -- (g total / mL total)
        source
      FROM
        aqueous_solution_densities
      WHERE (salt = ?
        OR salt = 'None')
      AND temperature <= 100;
    """
    conn = sqlite3.connect('physical_data.db')
    self.raw_data = pd.read_sql(
      query, conn, params=(self.salt,))
    conn.close()

    self.sources = self.raw_data.source.unique()

  def _initialize_interpolator(self):
    """ Generate interpolator function.

    This creates an interpolator function of the LinearNDInterpolator type,
    which takes the data points (weight_percent, temperature) tuples as an input
    and outputs (density).
    """
    x = self.raw_data.weight_percent.to_numpy()
    y = self.raw_data.temperature.to_numpy()
    z = self.raw_data.density.to_numpy()
    self._interpolate = LinearNDInterpolator(
      list(zip(x, y)), z
    )

  # ---------------------------------
  # Weight percent section

  def weight_percent_to_density(self, weight_percent, temperature=20):
    """Weight percent to molality.

    Returns the density (g solution / mL solution) of the solution for a given
    weight percent (g solute / g solution x 100%) and temperature (deg C).
    """
    return self.interpolate(weight_percent, temperature)

  def weight_percent_to_mass_fraction(self, weight_percent, temperature=None):
    """Weight percent to mass fraction. (Temperature not required)

    Returns the mass fraction (g solute / g solution) of the solution for a
    given weight percent (g solute / g solution x 100%).

    The temperature is not required, but is included as an allowed argument to
    allow the argument to be passed to (and ignored by) the function).
    """
    return weight_percent / 100

  def weight_percent_to_molar(self, weight_percent, temperature=20):
    """Weight percent to molality.

    Returns the molality (mol solute / kg solvent) of the solution for a given
    weight percent (g solute / g solution) and temperature (deg C).
    """
    w = self.weight_percent_to_mass_fraction(weight_percent, temperature)
    rho = self.weight_percent_to_density(weight_percent, temperature)   # g/mL
    Mw = self.solute_molar_mass                                         # g/mol
    molar_concentration = ((w * rho) / Mw)                              # mol/mL
    molar_concentration = molar_concentration * 1000                    # mol/L
    return molar_concentration

  def weight_percent_to_molal(self, weight_percent, temperature=20):
    """Weight percent to molality.

    Returns the molality (mol solute / kg solvent) of the solution for a given
    weight percent (g solute / g solution x 100%) and temperature (deg C).
    """
    w = self.weight_percent_to_mass_fraction(weight_percent, temperature)
    Mw = self.solute_molar_mass / 1000  # kg/mol
    molality = w / ((1 - w) * Mw)      # mol solute / kg solvent
    return molality

  def weight_percent_to_mole_fraction(self, weight_percent, temperature=None):
    """Weight percent to mole fraction.

    Returns the mole fraction (mol solute / mol total) of the solution for a
    given weight percent (g solute / g solution x 100%).

    The temperature is not required, but is included as an allowed argument to
    allow the argument to be passed to (and ignored by) the function).
    """
    w = weight_percent                         # Weight percent of solute
    Mw = self.solute_molar_mass                # Molar mass of solute, g/mol
    M_aq = self.solvent_molar_mass             # Molar mass of water,  g/mol
    M_bar = ((w / Mw) + ((1-w) / M_aq))**-1    # Average molar mass,   g/mol
    x = (w * M_bar) / Mw                       # Mole fraction
    return x

  # ---------------------------------
  # Molar section

  def molar_to_weight_percent(self, molarity, temperature=20):
    """Molarity to weight percent.

    Returns the weight percent (g solute / g solution * 100%) of the solution
    for a given molarity (mol solute / L solution) and temperature (deg C).
    """
    molarity_per_mL = molarity / 1000    # mol/mL

    def zero_function(weight_percent):
      zero = (
        ((molarity_per_mL * self.solute_molar_mass) /
         self.weight_percent_to_density(
           weight_percent, temperature)) - weight_percent / 100
      )
      return zero
    # Initial guess based on density of 1 g/mL
    initial_guess = molarity_per_mL * self.solute_molar_mass * 100
    weight_percent = fsolve(zero_function, [initial_guess])
    return weight_percent

  def molar_to_mass_fraction(self, molarity, temperature=20):
    """Molar concentration to mass fraction.

    Returns the mass fraction (g solute / g solvent) of the solution for a given
    molarity (mol solute / L solution) and temperature (deg C).
    """
    return self.molar_to_weight_percent(molarity, temperature) / 100

  def molar_to_density(self, molarity, temperature=20):
    """Molar concentration to density.

    Returns the solution density (g /mL) of the solution for a given molarity
    (mol solute / L solution) and temperature (deg C).
    """
    return self.weight_percent_to_density(
      self.molar_to_weight_percent(molarity, temperature), temperature
    )

  def molar_to_molal(self, molarity, temperature=20):
    """Molar concentration to molality.

    Returns the molality (mol solute / kg solvent) of the solution for a given
    molarity (mol solute / L solution) and temperature (deg C).
    """
    return self.weight_percent_to_molal(
      self.molar_to_weight_percent(molarity, temperature), temperature
    )

  def molar_to_mole_fraction(self, molarity, temperature=20):
    """Molar concentration to mole fraction.

    Returns the mole fraction (mol solute / mol total) of the solution for a
    given molarity (mol solute / L solution) and temperature (deg C).
    """
    return self.weight_percent_to_mole_fraction(
      self.molar_to_weight_percent(molarity, temperature), temperature
    )

  # ---------------------------------
  # Molal section

  def molal_to_weight_percent(self, molality, temperature=None):
    """Molal concentration to weight percent.

    Returns weight percent of the solution for a given molality (mol solute / kg
    solvent).

    The temperature is not required, but is included as an allowed argument to
    allow the argument to be passed to (and ignored by) the function).
    """
    b = molality                          # mol solute / kg solution
    Mw = self.solute_molar_mass / 1000    # kg solute / mol solute
    w = (1 + 1 / (b * Mw))**-1
    return w

  def molal_to_mass_fraction(self, molality, temperature=None):
    """Molal concentration to mass fraction.

    Returns mass fraction of the solution for a given molality (mol solute / kg
    solvent).

    The temperature is not required, but is included as an allowed argument to
    allow the argument to be passed to (and ignored by) the function).
    """
    return self.molal_to_weight_percent(molality) / 100

  def molal_to_density(self, molality, temperature=20):
    """Molal concentration to density.

    Returns density (g solution / mL solution) of the solution for a given
    molality (mol solute / kg solvent) and temperature (deg C).
    """
    return self.weight_percent_to_density(
      self.molal_to_weight_percent(molality), temperature
    )

  def molal_to_molar(self, molality, temperature=20):
    """Molal concentration to molar concentration.

    Returns molarity (mol solute / L solution) of the solution for a given
    molality (mol solute / kg solvent) and temperature (deg C).
    """
    return self.weight_percent_to_molar(
      self.molal_to_weight_percent(molality), temperature
    )

  def molal_to_mole_fraction(self, molarity, temperature=None):
    """Molal concentration to mole fraction.

    Returns the mole fraction (mol solute / mol total) of the solution for a
    given molality (mol solute / kg solvent) and temperature (deg C).

    The temperature is not required, but is included as an allowed argument to
    allow the argument to be passed to (and ignored by) the function).
    """
    return self.weight_percent_to_mole_fraction(
      self.molal_to_weight_percent(molarity)
    )
