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
    self._define_function_dictionaries()

  def __call__(self, *args):
    input_value = args[0]
    from_ = self.abbreviations[args[-2]]
    to_ = self.abbreviations[args[-1]]
    # If no temperature given, default to 20 °C
    if len(args)==3:
      temperature = 20
    else:
      temperature = args[1]
    result = self._from_weight_percent_functions[to_](
      self._to_weight_percent_functions[from_](
        input_value, temperature
      ), temperature
    )
    return result

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

  def _define_function_dictionaries(self):
    self._to_weight_percent_functions = {
      'molar':          self.molar_to_weight_percent,
      'molal':          self.molal_to_weight_percent,
      'mass_fraction':  self.mass_fraction_to_weight_percent,
      'weight_percent': self.weight_percent_to_weight_percent,
      'mole_fraction':  self.mole_fraction_to_weight_percent,
      'mole_percent':   self.mole_percent_to_weight_percent,
      'density':        self.density_to_weight_percent
      }

    self._from_weight_percent_functions = {
      'molar':          self.weight_percent_to_molar,
      'molal':          self.weight_percent_to_molal,
      'mass_fraction':  self.weight_percent_to_mass_fraction,
      'weight_percent': self.weight_percent_to_weight_percent,
      'mole_fraction':  self.weight_percent_to_mole_fraction,
      'mole_percent':   self.weight_percent_to_mole_percent,
      'density':        self.weight_percent_to_density
      }

    self.abbreviations = {
      'molar':          'molar',
      'M':              'molar',
      'c':              'molar',
      'molal':          'molal',
      'm':              'molal',
      'b':              'molal',
      'mass_fraction':  'mass_fraction',
      'w':              'mass_fraction',
      'weight_percent': 'weight_percent',
      'wt%':            'weight_percent',
      'mole_fraction':  'mole_fraction',
      'x':              'mole_fraction',
      'mole_percent':   'mole_percent',
      'mol%':           'mole_percent',
      'density':        'density',
      'rho':            'density',
      'ρ':              'density'
    }

  # ---------------------------------
  # To and from section
  def weight_percent_to_weight_percent(self, weight_percent, temperature=None):
    return weight_percent

  # ---------------------------------
  # From section

  def weight_percent_to_density(self, weight_percent, temperature=20):
    """Weight percent to molality.

    Returns the density (g solution / mL solution) of the solution for a given
    weight percent (g solute / g solution x 100%) and temperature (deg C).
    """
    return self._interpolate(weight_percent, temperature)

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

  def weight_percent_to_mole_percent(self, weight_percent, temperature=None):
    return self.weight_percent_to_mole_fraction(weight_percent) * 100

  # ---------------------------------
  # To section

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

  def mass_fraction_to_weight_percent(self, mass_fraction, temperature=None):
    return mass_fraction * 100

  def mole_fraction_to_weight_percent(self, mole_fraction, temperature=None):
    x = mole_fraction
    Mw = self.solute_molar_mass
    M_aq = self.solvent_molar_mass
    M_bar = x * Mw + (1 - x) * M_aq
    w = (x * Mw) / M_bar
    return w * 100

  def mole_percent_to_weight_percent(self, mole_percent, temperature=None):
    mole_fraction = mole_percent / 100
    return self.mole_fraction_to_weight_percent(mole_fraction)

  def density_to_weight_percent(self, density, temperature=20):
    def zero_function(weight_percent):
      return (self.weight_percent_to_density(weight_percent, temperature) -
              density)
    weight_percent = fsolve(zero_function, [5])
    return weight_percent
