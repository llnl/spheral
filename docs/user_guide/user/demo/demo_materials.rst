.. _material_properties_generation:

=============================
Units and material properties
=============================

The next section of ``Sedov-demo.py`` describes the material models we are going to use in this calculation::

  #-------------------------------------------------------------------------------
  # Material properties.
  #-------------------------------------------------------------------------------
  units = MKS()
  eos = GammaLawGas(gamma = gamma,
                    mu = mu,
                    constants = units)

The first line descibes our physical units choice: here we are going to use standard SI units (meters, kilograms, and seconds).  Spheral is different from many physics modeling codes in that we do not use a single set of physical units in the code, but rather allow the user to specify the units of their choice appropriate for the problem at hand.  In general you may specify any choice by explicitly setting the simulation units of length, mass, and time to a ``PhysicalConstants`` object like so::

  units = PhysicalConstants(1.0,   # unit length in meters
                            1.0,   # unit mass in kilograms
                            1.0,   # unit time in seconds)

In ``Sedov-demo.py`` we have chosen the standard SI (International System of Units) units, which is the equivalent of explicitly setting each of the unit length, mass, and time to unity as in the example above. We have used a Spheral shorthand for this choice (``MKS`` units, for [meter, kilogram, seconds]). A generally good rule of thumb for is that you want to choose units such that the scale of most of the quantities in your model are roughly unity. See the description of :ref:`PhysicalConstants` for more information.

The next line constructs the equation of state for the single material in this problem; in this case we using a simple gamma-law gas. The arguments are the choice of :math:`\gamma` (ratio of specific heats), :math:`\mu` the mean molecular weight, and the choice of units.  Note as mentioned earlier, one consequence of Spheral being a Python module and not watching everything entered into the Python prompt is that we need to pass the choice of units into the equation of state as shown here.
