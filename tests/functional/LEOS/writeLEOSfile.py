#-------------------------------------------------------------------------------
# Sample a Spheral EOS and write a LEOS tabulated file from it
#-------------------------------------------------------------------------------
def writeLEOSfile(eos,                 # A Spheral EquationOfState to be sampled
                  rhoMin,              # Minimum mass density for table
                  rhoMax,              # Maximum mass density for table
                  Tmin,                # Minimum temperature for table
                  Tmax,                # Maximum temperature for table
                  nrho,                # Number of samples along mass density axis
                  nT,                  # Number of samples along temperature axis
                  basename,            # Output file base name (makes basename.ascii and basename.data)
                  name,                # String name for EOS provided in file header section
                  eosNumber,           # Number for EOS in LEOS
                  rho0 = None,         # Reference density -- if None attempt to query eos
                  T0 = None,           # Reference temperature -- if None attempt to query eos
                  atomicWeight = None, # Atomic weight -- if None attempt to query eos
                  atomicNumber = None, # Atomic number -- if None attempt to query eos
                  atomicFrac = None,   # Atomic fraction -- if None assume monospecies
                  version = "1.0",     # Optional version string (header)
                  date = None,         # If None, current date is used (header)
                  format = "ASCII"):   # Currently only ASCII allowed

    import datetime
    import numpy as np
    from Spheral1d import CGS, FluidNodeList, ScalarField, makeFluidNodeList

    format = format.lower()
    assert format in ("ascii",)

    # The two files we're going to write
    metadatafile = basename + ".ascii"
    datafile = basename + ".data"

    # Figure our our units conversion (we write LEOS data in CGS units)
    units = eos.constants
    cgs = CGS()
    lconv = cgs.unitLengthMeters/units.unitLengthMeters
    mconv = cgs.unitMassKg/units.unitMassKg
    tconv = cgs.unitTimeSec/units.unitTimeSec
    rhoConv = mconv/(lconv**3)
    Pconv = mconv/(lconv*tconv**2)
    Tconv = 1.0  # We all work in K
    Econv = (lconv/tconv)**2
    CVconv = Econv/Tconv
    velConv = lconv/tconv
    Sconv = Econv/(mconv*Tconv)

    # Fill in any default values
    if rho0 is None:
        assert hasattr(eos, "referenceDensity"), "ERROR: eos does not provide reference density"
        rho0 = eos.referenceDensity
    if T0 is None:
        assert hasattr(eos, "referenceTemperature"), "ERROR: eos does not provide reference temperature"
        T0 = eos.referenceTemperature
    if atomicWeight is None:
        assert hasattr(eos, "atomicWeight"), "ERROR: eos does not provide atomic weight"
        atomicWeight = eos.atomicWeight
    if atomicNumber is None:
        assert hasattr(eos, "atomicNumber"), "ERROR: eos does not provide atomic number"
        atomicNumber = eos.atomicNumber
    if atomicFrac is None:
        assert isinstance(atomicNumber, float) and isinstance(atomicWeight, float), "ERROR: atomicFrac must be the same dimension as atomicWeight & atomicNumber"
        atomicFrac = 1.0
    if date is None:
        date = str(datetime.date.today())

    # Make fields of all the data we're going to write
    rhovals = np.linspace(rhoMin, rhoMax, nrho)
    Tvals = np.linspace(Tmin, Tmax, nT)
    nodes = makeFluidNodeList("stuff", eos,
                              numInternal = nrho*nT)
    rho = ScalarField("rho", nodes)
    T = ScalarField("T", nodes)
    eps = ScalarField("eps", nodes)
    P = ScalarField("P", nodes)
    cs = ScalarField("Cs", nodes)
    S = ScalarField("S", nodes)
    i = 0
    for Ti in Tvals:
        for rhoi in rhovals:
            T[i] = Ti
            rho[i] = rhoi
            i += 1
    eos.setSpecificThermalEnergy(eps, rho, T)
    eos.setPressure(P, rho, eps)
    eos.setSoundSpeed(cs, rho, eps)
    eos.setEntropy(S, rho, eps)

    # Convert Fields to LEOS (CGS) units
    rho /= rhoConv
    eps /= Econv
    P /= Pconv
    cs /= velConv
    S /= Econv

    # Write the meta data file
    with open(metadatafile, "w") as f:
        f.write(f"""data file:
    name: {basename}
    version: {version}
    date: {date}

materials:
    material_name:   {name}
    eosnum:          {eosNumber}
    comp_z:          {atomicNumber}
    comp_a:          {atomicWeight}
    comp_frac:       {atomicFrac}
    rho0:            {rho0}
    t0:              {T0}

    functions: [
""")
        for funcName, unitsstr in [("Et", "erg/g"),
                                   ("Pt", "erg/cc"),
                                   ("Cs", "cm/s"),
                                   ("St", "erg/g-K")]:
            f.write(f"""            {{
            name: {funcName},
            numrho: {nrho},
            numtemp: {nT},
            units_rho: g/cc,
            units_temp: K,
            units_table: {unitsstr},
            }},""")

        f.write(f"""
           ]

    layout:
        type: NCOL
        columns: ["rho", "temp", "Et", "Pt", "Cs", "St"]
        header: True
        headerSize: 1
        table: {datafile}
""")
        
    # Write the tabulated data
    with open(datafile, "w") as f:
        f.write(("#" + 6*" {:16}" + "\n").format("rho", "temp", "Et", "Pt", "Cs", "St"))
        for i in range(nodes.numInternalNodes):
            f.write((6*"{:.10e} " + "\n").format(rho[i], T[i], eps[i], P[i], cs[i], S[i]))

    return
