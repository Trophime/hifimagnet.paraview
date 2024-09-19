import os
import json
import copy


def dictTypeUnits(ureg, distance_unit: str) -> dict:
    """create dict of units per Type for 2D

    Args:
        ureg: pint unit registry
        distance_unit (str): unit of distance

    Returns:
        dict: dict of unit per type
    """

    TypeUnits = {
        "ThermalConductivity": {
            "Symbol": "k",
            "Units": [
                ureg.watt / ureg.meter / ureg.kelvin,
                ureg.watt / ureg.Unit(distance_unit) / ureg.kelvin,
            ],
            "Exclude": ["Air"],
        },
        "ElectricConductivity": {
            "Symbol": "Sigma",
            "mSymobol": r"$\sigma$",
            "Units": [
                ureg.siemens / ureg.meter,
                ureg.siemens / ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "PoissonCoef": {
            "Symbol": "Poisson",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "YoungModulus": {
            "Symbol": "YoungModulus",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "BackgroundMagneticField": {
            "Symbol": "B_Bg",
            "mSymbol": r"$B_{bg}$",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "MagneticField": {
            "Symbol": "B",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "Temperature": {
            "Symbol": "T",
            "Units": [ureg.degK, ureg.degC],
            "Exclude": ["Air"],
        },
        "TemperatureCoefficient": {
            "Symbol": "alpha",
            "mSymbol": r"$\alpha$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Stress_T": {
            "Symbol": "stress_T",
            "mSymbol": r"$\bar{\bar{\sigma}}_{T}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "ElectricField": {
            "Symbol": "E",
            "Units": [ureg.volt / ureg.meter, ureg.volt / ureg.Unit(distance_unit)],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricField_ur": {
            "Symbol": "Er",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricField_ut": {
            "Symbol": "Et",
            "mSymbol": r"$E_{\theta}$",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricFieldnorm": {
            "Symbol": "E",
            "mSymbol": r"$\| E \|$",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricPotential": {
            "Symbol": "V",
            "Units": [ureg.volt, ureg.volt],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity": {
            "Symbol": "J",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity_ur": {
            "Symbol": "Jr",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity_ut": {
            "Symbol": "Jt",
            "mSymbol": r"$J_{\theta}$",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensitynorm": {
            "Symbol": "J",
            "mSymbol": r"$\| J \|$",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace": {
            "Symbol": "F",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace_ur": {
            "Symbol": "Fr",
            "Units": [
                ureg.newton / ureg.meter**2,
                ureg.newton / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace_ut": {
            "Symbol": "Ft",
            "mSymbol": r"$F_{\theta}$",
            "Units": [
                ureg.newton / ureg.meter**2,
                ureg.newton / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplacenorm": {
            "Symbol": "F",
            "mSymbol": r"$\| F \|$",
            "Units": [
                ureg.newton / ureg.meter**2,
                ureg.newton / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "Displacement": {
            "Symbol": "u",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Displacement_ur": {
            "Symbol": "ur",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Displacement_ut": {
            "Symbol": "ut",
            "mSymbol": r"$u_{\theta}$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Displacementnorm": {
            "Symbol": "u",
            "mSymbol": r"$\| u \|$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Strain_00": {
            "Symbol": "strain_00",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{00}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Strain_01": {
            "Symbol": "strain_01",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{01}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Strain_10": {
            "Symbol": "strain_10",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{10}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Strain_11": {
            "Symbol": "strain_11",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{11}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Stress_00": {
            "Symbol": "stress_00",
            "mSymbol": r"$\bar{\bar{\sigma}}_{00}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Stress_01": {
            "Symbol": "stress_01",
            "mSymbol": r"$\bar{\bar{\sigma}}_{01}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Stress_10": {
            "Symbol": "stress_10",
            "mSymbol": r"$\bar{\bar{\sigma}}_{10}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Stress_11": {
            "Symbol": "stress_11",
            "mSymbol": r"$\bar{\bar{\sigma}}_{11}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "VonMises": {
            "Symbol": "VonMises",
            "mSymbol": r"$\bar{\bar{\sigma}}_{VonMises}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Q": {
            "Symbol": "Qth",
            "Units": [
                ureg.watt / ureg.meter**3,
                ureg.megawatt / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricPotential": {
            "Symbol": "V",
            "Units": [ureg.volt, ureg.volt],
            "Exclude": ["Air", "Isolant"],
        },
        "MagneticPotential": {
            "Symbol": "A",
            "Units": [
                ureg.ampere / ureg.meter,
                ureg.ampere / ureg.Unit(distance_unit),
            ],
            "Exclude": [],
        },
        "HoopStrain": {
            "Symbol": "HoopStrain",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{Hoop}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "HoopStress": {
            "Symbol": "HoopStress",
            "mSymbol": r"$\bar{\bar{\sigma}}_{Hoop}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
    }

    return TypeUnits


def addFieldToFieldunits(
    fieldunits: dict, name: str, Type: str, Exclude: list[str], TypeUnits: dict
) -> dict:
    """add field to fieldunits dict with units corresponding to its type

    Args:
        fieldunits (dict): dict of field units
        name (str): name of field
        Type (str): type of field
        Exclude (list[str]): list of excluded marker of field
        TypeUnits (dict): dict of Type units

    Returns:
        dict: updated fieldunits
    """

    if Type in ["Displacement", "ElectricField", "CurrentDensity", "ForceLaplace"]:
        for suffix in ["", "_ur", "_ut", "norm"]:
            fieldunits[f"{name}{suffix}"] = TypeUnits[f"{Type}{suffix}"]
            if Exclude:
                fieldunits[f"{name}{suffix}"]["Exclude"] = Exclude
    elif Type in ["Stress", "Strain"]:
        for i in range(2):
            for j in range(2):
                fieldunits[f"{name}_{i}{j}"] = TypeUnits[f"{Type}_{i}{j}"]
                if Exclude:
                    fieldunits[f"{name}_{i}{j}"]["Exclude"] = Exclude
    else:
        fieldunits[name] = TypeUnits[Type]
        if Exclude:
            fieldunits[name]["Exclude"] = Exclude

    return fieldunits


def create_dicts_fromjson(field_dict: dict, ureg, distance_unit: str, basedir: str):
    """create fieldunits dict for 2D from json dict

    Args:
        field_dict (dict): dictionnary of exported fields
        ureg: pint unit registry
        distance_unit (str): unit of distance
        basedir (str): result directory

    Returns:
        fieldunits, ignored_keys
    """

    # use r"$\theta$" for displaying units in mathplotlib
    fieldunits = {
        "coord": {
            "Symbol": "r",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
        },
        "Length": {"Symbol": "", "Units": [ureg.meter, ureg.Unit(distance_unit)]},
        "Area": {
            "Symbol": "S",
            "Units": [
                ureg.meter**2,
                ureg.Unit(distance_unit) ** 2,
            ],
        },
        "Volume": {
            "Symbol": "V",
            "Units": [
                ureg.meter**3,
                ureg.Unit(distance_unit) ** 3,
            ],
        },
        "Current": {
            "Val": None,
            "Units": [
                ureg.amperes,
                ureg.amperes,
            ],
        },
        "B0": {
            "Val": None,
            "Units": [
                ureg.tesla,
                ureg.tesla,
            ],
        },
        "Bbg": {
            "Val": None,
            "Units": [
                ureg.tesla,
                ureg.tesla,
            ],
        },
    }

    TypeUnits = dictTypeUnits(ureg, distance_unit)

    for f in field_dict:
        fieldunits = addFieldToFieldunits(
            fieldunits, f, field_dict[f]["Type"], field_dict[f]["Exclude"], TypeUnits
        )

    tmp = copy.deepcopy(fieldunits)
    for field, values in tmp.items():
        del values["Units"]
    with open(f"{basedir}/fieldunits.json", "w") as fp:
        json.dump(tmp, fp, indent=4)
    del tmp

    ignored_keys = [
        "Area",
        "Volume",
        "r",
        "Cos",
        "Sin",
        "cfpdes.pid",
        "thermo_electric.pid",
        "pid",
        "coord",
        "Current",
        "B0",
        "Bbg",
    ]

    with open(f"{basedir}/ignored_keys.json", "w") as fp:
        json.dump({"ignored_keys": ignored_keys}, fp, indent=4)

    return fieldunits, ignored_keys


def create_dicts(ureg, distance_unit: str, basedir: str):
    """create fieldunits dict for 2D from nothing (old version)

    Args:
        ureg: pint unit registry
        distance_unit (str): unit fo distance
        basedir (str): result directory

    Returns:
        fieldunits, ignored_keys
    """
    # use r"$\theta$" for displaying units in mathplotlib
    fieldunits = {
        "coord": {
            "Symbol": "r",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
        },
        "VolumicMass": {
            "Symbol": "rho",
            "mSymbol": r"$\rho$",
            "Units": [
                ureg.kilogram / ureg.meter**3,
                ureg.kilogram / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air"],
        },
        "k": {
            "Symbol": "k",
            "Units": [
                ureg.watt / ureg.meter / ureg.kelvin,
                ureg.watt / ureg.Unit(distance_unit) / ureg.kelvin,
            ],
            "Exclude": ["Air"],
        },
        "sigma": {
            "Symbol": "Sigma",
            "mSymobol": r"$\sigma$",
            "Units": [
                ureg.siemens / ureg.meter,
                ureg.siemens / ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "nu": {
            "Symbol": "Poisson",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "EE": {
            "Symbol": "YoungModulus",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Length": {"Symbol": "", "Units": [ureg.meter, ureg.Unit(distance_unit)]},
        "Area": {
            "Symbol": "S",
            "Units": [
                ureg.meter**2,
                ureg.Unit(distance_unit) ** 2,
            ],
        },
        "Volume": {
            "Symbol": "V",
            "Units": [
                ureg.meter**3,
                ureg.Unit(distance_unit) ** 3,
            ],
        },
        "mu0": {
            "Symbol": "mu",
            "Units": [ureg.henry / ureg.meter, ureg.henry / ureg.Unit(distance_unit)],
        },
        "h": {
            "Symbol": "hw",
            "Units": [
                ureg.watt / ureg.meter**2 / ureg.kelvin,
                ureg.watt / ureg.Unit(distance_unit) ** 2 / ureg.kelvin,
            ],
        },
        "Flow": {
            "Symbol": "Q",
            "Units": [
                ureg.liter / ureg.second,
                ureg.Unit(distance_unit)
                * ureg.Unit(distance_unit)
                * ureg.Unit(distance_unit)
                / ureg.second,
            ],
        },
        "Current": {"Symbol": "I", "Units": [ureg.ampere, ureg.ampere]},
        "Power": {"Symbol": "W", "Units": [ureg.watt, ureg.megawatt]},
        "Qth": {
            "Symbol": "Qth",
            "Units": [
                ureg.watt / ureg.meter**3,
                ureg.megawatt / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "temperature": {
            "Symbol": "T",
            "Units": [ureg.degK, ureg.degC],
            "Exclude": ["Air"],
        },
        "U": {
            "Symbol": "V",
            "Units": [ureg.volt, ureg.volt],
            "Exclude": ["Air", "Isolant"],
        },
        "E": {
            "Symbol": "E",
            "Units": [ureg.volt / ureg.meter, ureg.volt / ureg.Unit(distance_unit)],
            "Exclude": ["Air", "Isolant"],
        },
        "E_ur": {
            "Symbol": "Er",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "E_ut": {
            "Symbol": "Et",
            "mSymbol": r"$E_{\theta}$",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "Enorm": {
            "Symbol": "E",
            "mSymbol": r"$\| E \|$",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "J": {
            "Symbol": "J",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "J_ur": {
            "Symbol": "Jr",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "J_ut": {
            "Symbol": "Jt",
            "mSymbol": r"$J_{\theta}$",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "Jnorm": {
            "Symbol": "J",
            "mSymbol": r"$\| J \|$",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "F_laplace": {
            "Symbol": "F",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "F_laplace_ur": {
            "Symbol": "Fr",
            "Units": [
                ureg.newton / ureg.meter**2,
                ureg.newton / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "F_laplace_ut": {
            "Symbol": "Ft",
            "mSymbol": r"$F_{\theta}$",
            "Units": [
                ureg.newton / ureg.meter**2,
                ureg.newton / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "F_laplacenorm": {
            "Symbol": "F",
            "mSymbol": r"$\| F \|$",
            "Units": [
                ureg.newton / ureg.meter**2,
                ureg.newton / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "Bz": {
            "Symbol": "B",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "displacement": {
            "Symbol": "u",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacement_ur": {
            "Symbol": "ur",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacement_ut": {
            "Symbol": "ut",
            "mSymbol": r"$u_{\theta}$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacementnorm": {
            "Symbol": "u",
            "mSymbol": r"$\| u \|$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Lame1": {
            "Symbol": "Lame1",
            "Units": [
                ureg.pascal / ureg.meter,
                ureg.megapascal / ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Lame2": {
            "Symbol": "Lame2",
            "Units": [
                ureg.pascal / ureg.meter,
                ureg.megapascal / ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        # "PoissonCoefficient": {}, # no units
        "strain_00": {
            "Symbol": "strain_00",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{00}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "strain_01": {
            "Symbol": "strain_01",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{01}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "strain_10": {
            "Symbol": "strain_10",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{10}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "strain_11": {
            "Symbol": "strain_11",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{11}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "stress_00": {
            "Symbol": "stress_00",
            "mSymbol": r"$\bar{\bar{\sigma}}_{00}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "stress_T": {
            "Symbol": "stress_T",
            "mSymbol": r"$\bar{\bar{\sigma}}_{T}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "stress_01": {
            "Symbol": "stress_01",
            "mSymbol": r"$\bar{\bar{\sigma}}_{01}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "stress_10": {
            "Symbol": "stress_10",
            "mSymbol": r"$\bar{\bar{\sigma}}_{10}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "stress_11": {
            "Symbol": "stress_11",
            "mSymbol": r"$\bar{\bar{\sigma}}_{11}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Vonmises": {
            "Symbol": "VonMises",
            "mSymbol": r"$\bar{\bar{\sigma}}_{VonMises}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
    }

    tmp = copy.deepcopy(fieldunits)
    for field, values in tmp.items():
        del values["Units"]
    with open(f"{basedir}/fieldunits.json", "w") as fp:
        json.dump(tmp, fp, indent=4)
    del tmp

    ignored_keys = [
        "cfpdes.expr.nu",
        "cfpdes.expr.EE",
        "Area",
        "Volume",
        "r",
        "Cos",
        "Sin",
        "cfpdes.pid",
        "coord",
    ]

    with open(f"{basedir}/ignored_keys.json", "w") as fp:
        json.dump({"ignored_keys": ignored_keys}, fp, indent=4)

    return fieldunits, ignored_keys
