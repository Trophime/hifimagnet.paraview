import os
import json
import copy


from pint import Quantity

# Ignore warning for pint
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])


def create_dicts(ureg, distance_unit: str, basedir: str):
    # use r"$\theta$" for displaying units in mathplotlib
    fieldunits = {
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
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacement_ur": {
            "Symbol": "ur",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacement_ut": {
            "Symbol": "ut",
            "mSymbol": r"$u_{\theta}$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacementnorm": {
            "Symbol": "u",
            "mSymbol": r"$\| u \|$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
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
    ]

    with open(f"{basedir}/ignored_keys.json", "w") as fp:
        json.dump({"ignored_keys": ignored_keys}, fp, indent=4)

    return fieldunits, ignored_keys
