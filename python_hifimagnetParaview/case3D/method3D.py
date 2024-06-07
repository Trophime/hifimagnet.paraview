import os
import json
import copy
from typing import List


def dictTypeUnits(ureg, distance_unit: str):
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
        "ElectricField_x": {
            "Symbol": "Ex",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricField_y": {
            "Symbol": "Ey",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricField_z": {
            "Symbol": "Ez",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricField_r": {
            "Symbol": "Er",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricField_t": {
            "Symbol": "Et",
            "mSymbol": r"$E_{\theta}$",
            "Units": [
                ureg.volt / ureg.meter**2,
                ureg.volt / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ElectricField_norm": {
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
        "CurrentDensity_x": {
            "Symbol": "Jx",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity_y": {
            "Symbol": "Jy",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity_z": {
            "Symbol": "Jz",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity_r": {
            "Symbol": "Jr",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity_t": {
            "Symbol": "Jt",
            "mSymbol": r"$J_{\theta}$",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "CurrentDensity_norm": {
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
        "ForceLaplace_ext": {
            "Symbol": "Fext",
            "Units": [
                ureg.newton / ureg.meter**2,
                ureg.newton / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "Displacement": {
            "Symbol": "u",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "Displacement_x": {
            "Symbol": "ux",
            "mSymbol": r"$u_x$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "Displacement_y": {
            "Symbol": "uy",
            "mSymbol": r"$u_y$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "Displacement_z": {
            "Symbol": "uz",
            "mSymbol": r"$u_z$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "Displacement_r": {
            "Symbol": "ur",
            "mSymbol": r"$u_r$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "Displacement_t": {
            "Symbol": "ut",
            "mSymbol": r"$u_{\theta}$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "Displacement_norm": {
            "Symbol": "u",
            "mSymbol": r"$\| u \|$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "Stress_0": {
            "Symbol": "principal_stress_0",
            "mSymbol": r"$\bar{\bar{\sigma}}_0$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Stress_1": {
            "Symbol": "principal_stress_1",
            "mSymbol": r"$\bar{\bar{\sigma}}_1$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Stress_2": {
            "Symbol": "principal_stress_2",
            "mSymbol": r"$\bar{\bar{\sigma}}_2$",
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
    fieldunits: dict, name: str, Type: str, Exclude: List[str], TypeUnits: dict
):
    if Type in ["Displacement", "ElectricField", "CurrentDensity"]:
        fieldunits[name] = TypeUnits[Type]
        fieldunits[f"{name}_x"] = TypeUnits[f"{Type}_x"]
        fieldunits[f"{name}_y"] = TypeUnits[f"{Type}_y"]
        fieldunits[f"{name}_z"] = TypeUnits[f"{Type}_z"]
        fieldunits[f"{name}_ur"] = TypeUnits[f"{Type}_r"]
        fieldunits[f"{name}_ut"] = TypeUnits[f"{Type}_t"]
        fieldunits[f"{name}norm"] = TypeUnits[f"{Type}_norm"]
        if Exclude:
            fieldunits[name]["Exclude"] = Exclude
            fieldunits[f"{name}_x"]["Exclude"]  = Exclude
            fieldunits[f"{name}_y"]["Exclude"]  = Exclude
            fieldunits[f"{name}_z"]["Exclude"]  = Exclude
            fieldunits[f"{name}_ur"]["Exclude"] = Exclude
            fieldunits[f"{name}_ut"]["Exclude"] = Exclude
            fieldunits[f"{name}norm"]["Exclude"] = Exclude
    elif Type in ["Stress"]:
        fieldunits[f"{name}_0"] = TypeUnits[f"{Type}_0"]
        fieldunits[f"{name}_1"] = TypeUnits[f"{Type}_1"]
        fieldunits[f"{name}_2"] = TypeUnits[f"{Type}_2"]
        if Exclude:
            fieldunits[f"{name}_0"]["Exclude"] = Exclude
            fieldunits[f"{name}_1"]["Exclude"] = Exclude
            fieldunits[f"{name}_2"]["Exclude"] = Exclude
    else:
        fieldunits[name] = TypeUnits[Type]
        if Exclude:
            fieldunits[name]["Exclude"] = Exclude

    return fieldunits


def create_dicts_fromjson(field_dict: dict, ureg, distance_unit: str, basedir: str):
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
            "Symbol": "A",
            "Units": [
                ureg.meter**2,
                ureg.Unit(distance_unit) ** 2,
            ],
        },
        "Volume": {
            "Symbol": "V",
            "Units": [ureg.meter**3, ureg.Unit(distance_unit) ** 3],
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
        "elasticity.Lame1",
        "elasticity.Lame2",
        "elasticity.PoissonCoefficient",
        "elasticity.YoungModulus",
        "Volume",
        "r",
        "Cos",
        "Sin",
        "coord",
    ]

    with open(f"{basedir}/ignored_keys.json", "w") as fp:
        json.dump({"ignored_keys": ignored_keys}, fp, indent=4)

    return fieldunits, ignored_keys


def create_dicts(ureg, distance_unit: str, basedir: str):
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
        "ThermalConductivity": {
            "Symbol": "k",
            "Units": [
                ureg.watt / ureg.meter / ureg.kelvin,
                ureg.watt / ureg.Unit(distance_unit) / ureg.kelvin,
            ],
            "Exclude": ["Air"],
        },
        "ElectricalConductivity": {
            "Symbol": "Sigma",
            "mSymobol": r"$\sigma$",
            "Units": [
                ureg.siemens / ureg.meter,
                ureg.siemens / ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "YoungModulus": {
            "Symbol": "E",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "Length": {"Symbol": "", "Units": [ureg.meter, ureg.Unit(distance_unit)]},
        "Area": {
            "Symbol": "A",
            "Units": [
                ureg.meter**2,
                ureg.Unit(distance_unit) ** 2,
            ],
        },
        "Volume": {
            "Symbol": "V",
            "Units": [ureg.meter**3, ureg.Unit(distance_unit) ** 3],
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
        "Power": {"Symbol": "W", "Units": [ureg.watt, ureg.watt]},
        "temperature": {
            "Symbol": "T",
            "Units": [ureg.degK, ureg.degC],
            "Exclude": ["Air"],
        },
        "electric_potential": {
            "Symbol": "V",
            "Units": [ureg.volt, ureg.volt],
            "Exclude": ["Air", "Isolant"],
        },
        "current_density": {
            "Symbol": "J",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "current_densitynorm": {
            "Symbol": "J",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "current_density_x": {
            "Symbol": "Jx",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "current_density_y": {
            "Symbol": "Jy",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "current_density_z": {
            "Symbol": "Jz",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "current_density_ur": {
            "Symbol": "Jr",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "current_density_ut": {
            "Symbol": "Jt",
            "Units": [
                ureg.ampere / ureg.meter**2,
                ureg.ampere / ureg.Unit(distance_unit) ** 2,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "Fext": {
            "Symbol": "F_ext",
            "mSymbol": r"$F_{ext}$",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "Fh": {
            "Symbol": "F",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "magnetic_potential": {
            "Symbol": "A",
            "Units": [
                ureg.ampere / ureg.meter,
                ureg.ampere / ureg.Unit(distance_unit),
            ],
        },
        "Bg_MagField": {
            "Symbol": "B_Bg",
            "mSymbol": r"$B_{bg}$",
            "Units": [ureg.tesla, ureg.tesla],
        },
        "magnetic_field": {"Symbol": "B", "Units": [ureg.tesla, ureg.tesla]},
        "displacement": {
            "Symbol": "u",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacement_x": {
            "Symbol": "ux",
            "mSymbol": r"$u_x$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacement_y": {
            "Symbol": "uy",
            "mSymbol": r"$u_y$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacement_z": {
            "Symbol": "uz",
            "mSymbol": r"$u_z$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacement_ur": {
            "Symbol": "ur",
            "mSymbol": r"$u_r$",
            "Units": [
                ureg.meter / ureg.second,
                ureg.Unit(distance_unit) / ureg.second,
            ],
            "Exclude": ["Air"],
        },
        "displacement_ut": {
            "Symbol": "ut",
            "mSymbol": r"$u_\theta$",
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
        "princial_stress_0": {
            "Symbol": "principal_stress_0",
            "mSymbol": r"$\bar{\bar{\sigma}}_0$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "princial_stress_1": {
            "Symbol": "principal_stress_1",
            "mSymbol": r"$\bar{\bar{\sigma}}_1$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "princial_stress_2": {
            "Symbol": "principal_stress_2",
            "mSymbol": r"$\bar{\bar{\sigma}}_2$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
        "von_mises_criterions": {
            "Symbol": "VonMises",
            "mSymbol": r"$\bar{\bar{\sigma}}_{VonMises}$",
            "Units": [
                ureg.pascal,
                ureg.megapascal,
            ],
            "Exclude": ["Air"],
        },
    }

    # build fieldunits when creating setup and store as a json
    # load json
    # how to dump and load pint unit??
    # could the units be stored in the python code?
    # how to attach a field to a unit?
    # see: https://pint.readthedocs.io/en/0.10.1/serialization.html

    tmp = copy.deepcopy(fieldunits)
    for field, values in tmp.items():
        del values["Units"]
    with open(f"{basedir}/fieldunits.json", "w") as fp:
        json.dump(tmp, fp, indent=4)
    del tmp

    ignored_keys = [
        "elasticity.Lame1",
        "elasticity.Lame2",
        "elasticity.PoissonCoefficient",
        "elasticity.YoungModulus",
        "Volume",
        "r",
        "Cos",
        "Sin",
        "coord",
    ]

    with open(f"{basedir}/ignored_keys.json", "w") as fp:
        json.dump({"ignored_keys": ignored_keys}, fp, indent=4)

    return fieldunits, ignored_keys
