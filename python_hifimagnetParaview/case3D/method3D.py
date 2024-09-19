import os
import json
import copy


def dictTypeUnits(ureg, distance_unit: str):
    """create dict of units per Type for 3D

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
        "MagneticField_x": {
            "Symbol": "Bx",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "MagneticField_y": {
            "Symbol": "By",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "MagneticField_z": {
            "Symbol": "Bz",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "MagneticField_ur": {
            "Symbol": "Br",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "MagneticField_ut": {
            "Symbol": "Bt",
            "mSymbol": r"$B_{\theta}$",
            "Units": [ureg.tesla, ureg.tesla],
            "Exclude": [],
        },
        "MagneticFieldnorm": {
            "Symbol": "B",
            "mSymbol": r"$\| B \|$",
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
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Displacement_x": {
            "Symbol": "ux",
            "mSymbol": r"$u_x$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Displacement_y": {
            "Symbol": "uy",
            "mSymbol": r"$u_y$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Displacement_z": {
            "Symbol": "uz",
            "mSymbol": r"$u_z$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "Displacement_ur": {
            "Symbol": "ur",
            "mSymbol": r"$u_r$",
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
        "ForceLaplace": {
            "Symbol": "F",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace_x": {
            "Symbol": "Fx",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace_y": {
            "Symbol": "Fy",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace_z": {
            "Symbol": "Fz",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace_ur": {
            "Symbol": "Fr",
            "mSymbol": r"$F_r$",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplace_ut": {
            "Symbol": "Ft",
            "mSymbol": r"$F_{\theta}$",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
        },
        "ForceLaplacenorm": {
            "Symbol": "F",
            "mSymbol": r"$\| F \|$",
            "Units": [
                ureg.newton / ureg.meter**3,
                ureg.newton / ureg.Unit(distance_unit) ** 3,
            ],
            "Exclude": ["Air", "Isolant"],
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
        "Strain_02": {
            "Symbol": "strain_02",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{02}$",
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
        "Strain_12": {
            "Symbol": "strain_12",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{12}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Strain_20": {
            "Symbol": "strain_20",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{20}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Strain_21": {
            "Symbol": "strain_21",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{21}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
            ],
            "Exclude": ["Air"],
        },
        "Strain_22": {
            "Symbol": "strain_22",
            "mSymbol": r"$\bar{\bar{\epsilon}}_{22}$",
            "Units": [
                ureg.dimensionless,
                ureg.dimensionless,
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
        "MagneticPotential_x": {
            "Symbol": "Ax",
            "Units": [
                ureg.ampere / ureg.meter,
                ureg.ampere / ureg.Unit(distance_unit),
            ],
            "Exclude": [],
        },
        "MagneticPotential_y": {
            "Symbol": "Ay",
            "Units": [
                ureg.ampere / ureg.meter,
                ureg.ampere / ureg.Unit(distance_unit),
            ],
            "Exclude": [],
        },
        "MagneticPotential_z": {
            "Symbol": "Az",
            "Units": [
                ureg.ampere / ureg.meter,
                ureg.ampere / ureg.Unit(distance_unit),
            ],
            "Exclude": [],
        },
        "MagneticPotential_ur": {
            "Symbol": "Ar",
            "mSymbol": r"$A_r$",
            "Units": [
                ureg.ampere / ureg.meter,
                ureg.ampere / ureg.Unit(distance_unit),
            ],
            "Exclude": [],
        },
        "MagneticPotential_ut": {
            "Symbol": "At",
            "mSymbol": r"$A_{\theta}$",
            "Units": [
                ureg.ampere / ureg.meter,
                ureg.ampere / ureg.Unit(distance_unit),
            ],
            "Exclude": [],
        },
        "MagneticPotentialnorm": {
            "Symbol": "A",
            "mSymbol": r"$\| A \|$",
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
):
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

    if Type in [
        "Displacement",
        "ElectricField",
        "CurrentDensity",
        "MagneticField",
        "ForceLaplace",
        "MagneticPotential",
    ]:
        for suffix in ["", "_x", "_y", "_z", "_ur", "_ut", "norm"]:
            fieldunits[f"{name}{suffix}"] = TypeUnits[f"{Type}{suffix}"]
            if Exclude:
                fieldunits[f"{name}{suffix}"]["Exclude"] = Exclude
    elif Type in ["Stress"]:
        for suffix in range(3):
            fieldunits[f"{name}_{suffix}"] = TypeUnits[f"{Type}_{suffix}"]
            if Exclude:
                fieldunits[f"{name}_{suffix}"]["Exclude"] = Exclude
    elif Type in ["Strain"]:
        for i in range(3):
            for j in range(3):
                fieldunits[f"{name}_{i}{j}"] = TypeUnits[f"{Type}_{i}{j}"]
                if Exclude:
                    fieldunits[f"{name}_{i}{j}"]["Exclude"] = Exclude
    else:
        fieldunits[name] = TypeUnits[Type]
        if Exclude:
            fieldunits[name]["Exclude"] = Exclude

    return fieldunits


def create_dicts_fromjson(field_dict: dict, ureg, distance_unit: str, basedir: str):
    """create fieldunits dict for 3D from json dict

    Args:
        field_dict (dict): dictionnary of exported fields
        ureg: pint unit registry
        distance_unit (str): unit of distance
        basedir (str): result directory

    Returns:
        fieldunits, ignored_keys
    """

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
        "elasticity.Lame1",
        "elasticity.Lame2",
        "elasticity.PoissonCoefficient",
        "elasticity.YoungModulus",
        "Volume",
        "r",
        "Cos",
        "Sin",
        "coord",
        "Current",
        "B0",
        "Bbg",
        "thermo_electric.pid",
        "pid",
        "cfpdes.pid",
    ]

    with open(f"{basedir}/ignored_keys.json", "w") as fp:
        json.dump({"ignored_keys": ignored_keys}, fp, indent=4)

    return fieldunits, ignored_keys


def create_dicts(ureg, distance_unit: str, basedir: str):
    """create fieldunits dict for 3D from nothing (old version)

    Args:
        ureg: pint unit registry
        distance_unit (str): unit fo distance
        basedir (str): result directory

    Returns:
        fieldunits, ignored_keys
    """

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
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacement_x": {
            "Symbol": "ux",
            "mSymbol": r"$u_x$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacement_y": {
            "Symbol": "uy",
            "mSymbol": r"$u_y$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacement_z": {
            "Symbol": "uz",
            "mSymbol": r"$u_z$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacement_ur": {
            "Symbol": "ur",
            "mSymbol": r"$u_r$",
            "Units": [
                ureg.meter,
                ureg.Unit(distance_unit),
            ],
            "Exclude": ["Air"],
        },
        "displacement_ut": {
            "Symbol": "ut",
            "mSymbol": r"$u_\theta$",
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
