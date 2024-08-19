import json


def get_materials_markers(materials: dict) -> list[str]:
    """find all materials markers

    Args:
        materials (dict): dict of all the materials

    Returns:
        list[str]: list all the markers in materials
    """
    markers_list = []
    for mat in materials.keys():
        markers = json_get(materials, mat, "markers")
        if markers:
            if isinstance(markers, str):
                markers = [markers]
        else:
            markers = [mat]
        markers_list.extend(markers)

    return markers_list


def json_get(data: dict, *keys):
    """find part of dict from keys, return None if it doesn't exist

    Args:
        data (dict): dict to explore
        keys (str): keys that form the path in the dict (data[key1][key2][key3]....)

    Returns:
        data[key1][key2][key3]... or None if doesn't exist
    """
    current_data = data
    if isinstance(current_data, dict):
        for key in keys:
            current_data = current_data.get(key)
            if not current_data:
                break
    else:
        return None

    return current_data


def jsonFeel_to_fieldDict(data: dict, PostProcess: dict, dictType: dict) -> dict:
    """From Feelpp json config file, create a fieldType dict

    Args:
        data (dict): dict of all the Feelpp json config file
        PostProcess (dict): dict post-processing part
        dictType (dict): dict translating name of field in feelpp json to type

    Returns:
        dict: fieldType dict
    """

    dict = {}

    exportcfpdes = json_get(PostProcess, "cfpdes", "Exports")
    exportsolid = json_get(PostProcess, "solid", "Exports")
    exportthelec = json_get(PostProcess, "thermo-electric", "Exports")
    exportCG = json_get(PostProcess, "thermoelectric", "Exports")
    exportmaxwell = json_get(PostProcess, "maxwell", "Exports")

    exports = [exportcfpdes, exportsolid, exportthelec, exportmaxwell, exportCG]
    exports = [x for x in exports if x is not None]

    allmarkers = get_materials_markers(data["Materials"])

    for export in exports:
        if "fields" in export:
            if isinstance(export["fields"], str):
                export["fields"] = [export["fields"]]
            for f in export["fields"]:
                if f == "all":
                    for m in data["Models"]:
                        field = json_get(
                            data, "Models", m, "common", "setup", "unknown", "name"
                        )
                        if field:
                            dict[field.replace("-", "_")] = {
                                "Type": dictType[field],
                                "Exclude": [],
                            }

                elif f not in ["pid", "tresca", "material-properties"]:
                    field = (
                        f.replace("heat.", "")
                        .replace("elastic.", "")
                        .replace("magnetic.", "")
                        .replace("electric.", "")
                    )
                    dict[field.replace("-", "_")] = {
                        "Type": dictType[field],
                        "Exclude": [],
                    }

        if "expr" in export:
            for f in export["expr"]:
                exclude = []
                include = json_get(export["expr"][f], "markers")
                if include:
                    if isinstance(include, str):
                        include = [include]
                    exclude = list(set(allmarkers) - set(include))
                dict[f] = {"Type": dictType[f], "Exclude": exclude}
    return dict


def returnExportFields(jsonmodel: str, basedir: str) -> dict:
    """create FieldType.json, with all exported fields and their type

    Args:
        jsonmodel (str): json file from feelpp or with fields and type
        basedir (str): result directory

    Returns:
        dict: fieldType dict
    """
    with open(jsonmodel, "r") as jsonfile:
        data = json.load(jsonfile)

    PostProcess = json_get(data, "PostProcess")
    if PostProcess:
        with open("./python_hifimagnetParaview/FeelppType.json", "r") as jsonfile:
            feel_dictType = json.load(jsonfile)
        data = jsonFeel_to_fieldDict(data, PostProcess, feel_dictType)

    with open(f"{basedir}/FieldType.json", "w") as fp:
        json.dump(data, fp, indent=4)

    return data
