import json


def get_materials_markers(materials: dict):
    """
    find all materials marker

    model_mat: dict of the specific model wanted
    model_mat: dict of  all the materials
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
    """
    find part of dict from keys, return None if it doesn't exist

    data: dict to explore
    keys: keys that form the path in the dict (data[key1][key2][key3]....)
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
    dict = {}

    exportcfpdes = json_get(PostProcess, "cfpdes", "Exports")
    exportsolid = json_get(PostProcess, "solid", "Exports")
    exportthelec = json_get(PostProcess, "thermo-electric", "Exports")
    exportmaxwell = json_get(PostProcess, "maxwell", "Exports")

    exports = [exportcfpdes, exportsolid, exportthelec, exportmaxwell]
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
                            dict[field] = {"Type": dictType[field], "Exclude": []}

                elif f not in ["pid", "tresca", "material-properties"]:
                    field = (
                        f.replace("heat.", "")
                        .replace("elastic.", "")
                        .replace("magnetic.", "")
                        .replace("electric.", "")
                    )
                    dict[field] = {"Type": dictType[field], "Exclude": []}

        if "expr" in export:
            for f in export["expr"]:
                exclude = []
                include = json_get(export["expr"][f], "markers")
                if include:
                    exclude = list(set(allmarkers) - set(include))
                dict[f] = {"Type": dictType[f], "Exclude": exclude}
    return dict


def returnExportFields(jsonmodel: str, basedir: str):
    with open(jsonmodel, "r") as jsonfile:
        data = json.load(jsonfile)

    with open("./python_hifimagnetParaview/FeelppType.json", "r") as jsonfile:
        feel_dictType = json.load(jsonfile)
    PostProcess = json_get(data, "PostProcess")
    if PostProcess:
        data = jsonFeel_to_fieldDict(data, PostProcess, feel_dictType)

    with open(f"{basedir}/FieldType.json", "w") as fp:
        json.dump(data, fp, indent=4)

    return data
