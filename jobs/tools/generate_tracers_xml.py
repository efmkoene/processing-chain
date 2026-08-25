from pathlib import Path
import re
import xml.etree.ElementTree as ET
import xml.dom.minidom
import numpy as np


def generate_tracers_xml(data,
                         nens=-1,
                         n_bg_ens=-1,
                         restart=False,
                         runthrough=False,
                         propagate_bg=False):
    """
    Generate an XML representation for chemtracers.

    Args:
        data (dict):
            A dictionary containing details for chemtracers. Example structure:
            {
                "TRCO2_A": {
                    "oem_cat": "A-CO2, ...",
                    "oem_vp": "GNFR_A, ...",
                    "oem_tp": "GNFR_A-CO2, ..."
                },
                "TRCO2_BG": {
                    "init_name": "CO2"
                },
                "CO2_RA": {},
                "CO2_GPP": {},
                "TRCO2_A-XXX": {"bg": "TRCO2_BG", "ra": "CO2_RA", "gpp": "CO2_GPP"}
            }

    Returns:
        str: The prettyfied XML string.
    """
    tracers = ET.Element("tracers")

    # Iterate over all items in data
    for item_id, item_data in data.items():
        if any(key == "oem_cat" for key in item_data):
            # Make an OEM tracer
            tracer = ET.SubElement(tracers, "chemtracer", id=item_id)
            ET.SubElement(
                tracer, "transport", type="char"
            ).text = "stdaero" if not item_id.startswith("EM_") else "off"
            if item_id.startswith("EM_"):
                ET.SubElement(tracer, "iconv", type="int").text = "0"
                ET.SubElement(tracer, "iturb", type="int").text = "0"
                ET.SubElement(tracer, "output_emis", type="int").text = "1"
            ET.SubElement(tracer, "c_solve", type="char").text = "passive"
            ET.SubElement(tracer, "init_mode", type="int").text = "0"
            ET.SubElement(tracer, "unit", type="char").text = "none"
            ET.SubElement(tracer, "oem_tscale", type="int").text = "2"
            ET.SubElement(tracer, "oem_type", type="char").text = "emis"
            for key, value in item_data.items():
                if key.startswith("oem_"):
                    ET.SubElement(tracer, key, type="char").text = value
            if restart and not item_id.startswith("EM_"):
                ET.SubElement(tracer, "oem_restart", type="char").text = "file"
        if item_id.endswith("BG"):
            # Make a background tracer
            tracer_bg = ET.SubElement(tracers, "chemtracer", id=item_id)
            ET.SubElement(tracer_bg, "transport", type="char").text = "stdaero"
            ET.SubElement(tracer_bg, "c_solve", type="char").text = "passive"
            ET.SubElement(tracer_bg, "init_mode", type="int").text = "1"
            ET.SubElement(tracer_bg, "unit", type="char").text = "none"
            ET.SubElement(tracer_bg, "init_name",
                          type="char").text = item_data["init_name"]
            ET.SubElement(tracer_bg, "oem_type", type="char").text = "bg"
            if restart:
                ET.SubElement(tracer_bg, "oem_restart",
                              type="char").text = "file"
            ET.SubElement(tracer_bg, "latbc", type="char").text = "file"
        if any(key == "oem_ftype" for key in item_data):
            # Make a VPRM tracer
            tracer_ra = ET.SubElement(tracers, "chemtracer", id=item_id)
            ET.SubElement(
                tracer_ra, "transport", type="char"
            ).text = "stdaero" if not item_id.startswith("EM_") else "off"
            if item_id.startswith("EM_"):
                ET.SubElement(tracer_ra, "iconv", type="int").text = "0"
                ET.SubElement(tracer_ra, "iturb", type="int").text = "0"
                ET.SubElement(tracer_ra, "output_emis", type="int").text = "1"
            ET.SubElement(tracer_ra, "c_solve", type="char").text = "passive"
            ET.SubElement(tracer_ra, "init_mode", type="int").text = "0"
            ET.SubElement(tracer_ra, "unit", type="char").text = "none"
            ET.SubElement(tracer_ra, "oem_type", type="char").text = "vprm"
            ET.SubElement(tracer_ra, "oem_ftype",
                          type="char").text = item_data["oem_ftype"]
            if restart and not item_id.startswith("EM_"):
                ET.SubElement(tracer_ra, "oem_restart",
                              type="char").text = "file"
        if not runthrough:
            if item_id.endswith("XXX"):
                # Make a set of ensemble tracers
                for i in np.arange(nens) + 1:
                    tracer_xxx = ET.SubElement(tracers,
                                               "chemtracer",
                                               id=f"{item_id[:-4]}-{i:03}")
                    ET.SubElement(tracer_xxx, "transport",
                                  type="char").text = "stdaero"
                    ET.SubElement(tracer_xxx, "oem_type",
                                  type="char").text = "ens"
                    ET.SubElement(tracer_xxx, "c_solve",
                                  type="char").text = "passive"
                    ET.SubElement(tracer_xxx, "init_mode",
                                  type="int").text = "0"
                    if "bg" in item_data:
                        ET.SubElement(tracer_xxx, "oem_bg_ens",
                                      type="char").text = item_data["bg"]
                    if "ra" in item_data and "gpp" in item_data:
                        ET.SubElement(
                            tracer_xxx, "oem_vprm_bg_ens", type="char"
                        ).text = f"{item_data['ra']}, {item_data['gpp']}"
                    if restart:
                        ET.SubElement(tracer_xxx, "oem_restart",
                                      type="char").text = "file"
                    ET.SubElement(tracer_xxx, "unit",
                                  type="char").text = "none"
                if propagate_bg:
                    tracer_xxx = ET.SubElement(
                        tracers,
                        "chemtracer",
                        id=f"{item_id[:-4]}-{nens+1:03}")
                    ET.SubElement(tracer_xxx, "transport",
                                  type="char").text = "stdaero"
                    ET.SubElement(tracer_xxx, "oem_type",
                                  type="char").text = "ens"
                    ET.SubElement(tracer_xxx, "c_solve",
                                  type="char").text = "passive"
                    ET.SubElement(tracer_xxx, "init_mode",
                                  type="int").text = "0"
                    if "bg" in item_data:
                        ET.SubElement(tracer_xxx, "oem_bg_ens",
                                      type="char").text = item_data["bg"]
                    if restart:
                        ET.SubElement(tracer_xxx, "oem_restart",
                                      type="char").text = "file"
                    ET.SubElement(tracer_xxx, "unit",
                                  type="char").text = "none"
        else:
            if item_id.endswith("XXX"):
                for i in np.arange(n_bg_ens+1) + 1:
                    tracer_bg_xxx = ET.SubElement(tracers,
                                                  "chemtracer",
                                                  id=f"{item_id[:-4]}-{i:03}")
                    ET.SubElement(tracer_bg_xxx, "transport",
                                  type="char").text = "stdaero"
                    ET.SubElement(tracer_bg_xxx, "oem_type",
                                  type="char").text = "ens"
                    ET.SubElement(tracer_bg_xxx, "c_solve",
                                  type="char").text = "passive"
                    ET.SubElement(tracer_bg_xxx, "init_mode",
                                  type="int").text = "0"
                    if "bg" in item_data:
                        ET.SubElement(tracer_bg_xxx, "oem_bg_ens",
                                      type="char").text = item_data["bg"]
                    if i == n_bg_ens + 1:
                        if "ra" in item_data and "gpp" in item_data:
                            ET.SubElement(
                                tracer_bg_xxx, "oem_vprm_bg_ens", type="char"
                            ).text = f"{item_data['ra']}, {item_data['gpp']}"
                    if restart:
                        ET.SubElement(tracer_bg_xxx,
                                      "oem_restart",
                                      type="char").text = "file"
                    ET.SubElement(tracer_bg_xxx, "unit",
                                  type="char").text = "none"
    # Convert to string
    xml_declaration = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE tracers SYSTEM \"tracers.dtd\">\n"
    xml_string = ET.tostring(tracers, encoding="unicode")
    return xml.dom.minidom.parseString(xml_declaration +
                                       xml_string).toprettyxml()


def reduce_tracers_xml(tracers_xml_in: Path,
                       tracers_xml_out: Path,
                       cfg):
    """
    Reduce tracer ensembles in ICON tracer XML
    Keeps:
        - first ensemble member  → renamed to -001
        - if cfg.CTDAS_propagate_bg: last ensemble member → renamed to -002
    Removes all other ensemble members.

    Args:
        tracers_xml_in: Path to input tracers.xml
        tracers_xml_out: Path to write reduced tracers.xml
        cfg: config object with:
            - cfg.tracers: dict of tracer keys, one ending in -XXX
            - cfg.CTDAS_nensembles
            - cfg.CTDAS_propagate_bg
    """

    # ----------------------------------------------------
    # 1. Identify prefix from tracers ending in "-XXX"
    # ----------------------------------------------------
    ens_keys = [k for k in cfg.tracers.keys() if k.endswith("-XXX")]
    if len(ens_keys) != 1:
        raise ValueError("Expected exactly one tracer key ending in -XXX")
    template_key = ens_keys[0]
    prefix = template_key[:-4]    # remove "-XXX"

    # Pattern to match ensemble tracers in XML, e.g. TRCO2_A-001
    pattern = re.compile(rf"^{re.escape(prefix)}-(\d{{3}})$")

    # ----------------------------------------------------
    # 2. Parse XML using ElementTree
    # ----------------------------------------------------
    tree = ET.parse(str(tracers_xml_in))
    root = tree.getroot()

    # ----------------------------------------------------
    # 3. Extract all ensemble tracer nodes
    # ----------------------------------------------------
    ens_nodes = []
    for node in list(root.findall("chemtracer")):   # list() because we may remove
        tid = node.get("id")
        m = pattern.match(tid)
        if m:
            num = int(m.group(1))
            ens_nodes.append((num, node))

    if not ens_nodes:
        raise ValueError("No ensemble tracer entries found in tracers.xml")

    # Sort by numeric suffix
    ens_nodes.sort(key=lambda x: x[0])

    # ----------------------------------------------------
    # 4. Decide which nodes to keep
    # ----------------------------------------------------
    first_num, first_node = ens_nodes[0]
    to_keep = [(first_node, "001")]   # always keep first

    if cfg.CTDAS_propagate_bg:
        last_num, last_node = ens_nodes[-1]
        to_keep.append((last_node, "002"))

    # ----------------------------------------------------
    # 5. Remove all ensemble nodes
    # ----------------------------------------------------
    for _, node in ens_nodes:
        root.remove(node)

    # ----------------------------------------------------
    # 6. Reinsert the kept nodes with renamed IDs
    # ----------------------------------------------------
    for node, suffix in to_keep:
        new_id = f"{prefix}-{suffix}"
        node.set("id", new_id)
        root.append(node)

    # ----------------------------------------------------
    # 7. Pretty-print and write output (minidom)
    # ----------------------------------------------------
    def strip_whitespace_nodes(elem):
        """Recursively remove text nodes that are only whitespace"""
        if elem.text:
            elem.text = elem.text.strip()
        if elem.tail:
            elem.tail = elem.tail.strip()
        for child in elem:
            strip_whitespace_nodes(child)

    strip_whitespace_nodes(root)

    rough = ET.tostring(root, encoding="utf-8")
    reparsed = xml.dom.minidom.parseString(rough)
    pretty = reparsed.toprettyxml(indent="  ")

    tracers_xml_out.parent.mkdir(parents=True, exist_ok=True)
    with open(tracers_xml_out, "w", encoding="utf-8") as f:
        f.write(pretty)
