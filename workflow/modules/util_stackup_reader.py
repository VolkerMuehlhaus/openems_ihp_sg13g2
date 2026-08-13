#
# Copyright 2025 Volker Muehlhaus and IHP PDK Authors
#
# Licensed under the GNU General Public License, Version 3.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    https://www.gnu.org/licenses/gpl-3.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
########################################################################


# Read XML file with SG13G2 stackup

# File history: 
# Initial version 20 Nov 2024  Volker Muehlhaus 
# Added support for sheet resistance 07 Oct 2025 Volker Muehlhaus 
# Added docstrings 
# 20 Nov 2025: added functionality to get relative positions between metals
# 10 Aug 2026: added derived layer option
# 10 Aug 2026: added Oversize option for derived layers (grow/shrink outline)

__version__ = "1.3.0"

import os
import math
import xml.etree.ElementTree 


def safe_get (data, key, default):
  val = data.get(key)
  if val is not None:
    return val
  else:  
    return default

def _make_comment_preserving_parser():
  """XML parser that keeps <!-- comments --> as Comment nodes in the tree, instead of
     silently dropping them (the default xml.etree behavior). Needed so that editors
     which load, edit, and save a stackup file back to disk don't lose comments that
     live inside sections they don't touch (e.g. DerivedLayers).
  """
  target = xml.etree.ElementTree.TreeBuilder(insert_comments=True)
  return xml.etree.ElementTree.XMLParser(target=target)


# -------------------- material types ---------------------------

class stackup_material:
  """
    stackup material object, can be dielectric or metal with conductivity or sheet with Ohm/square
  """
    
  def __init__ (self, data):
    """create stackup material object from XML data line

    Args:
        data (string): line from XML data, required parameters are "Name" and "Type" strings. Optional: "Permittivity","DielectricLossTangent","Conductivity","Rs","Color"
    """


    self.name = data.get("Name")
    self.type = data.get("Type").upper()
    
    self.eps   = float(safe_get(data, "Permittivity", 1))
    self.tand  = float(safe_get(data, "DielectricLossTangent", 0))
    self.sigma = float(safe_get(data, "Conductivity", 0))
    self.Rs    = float(safe_get(data, "Rs", 0))
    self.density = float(safe_get(data, "Density", 1))
    self.color = data.get("Color")  # no default here, will be handled later 

    self.thermalcond = float(safe_get(data, "ThermalConductivity", 0))
    self.thermaltablename = safe_get(data, "ThermalConductivityTable", "")
    self.thermaltable = None


  def __str__ (self):
    """String representation of stackup_material, useful for debugging

    Returns:
        string: String representation of stackup_material
    """
    # string representation 
    mystr = '      Material Name=' + self.name + ' Type=' + self.type +' Permittivity=' + str(self.eps) + ' DielectricLossTangent=' +  str(self.tand) + ' Conductivity=' +  str(self.sigma) + ' ThermalConductivity=' +  str(self.thermalcond)  + ' Color = ' + self.color
    return mystr



class stackup_materials_list:
  """
    structure with list of stackup material objects (.materials) and value of maximum permittivy (.eps_max)
  """

  def __init__ (self):
    """Create empty structure. 
    """
    self.materials = []      # list with material objects
    self.eps_max   = 0
    
  def append (self, material):
    """Append one material
    Args:
        material (stackup_material): material to add
    """

    # append material
    self.materials.append (material)
    # set maximum permittivity in model
    self.eps_max = max(self.eps_max, material.eps)
  

  def get_by_name (self, materialname):
    """find material object from materialname
    Args:
        materialname (string): Name as specified in XML data line
    Returns:
        stackup_material: the material with that name
    """
  
    # find material object from materialname
    found = None
    for material in self.materials:
      if material.name.upper() == materialname.upper():
        found = material
    return found    


# -------------------- dielectrics ---------------------------

class dielectric_layer:
  """
    dielectric layer object. Holds information on stackup layers that are always there, without drawing them explicitely in GDSII
  """
    
  def __init__ (self, data):
    """create stackup layer object (usually dielectric or semiconductor) from XML data line

    Args:
        data (string): line from XML data, required parameters: "Name","Material","Thickness", optional parameter "Boundary" for bounding layer number
    """
    self.name = data.get("Name")
    self.material = data.get("Material")

    # dielectrics can be specified by thickness when stacked on top of each other, or by absolute zmin/zmax otherwise
    self.zmin = safe_get(data, "Zmin", None)
    self.zmax = safe_get(data, "Zmax", None)
    if not (self.zmin is None or self.zmax is None):
      # we have a valid position, use that instead of stacking everything one after another
      self.zmin = float(self.zmin)
      self.zmax = float(self.zmax)
      self.thickness = self.zmax - self.zmin
      self.absolute_zpos = True
    else:
      # No absolute zmin and zmax, position results from stacking dielectric by order in file, using their thickness values
      # z Position will be set later, by stacking dielectrics on top of each other
      self.zmin = None
      self.zmax = None
      self.thickness  = float(data.get("Thickness"))
      self.absolute_zpos = False

    self.is_top = False
    self.is_bottom = False
    self.gdsboundary = data.get("Boundary")  # optional entry in stackup file

    self.metals_inside = [] # metals that are located inside this dielectric, set by function 

  def get_planar_metals_inside (self):
    """evaluates metals_inside list and returns only items that are conductor or sheet (no via, not dielectric via)
    Returns:
        list of metal_layer: metals that are conductor or sheet
    """
    planar_metals = []
    for metal in self.metals_inside:
      if metal.is_metal or metal.is_sheet:
        planar_metals.append(metal)
    return planar_metals


  def __str__ (self):
    """String representation of dielectric_layer, useful for debugging
    Returns:
        string: String representation of stackup_material
    """
    enclosed_metal_names = []
    for metal in self.metals_inside:
      enclosed_metal_names.append(metal.name)

    mystr = '      Dielectric Name=' + self.name + ' Material=' + self.material +' Thickness=' \
            + str(self.thickness) + ' Zmin=' +  str(self.zmin) + ' Zmax=' +  str(self.zmax) \
            + ' Metals inside: ' + str(enclosed_metal_names)
            
    return mystr



class dielectric_layers_list:
  """
    list that holds all dielectric layer objects
  """

  def __init__ (self):
    """Initialize empty list
    """
    self.dielectrics = []      # list with dielectric objects
    
  def append (self, dielectric, materials_list ):
    """Append one dielectric to the list

    Args:
        dielectric (dielectric_layer): the dielectric that is appended
        materials_list (_type_): not used
    """

    self.dielectrics.append (dielectric)


  def calculate_zpositions (self):
    """dielectrics in XML are in reverse order, so we need to build position upside down
    """

    z = 0
    for dielectric in reversed(self.dielectrics):
      if (dielectric.zmin is None) and (dielectric.zmax is None):
        # only handle layers that don't have zmin, zmax specified in the XML file
        t = float(dielectric.thickness)
        dielectric.zmin = z
        dielectric.zmax = z + t
        z = dielectric.zmax


  def get_by_name (self, name_to_find):  
    """find material object from materialname
    Args:
        name_to_find (string): name of material to find
    Returns:
        dielectric_layer: dielectric with that name, otherwise None
    """

    found = None
    for dielectric in self.dielectrics:
      if dielectric.name.upper() ==  name_to_find.upper():
        found = dielectric
    return found    


  def get_boundary_layers (self):
    """For substrates where Boundary is specified in dielectric layers, return a list of those layers. This is required for the next step, GDSII reader, which needs to know the layers to read. 
    Returns:
        list of int: list of layer numbers specified as boundary
    """
    
    boundary_layer_list = []
    for dielectric in self.dielectrics:
      if dielectric.gdsboundary is not None:
        value = int(dielectric.gdsboundary)
        if value not in boundary_layer_list:
          boundary_layer_list.append(value) 
    return boundary_layer_list
  

  def register_metals_inside (self, metals_list):
    """iterates over dielectrics and metals, sets metals_inside property for each dielectric with list of metals within that z range
    Args:
        metals_list (metal_layers_list): metals read from stackup
    """
    for dielectric in self.dielectrics:
      enclosed = []
      for metal in metals_list.metals:
        # check if metal is enclosed in dielectric, excluding zmax exactly
        if (metal.zmin >= dielectric.zmin) and (metal.zmax < dielectric.zmax):
          enclosed.append(metal)          
        # also include metals that fit exactly the height of the dielectric
        if (metal.zmin == dielectric.zmin) and (metal.zmax == dielectric.zmax):
          enclosed.append(metal)          
      dielectric.metals_inside = enclosed    


# -------------------- conductor layers (metal and via) ---------------------------

class metal_layer:
  """
    drawing layer object ( name metal_layer is misleading, this drawn layer that uses material from the XML materials section)
  """
    
  def __init__ (self, data):
    """create metal layer object (planar metal, via, sheet or dielectric) from XML data line

    Args:
        data (string): line from XML data, required parameters: "Name","Layer","Type","Material","Zmin","Zmax"
       """
    self.name = data.get("Name")
    self.layernum = data.get("Layer")
    self.type = data.get("Type").upper()
    self.material = data.get("Material")
    self.zmin = float(data.get("Zmin"))
    self.zmax = float(data.get("Zmax"))
    
    # force to sheet if zero thickness
    if data.get("Zmin") == data.get("Zmax"):
      self.type = "SHEET"

    if self.type == "SHEET" and not self.zmin==self.zmax:
      print('ERROR: Layer ', self.name, ' is defined as sheet layer, but zmax is different from zmin. This is not valid!')
      exit(1)

    self.thickness = self.zmax-self.zmin
    self.is_via = (self.type=="VIA")
    self.is_metal = (self.type=="CONDUCTOR")
    self.is_dielectric = (self.type=="DIELECTRIC")
    self.is_sheet = (self.type=="SHEET")
    self.is_used = False

    # Metals directly above and below, this is set by metal_layers_list.sort_and_evaluate()
    self.above = []
    self.below = []


  def __str__ (self):
    """String representation of dielectric_layer, useful for debugging
    Returns:
        string: String representation of stackup_material
    """

    # convert list of layers above and below to layer names
    below_names = []
    for layer in self.below:
      below_names.append(layer.name)

    above_names = []
    for layer in self.above:
      above_names.append(layer.name)


    mystr = '      Metal Name=' + self.name + ' Layer=' + self.layernum  + \
            ' Type=' + self.type + ' Material=' + self.material + \
            ' Zmin=' +  str(self.zmin) + ' Zmax=' +  str(self.zmax) + \
            ' below=' + str(below_names) + ' above=' + str(above_names)
    
    return mystr
  



class metal_layers_list:
  """
    list of drawn layer objects (metal, via, dielectric brick)
  """


  def __init__ (self):
    """Initialize emptry list
    """
    self.metals = []      # list with conductor objects
    self.lowest = None    # metal with smallest zmin value
    self.orphan_layers = []  # list with layers that have no direct neighbor above or below
    self.derived_layers = derived_layers_list()  # empty by default, populated by read_substrate() if XML has a DerivedLayers section

  def append (self, metal):
    """Append one metal layer (drawn layer)
    Args:
        metal (metal_layer): metal layer to be added to list
    """
    self.metals.append (metal)

  
  def getbylayernumber (self, number_to_find):
    """Find metal layer by layer number, returns first match
    Args:
        number_to_find (int): layer number to find
    Returns:
        metal_layer: metal layer with that layer number
    """
    
    found = None
    for metal in self.metals:
      if metal.layernum == str(number_to_find):
        found = metal
        break 
    return found  


  def getallbylayernumber (self, number_to_find):
    """returns all metals by layer number as list, finds multiple metals mapped to same number
    Args:
        number_to_find (int): layer number to find
    Returns:
        list: list of metal_layer with that layer number, None if not found
    """
         
    found = []
    for metal in self.metals:
      if metal.layernum == str(number_to_find):
          found.append(metal)
    if found==[]:
      found = None
    return found  


  def getallplanarmetals (self):
    """returns all metals (conductor or sheet) as list, skip vias and dielectric vias
    Returns:
        list: list of metal_layer 
    """
         
    found = []
    for metal in self.metals:
      if metal.is_metal or metal.is_sheet:
          found.append(metal)
    return found  


  def getbylayername (self, name_to_find):
    """Find metal layer by layer number, returns first match
    Args:
        name_to_find (string): layer name to find
    Returns:
        metal_layer: metal layer with that layer name
    """

    found = None
    for metal in self.metals:
      if metal.name == str(name_to_find):
        found = metal
        break 
    return found  


  def getlayernumbers (self):
    """list of all metal and via layer numbers in technology
    Returns:
        list of int: all layer numbers in technology file
    """

    layernumbers = []
    for metal in self.metals:
      layernumbers.append(int(metal.layernum))
    return layernumbers 


  def add_offset (self, offset): 
    """Add offset in z position to all metal layers, used to add stackup height for final z position
    Args:
        offset (float): z offset in project units
    """
    for metal in self.metals:
      metal.zmin = metal.zmin + offset
      metal.zmax = metal.zmax + offset


  def sort_and_evaluate(self):
    """After reading all metals, sort them by position and detect the neighbors above/below
       This is set in each metal as .above and .below list
    """
    # sort the list by zmin of each metal
    self.metals.sort(key=lambda metal: metal.zmin)
    # metal with lowest zmin value
    self.lowest = self.metals[0]
    self.highest = self.metals[-1]

    # delta for comparison, i.e. what is considered equal
    delta = 1e-5

    # Build above/below relationships efficiently
    for i, layer in enumerate(self.metals):
        # Layers above: all layers with zmin >= current zmax
        for other in self.metals[i+1:]:
            if abs(other.zmin - layer.zmax) < delta:
                layer.above.append(other)
        # Layers below: all layers with zmax <= current zmin
        for other in self.metals[:i]:
            if abs(other.zmax - layer.zmin) < delta:
                layer.below.append(other)

    # Identify orphan layers (no above or below)
    self.orphan_layers = [layer for layer in self.metals if not layer.above and not layer.below]


# -------------------- derived layers (boolean operations on other layers) ---------------------------

class derived_layer:
  """
    derived layer object, defines a new layer number that is computed via boolean operation
    on other layers (native GDSII layers or other derived layers), instead of being read
    directly from GDSII. Optionally, the result can be grown or shrunk by a fixed distance
    (Oversize). Operation "SIZE" takes a single operand and applies only the Oversize, with
    no boolean operation.
  """

  VALID_OPERATIONS = ("AND", "OR", "XOR", "NOT", "SIZE")

  def __init__ (self, data):
    """create derived layer object from XML data line

    Args:
        data (xml.etree.ElementTree.Element): "DerivedLayer" XML element, required attributes
          "Name", "Layer", "Operation" ; optional attribute "Oversize" (grow outline by this
          distance in layout units, negative value shrinks); child "Operand" elements with
          "Layer" attribute, in order. For "NOT", first operand minus all following operands.
          Order does not matter for "AND", "OR", "XOR". "SIZE" takes exactly one operand and
          requires a non-zero Oversize; it just resizes that operand onto the new layer number.
    """
    self.name = data.get("Name")
    self.layernum = data.get("Layer")

    operation = data.get("Operation")
    self.operation = operation.upper() if operation is not None else None
    if self.operation not in self.VALID_OPERATIONS:
      print('ERROR: Derived layer ', self.name, ' has invalid Operation "', operation, '". Must be one of ', self.VALID_OPERATIONS)
      exit(1)

    oversize = data.get("Oversize")
    self.oversize = float(oversize) if oversize is not None else 0.0

    self.operands = []
    for operand in data.findall("Operand"):
      self.operands.append(operand.get("Layer"))

    if self.operation == "SIZE":
      if len(self.operands) != 1:
        print('ERROR: Derived layer ', self.name, ' has Operation="SIZE", which needs exactly 1 Operand entry, found ', len(self.operands))
        exit(1)
      if self.oversize == 0:
        print('ERROR: Derived layer ', self.name, ' has Operation="SIZE", which needs a non-zero Oversize value')
        exit(1)
    elif len(self.operands) < 2:
      print('ERROR: Derived layer ', self.name, ' needs at least 2 Operand entries, found ', len(self.operands))
      exit(1)


  def __str__ (self):
    """String representation of derived_layer, useful for debugging
    Returns:
        string: String representation of derived_layer
    """
    mystr = '      DerivedLayer Name=' + self.name + ' Layer=' + self.layernum + \
            ' Operation=' + self.operation + ' Operands=' + str(self.operands) + \
            ' Oversize=' + str(self.oversize)
    return mystr



class derived_layers_list:
  """
    list of derived layer objects. Operands of a derived layer can be native GDSII layers
    or other derived layers, so this class also provides dependency resolution (topological
    sort) to determine a safe processing order.
  """

  def __init__ (self):
    """Initialize empty list
    """
    self.derived_layers = []

  def append (self, derived):
    """Append one derived layer
    Args:
        derived (derived_layer): derived layer to add to list
    """
    self.derived_layers.append (derived)


  def getbylayernumber (self, number_to_find):
    """Find derived layer by layer number
    Args:
        number_to_find (int): layer number to find
    Returns:
        derived_layer: derived layer with that layer number, None if not found
    """
    found = None
    for derived in self.derived_layers:
      if derived.layernum == str(number_to_find):
        found = derived
        break
    return found


  def getlayernumbers (self):
    """list of all derived layer numbers
    Returns:
        list of int: all derived layer numbers
    """
    layernumbers = []
    for derived in self.derived_layers:
      layernumbers.append(int(derived.layernum))
    return layernumbers


  def get_ordered (self):
    """Resolve dependencies between derived layers (a derived layer can use another derived
       layer as operand) and return the list sorted so that a derived layer never appears
       before the derived layers it depends on. Native GDSII layer operands are not
       dependencies here, since they need no prior computation.
    Returns:
        list of derived_layer: topologically sorted derived layers
    Raises:
        SystemExit: if a circular or otherwise unresolvable dependency is detected
    """

    resolved = []
    resolved_layernums = set()
    remaining = list(self.derived_layers)

    while remaining:
      progress = False
      still_remaining = []

      for derived in remaining:
        # dependencies are operands that are themselves derived layers (not native GDSII layers)
        dependencies = [op for op in derived.operands if self.getbylayernumber(op) is not None]
        if all(dep in resolved_layernums for dep in dependencies):
          resolved.append(derived)
          resolved_layernums.add(derived.layernum)
          progress = True
        else:
          still_remaining.append(derived)

      if not progress:
        unresolved_names = [derived.name for derived in still_remaining]
        print('ERROR: Circular or unresolved dependency in DerivedLayers: ', unresolved_names)
        exit(1)

      remaining = still_remaining

    return resolved


# ----------- thermal tables -----------------

class thermal_table:
    def __init__(self, xml_data):
        self.name = xml_data.attrib["Name"]
        self.points = []

        for point in xml_data.iter("Point"):
            T = float(point.attrib["Temperature"])
            k = float(point.attrib["Value"])
            self.points.append((T, k))


class thermal_tables_list(list):
    def get(self, name):
        for table in self:
            if table.name == name:
                return table
        return None


# ----------- parse substrate file, get materials from list created before -----------

DESCRIPTION_COMMENT_PREFIX = "File description:"
# duplicated (not imported) from util_stackup_writer.py deliberately: that module
# is specific to the interactive XML editor, while this reader module is used by
# the whole gds2palace pipeline and is meant to stay independent of it.


def read_file_description (XML_filename):
  """
  Read the free-text file description previously stamped by the XML Stackup Editor
  (see util_stackup_writer.stamp_header_comments), if any. Cheap standalone lookup
  that does not build the full materials/dielectrics/metals object model.
  Args:
      XML_filename (string): filename of XML technology file
  Returns:
      string: the description text, or "" if the file has none / is missing / fails to parse
  """
  if not XML_filename or not os.path.isfile(XML_filename):
    return ""

  try:
    root = xml.etree.ElementTree.parse(XML_filename, parser=_make_comment_preserving_parser()).getroot()
  except xml.etree.ElementTree.ParseError:
    return ""

  for child in root:
    if child.tag is xml.etree.ElementTree.Comment:
      text = (child.text or "").strip()
      if text.startswith(DESCRIPTION_COMMENT_PREFIX):
        return text[len(DESCRIPTION_COMMENT_PREFIX):].strip()
  return ""


def read_substrate (XML_filename):
  """
  Read XML substrate and return materials_list, dielectrics_list, metals_list.
  Derived layer definitions (if any) are attached as metals_list.derived_layers,
  so this return signature stays backward compatible with existing 3-value unpacking.
  Args:
      XML_filename (string): filename of XML technology file
  """

  if os.path.isfile(XML_filename):
    # print('Reading XML stackup  file:', XML_filename)

    # data source is *.subst XML file; comment-preserving parser so that callers who
    # keep the parsed tree around for editing (e.g. a GUI editor) don't lose comments
    substrate_tree = xml.etree.ElementTree.parse(XML_filename, parser=_make_comment_preserving_parser())
    substrate_root = substrate_tree.getroot()

    return parse_substrate(substrate_root)

  else:
    print('XML stackup file not found: ', XML_filename)
    exit(1)


def parse_substrate (substrate_root):
  """
  Build materials_list, dielectrics_list, metals_list from an already-parsed XML
  <Stackup> root element. This is the part of read_substrate() that doesn't touch
  disk, split out so callers that already hold a parsed/edited tree in memory (e.g.
  a GUI editor re-deriving a live preview after an edit) can re-run it without a
  round trip through the filesystem.
  Args:
      substrate_root (xml.etree.ElementTree.Element): root element of a stackup XML tree
  """

  # get materials  from  XML
  materials_list = stackup_materials_list() # initialize empty list
  for data in  substrate_root.iter("Material"):
      materials_list.append (stackup_material(data))

  # get dielectric layers from  XML
  dielectrics_list = dielectric_layers_list() # initialize empty list
  for data in  substrate_root.iter("Dielectric"):
      dielectrics_list.append (dielectric_layer(data), materials_list)

  # get optional thermal tables from XML
  thermal_tables = thermal_tables_list()
  tables = substrate_root.find("Tables")
  if tables is not None:
      for data in tables.findall("Table"):
          thermal_tables.append(thermal_table(data))

  # iterate over materials to check if they refer to a thermal table
  for material in materials_list.materials:
    if material.thermaltablename != "":
      material.thermaltable = thermal_tables.get(material.thermaltablename)
    else:
      material.thermaltable = None

  # calculate z positions in dielectric layers, after reading all of them
  dielectrics_list.calculate_zpositions()

  # mark top and bottom, order from XML is top material first
  if len(dielectrics_list.dielectrics) > 0:
    zmin_overall = math.inf
    zmax_overall = -math.inf
    lowest = None
    highest = None

    for dielectric in dielectrics_list.dielectrics:
      if dielectric.zmin < zmin_overall:
        zmin_overall = dielectric.zmin
        lowest = dielectric
      if dielectric.zmax > zmax_overall:
        zmax_overall = dielectric.zmax
        highest = dielectric

    if highest is not None:
      highest.is_top = True
    if lowest is not None:
      lowest.is_bottom = True

  # get metal layers (metals + vias) from XML
  metals_list = metal_layers_list() # initialize empty list
  for data in  substrate_root.iter("Layer"):
      metals_list.append (metal_layer(data))

  # sort metals by zmin and detect their neighbors above/below
  metals_list.sort_and_evaluate()

  # get derived layers (boolean operations on other layers) from XML, if present
  # attached to metals_list instead of added as a new return value, so that existing
  # code doing "materials_list, dielectrics_list, metals_list = read_substrate(...)" keeps working
  derived_layers_section = substrate_root.find(".//DerivedLayers")
  if derived_layers_section is not None:
    for data in derived_layers_section.findall("DerivedLayer"):
      metals_list.derived_layers.append (derived_layer(data))

  # get substrate offset, required for v2 stackup file version
  offset = 0
  for data in substrate_root.iter("Substrate"):
      assert data!=None
      offset = float(data.get("Offset"))
  if offset > 0:
    metals_list.add_offset(offset)

  # register metals with the enclosing dielectrics
  dielectrics_list.register_metals_inside (metals_list)

  return materials_list, dielectrics_list, metals_list


# =======================================================================================
# Test code when running as standalone script
# =======================================================================================

if __name__ == "__main__":

  XML_filename = "SG13G2_200um.xml"
  materials_list, dielectrics_list, metals_list = read_substrate (XML_filename)
  derived_layers = metals_list.derived_layers

  for material in materials_list.materials:
    print(material)

  for dielectric in dielectrics_list.dielectrics:
    print(dielectric)

  for metal in metals_list.metals:
    print(metal)

  print('__________________________________________')

  # test finding a layer by layer number
  metal = metals_list.getbylayernumber (134)
  print('Layer 134 name => ', metal.name)

  print('Layer 134 thickness => ', metals_list.getbylayernumber (134).thickness)
  print('Test if Layer 134 is a via layer => ', metals_list.getbylayernumber(134).is_via)

  # test finding a layer by name
  metal = metals_list.getbylayername ("TopMetal1")
  print('TopMetal1 layer number => ', metal.layernum)

  # find orphaned layers that have no neighbor above or below
  orphan_names = []
  for layer in metals_list.orphan_layers:
    orphan_names.append(layer.name)
  print('Orphaned layers: ', orphan_names)  

  # get all planar metals in stackup (no via, no dielectric via)
  planar_metal_names = []
  for metal in metals_list.getallplanarmetals():
    planar_metal_names.append(metal.name)
  print('Planar metals: ', planar_metal_names)  

  # get all planar metals inside a dielectric
  DK = dielectrics_list.get_by_name('SiO2')
  if DK is not None:
    names = []
    metals = DK.get_planar_metals_inside()
    for metal in metals:
      names.append(metal.name)
    print('Planar metals inside ', DK.name, ': ', names)

  print('__________________________________________')

  for derived in derived_layers.derived_layers:
    print(derived)

  # test finding a derived layer by layer number
  derived = derived_layers.getbylayernumber (240)
  if derived is not None:
    print('Layer 240 name => ', derived.name)
    print('Layer 240 operation => ', derived.operation)
    print('Layer 240 operands => ', derived.operands)

  # dependency-resolved processing order (derived layers built from other derived layers come last)
  processing_order = [derived.name for derived in derived_layers.get_ordered()]
  print('Derived layer processing order: ', processing_order)
