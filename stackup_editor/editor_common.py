########################################################################
#
# Copyright 2025-2026 Volker Muehlhaus and IHP PDK Authors
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

"""
editor_common.py

The stackup cross-section preview widget (VectorWidget) and its default color/label
helpers, extracted from setupEM's setup_common.py for use by this repo's own,
independent Stackup Editor. setup_common.py is a large module with a
gds2palace/gdspy/scipy/requests-heavy dependency footprint that this editor doesn't
need - this file carries over only the 4 symbols stackupEditor.py actually uses:
VectorWidget, epsilon_to_color, default_stackup_dielectric_label,
default_stackup_metal_label.

VectorWidget draws directly from resolved dielectric_layer/metal_layer/
stackup_material objects (attributes like .zmin/.zmax/.thickness/.name/.material/
.eps/.sigma/.Rs/.type/.is_via/.is_dielectric/.is_sheet/.get_planar_metals_inside()) -
these come from gds2openEMS's own util_stackup_reader.py, whose object model is
identical to gds2palace's in this respect (both readers are, in fact, byte-for-byte
identical copies as of this port - see openems_ihp_sg13g2's util_stackup_reader.py),
so no attribute-name adaptation was needed here.

The one real dependency change: the original used scipy.interpolate.interp1d(kind=
"linear", fill_value="extrapolate") for the z-to-screen-y mapping. gds2openEMS does
not depend on scipy, and pulling it in just for this one bounded linear interpolation
is unwarranted, so _linear_interp_extrapolate() below reimplements exactly that
behavior (piecewise-linear between the two nearest stored points, linear extrapolation
using the nearest two points beyond either end) without adding a new dependency.
"""

import numpy as np
from PySide6.QtWidgets import QWidget
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtCore import Qt, QRect


# ---------- STACKUP PREVIEW COLOR/LABEL DEFAULTS (permittivity-based) ------------------

def epsilon_to_color(erel, transparency):
    # Compute raw float components
    red   = 250 - 30 * (erel - 1)
    green = 255 - 20 * (erel - 1) + (20 / erel) + 10 * erel
    blue  = 100 + 15 * erel + (250 / erel)

    # Extra adjustment
    if 3.8 < erel < 4.5:
        red   += 50 * (erel - 3.8)
        green -= 100 * (erel - 3.8)

    # Clamp to range 0–255
    red   = min(max(red,   0), 255)
    green = min(max(green, 0), 255)
    blue  = min(max(blue,  0), 255)

    # Convert to integer RGB
    r = int(round(red))
    g = int(round(green))
    b = int(round(blue))

    return QColor(r, g, b, transparency)


def default_stackup_dielectric_label(dielectric, material):
    material_string = f'εr={material.eps:.1f}'
    if material.sigma > 1e-3:
        material_string = material_string + f' σ={material.sigma:.1f}'
    material_string = material_string + f'\n{dielectric.thickness:.2f}µm'
    return material_string


def default_stackup_metal_label(metal, material, is_sheet):
    if is_sheet:
        return f'Rs={material.Rs*1e3:.1f}mΩ'
    else:
        if (material.sigma > 0) and (metal.thickness > 0):
            Rs = 1 / (material.sigma*metal.thickness*1e-6)
            if Rs < 1:
                return f'Rs={Rs*1e3:.1f} mΩ'
            else:
                return f'Rs={Rs:.2f} Ω'
        else:
            return '? ' + material.type + ' ?'


# ---------- linear interpolation/extrapolation (scipy-free replacement for
#             interp1d(kind="linear", fill_value="extrapolate")) ------------------

class _LinearInterpExtrapolate:
    """Callable piecewise-linear interpolator over sorted (x, y) points, with linear
       extrapolation beyond either end - a drop-in replacement for
       scipy.interpolate.interp1d(x, y, kind="linear", fill_value="extrapolate") for the
       one bounded use VectorWidget needs (mapping a stackup z position to a screen y
       coordinate), so this module has no scipy dependency.
    Args:
        x_sorted (numpy.ndarray): x values, strictly sorted ascending, len >= 2
        y_sorted (numpy.ndarray): corresponding y values, same length as x_sorted
    """

    def __init__(self, x_sorted, y_sorted):
        self.x_sorted = x_sorted
        self.y_sorted = y_sorted

    def __call__(self, x_query):
        x_sorted = self.x_sorted
        y_sorted = self.y_sorted
        if x_query <= x_sorted[0]:
            x0, x1, y0, y1 = x_sorted[0], x_sorted[1], y_sorted[0], y_sorted[1]
        elif x_query >= x_sorted[-1]:
            x0, x1, y0, y1 = x_sorted[-2], x_sorted[-1], y_sorted[-2], y_sorted[-1]
        else:
            idx = int(np.searchsorted(x_sorted, x_query))
            x0, x1, y0, y1 = x_sorted[idx - 1], x_sorted[idx], y_sorted[idx - 1], y_sorted[idx]
        if x1 == x0:
            return y0
        t = (x_query - x0) / (x1 - x0)
        return y0 + t * (y1 - y0)


# ---------- POP UP WINDOW TO SHOW STACKUP ------------------

class VectorWidget(QWidget):
    """This widget actually draws the stackup preview.

    The color/label logic for dielectrics and metals is injected as callables
    instead of being hardcoded here, so a host application can customize it
    (e.g. a permittivity-based preview vs. a thermal-conductivity-based one):

        dielectric_color_fn(material) -> QColor
        dielectric_label_fn(dielectric, material) -> str
        metal_label_fn(metal, material, is_sheet) -> str
        via_label_suffix_fn(metal, material) -> str
    """

    def __init__(self, materials_list, dielectrics_list, metals_list,
                 dielectric_color_fn, dielectric_label_fn,
                 metal_label_fn, via_label_suffix_fn):
        super().__init__()
        self.materials_list = materials_list
        self.dielectrics_list = dielectrics_list
        self.metals_list = metals_list
        self.dielectric_color_fn = dielectric_color_fn
        self.dielectric_label_fn = dielectric_label_fn
        self.metal_label_fn = metal_label_fn
        self.via_label_suffix_fn = via_label_suffix_fn

    def paintEvent(self, event):

        # utility: flip y to have y=0 at bottom
        def flipy(y):
            return self.height() - y

        # utility to draw text with alignment on right side
        def drawText_right(x, y, w, h, text):
            rect = QRect(x, y - h, w, h)
            painter.drawText(rect, Qt.AlignVCenter | Qt.AlignRight, text)

        def drawText_left(x, y, w, h, text):
            rect = QRect(x, y - h, w, h)
            painter.drawText(rect, Qt.AlignVCenter | Qt.AlignLeft, text)

        painter = QPainter(self)
        painter.fillRect(self.rect(), Qt.white)
        painter.setRenderHint(QPainter.Antialiasing)

        xmin = int(self.width() * 0.02)
        xmax = int(self.width() * 0.98)

        ymin = int(self.height() * 0.025)
        ymax = int(self.height() * 0.975)

        penBlack = QPen(Qt.black, 1)
        penGray = QPen(QColor(134, 132, 130))
        penDarkGray = QPen(QColor(53, 50, 47))

        # get total dielectric parts, where each metal in a dielectric adds one part
        dielectric_shapes = []
        total_parts = 0
        # sorted by resolved zmin, not just reversed file/array order: a Reference-based
        # dielectric's actual position comes from resolving its Reference by name (see
        # dielectric_layers_list.resolve_references()), entirely independent of where it
        # sits in the file - so reordering it there (e.g. Move Up/Down in the Dielectric
        # Stack tab) must not change where it's drawn here, even though it does change
        # self.dielectrics_list.dielectrics' own array order
        dielectrics_bottom_up = sorted(self.dielectrics_list.dielectrics, key=lambda d: d.zmin)
        for dielectric in dielectrics_bottom_up:  # bottom up
            painter.setPen(penBlack)

            metals_inside = dielectric.get_planar_metals_inside()
            # get number of unique zmin values in that list
            zmin_list = []
            for metal in metals_inside:
                if not metal.zmin in zmin_list:
                    zmin_list.append(metal.zmin)
            metals_count = len(zmin_list)

            # first metal not aligned with dielectric?
            if len(metals_inside) > 0:
                if metals_inside[0].zmin > dielectric.zmin:
                    metals_count = metals_count + 0.5

            parts = max(1, metals_count)
            dielectric_shape = {}
            dielectric_shape['name'] = dielectric.name
            dielectric_shape['dielectric'] = dielectric
            dielectric_shape['numparts'] = parts

            materialname = dielectric.material
            material = self.materials_list.get_by_name(materialname)
            # dielectric color/label are app-specific (permittivity vs. thermal conductivity)
            dielectric_shape['color'] = self.dielectric_color_fn(material)
            dielectric_shape['material'] = material

            total_parts = total_parts + parts
            dielectric_shapes.append(dielectric_shape)

        # calculate height of one dielectric shape
        total_parts = max(total_parts, 1)
        part_height = int((ymax - ymin) / (total_parts))

        y = ymin
        w = xmax - ymin

        # we need to store data for original z position and the displayed y position
        stored_z = np.array([0])
        stored_y = np.array([ymin])

        for dielectric_shape in dielectric_shapes:
            h = part_height * dielectric_shape['numparts']
            dielectric = dielectric_shape['dielectric']
            color = dielectric_shape['color']
            material = dielectric_shape['material']

            material_string = self.dielectric_label_fn(dielectric, material)

            painter.setPen(penBlack)
            painter.setBrush(color)
            painter.drawRect(xmin, flipy(y), w, -h)
            drawText_left(xmin + 5, flipy(y), w, h, dielectric.name)
            drawText_right(xmin, flipy(y), w - 5, h, material_string)

            if not dielectric.zmax in stored_z:
                stored_z = np.append(stored_z, dielectric.zmax)
                stored_y = np.append(stored_y, y + h)

            # get metals inside this dielectric
            metals_inside = dielectric.get_planar_metals_inside()
            # height for one dielectric segment including one metal is part_height
            if len(metals_inside) > 0:

                # there could be multiple metals starting at the same zmin

                # draw planar metals, one after another
                ymetal = y
                for n, metal in enumerate(metals_inside):

                    painter.setPen(penBlack)

                    # check if metal is aligned with dielectric zmin
                    elevation = metal.zmin - dielectric.zmin
                    if n == 0 and (abs(elevation) > 0.001):
                        # draw some vertical offset, not aligned with dielectric
                        ymetal = ymetal + part_height * 0.5  # slight offset

                    # check if next metal is at same zmin
                    next_at_same_zmin = False
                    previous_at_same_zmin = False
                    xmetal = xmin + 120
                    wmetal = w - 200

                    if n < len(metals_inside) - 1:
                        next_metal = metals_inside[n + 1]
                        if abs(next_metal.zmin - metal.zmin) < 0.001:
                            next_at_same_zmin = True
                            xmetal = xmin + 120
                            wmetal = int(w / 2) - 100
                    else:
                        next_metal = None

                    # for the "distance to metal above" label below: several metals
                    # can share this zmin (e.g. sheet resistors drawn side by side),
                    # so skip past all of them to the first one that's actually at a
                    # different (higher) zmin - next_metal above is only the very next
                    # list entry, which for a same-zmin sibling would wrongly give 0
                    next_metal_above = None
                    for candidate in metals_inside[n + 1:]:
                        if abs(candidate.zmin - metal.zmin) >= 0.001:
                            next_metal_above = candidate
                            break

                    if n > 0:
                        previous_metal = metals_inside[n - 1]
                        if abs(previous_metal.zmin - metal.zmin) < 0.001:
                            xmetal = xmin + int(w / 2) + 20
                            wmetal = int(w / 2) - 100
                            previous_at_same_zmin = True

                    material = self.materials_list.get_by_name(metal.material)
                    if material is not None:
                        if metal.is_sheet:
                            # sheet metal that is simulated with zero extrusion
                            height = 3
                            label_string = self.metal_label_fn(metal, material, True)
                        else:
                            # regular extruded metal
                            height = part_height / 2
                            label_string = self.metal_label_fn(metal, material, False)

                        # the box for this metal
                        if material.type.upper() == "CONDUCTOR":
                            painter.setBrush(QColor(230, 230, 230, 90))
                            painter.drawRect(xmetal, flipy(ymetal), wmetal, -int(height))
                        else:
                            painter.setBrush(QColor(230, 130, 130, 90))
                            painter.drawRect(xmetal, flipy(ymetal), wmetal, -int(height))
                    else:
                        # material assignment is invalid
                        height = part_height / 2
                        painter.setBrush(QColor(255, 0, 0, 80))
                        painter.drawRect(xmetal, flipy(ymetal), wmetal, -int(height))
                        label_string = 'INVALID MATERIAL REFERENCE: ' + metal.material

                    painter.setPen(penBlack)
                    drawText_left(xmetal + 10, flipy(ymetal), wmetal, part_height / 2, f"{metal.name} ({metal.layernum})")
                    painter.setPen(penGray)
                    drawText_right(xmetal, flipy(ymetal), wmetal - 10, part_height / 2, label_string)
                    # store the drawing position, because vias will refer to that
                    if not metal.zmin in stored_z:
                        stored_z = np.append(stored_z, metal.zmin)
                        stored_y = np.append(stored_y, ymetal)
                    if not metal.zmax in stored_z:
                        stored_z = np.append(stored_z, metal.zmax)
                        stored_y = np.append(stored_y, ymetal + height)

                    painter.setPen(penGray)
                    painter.drawLine(xmetal - 60, flipy(ymetal), xmetal - 10, flipy(ymetal))
                    # draw line at top side of metal
                    if not metal.is_sheet:
                        painter.drawLine(xmetal - 60, flipy(ymetal + height), xmetal - 10, flipy(ymetal + height))
                        heightstring = f'{metal.thickness:.3f}µm'
                        painter.setPen(penDarkGray)
                        drawText_left(xmetal - 60, flipy(ymetal), 50, height, heightstring)

                    if not previous_at_same_zmin:
                        # draw height to metal above
                        if next_metal_above is not None:
                            dz = abs(next_metal_above.zmin - metal.zmax)
                            heightstring = f'{dz:.3f}µm'
                            painter.setPen(penGray)
                            # sheet metals draw at height=3px, too short to fit this
                            # label without vertical clipping - give the text its own
                            # minimum box height, independent of the drawn box height
                            text_height = max(height, 14)
                            drawText_left(xmetal - 60, flipy(ymetal + height), 50, text_height, heightstring)

                    if n == len(metals_inside) - 1:
                        # last metal (top metal)
                        # place text for distance to dielectric boundary

                        painter.setPen(penBlack)
                        # a metal is registered "inside" a dielectric by its zmin alone
                        # (see util_stackup_reader.register_metals_inside()) - its zmax
                        # can legitimately extend past that dielectric's own zmax into
                        # the one(s) above (e.g. a thick metal sitting in a very thin
                        # dielectric slab), which would otherwise show as a negative,
                        # confusingly-worded "distance to the boundary above". Floor at
                        # 0 - the metal is still drawn at its correct position/height,
                        # this only affects this one label.
                        dz = max(0.0, dielectric.zmax - metal.zmax)
                        if dz > 10:
                            heightstring = f'{dz:.1f}µm'
                        else:
                            heightstring = f'{dz:.3f}µm'
                        painter.setPen(penGray)
                        painter.drawText(xmetal - 60, flipy(ymetal + height + 5), heightstring)

                    if n == 0 and elevation > 0.001:
                        # metal not aligned with bottom of dielectric, add a label for offset value
                        heightstring = f'{elevation:.3f}µm'
                        painter.setPen(penGray)
                        painter.drawText(xmetal - 60, flipy(ymetal - 10), heightstring)

                    if not next_at_same_zmin:
                        # increase screen y for next metal
                        ymetal = ymetal + part_height

            y = y + h

        # sort stored positions
        if len(stored_z) > 2:
            idx = np.argsort(stored_z)
            y_sorted = stored_y[idx]
            z_sorted = stored_z[idx]
            # linear, not cubic: the z->y mapping is a layout position (screen height
            # per dielectric is set by how many metals are stacked inside it, not by
            # its physical thickness), so slope can change drastically between
            # consecutive stored points - e.g. a thick, metal-free substrate maps to
            # almost no screen height while a thin, via-packed dielectric maps to a
            # lot. A cubic spline through data like that readily overshoots (Runge's
            # phenomenon), and with unbounded extrapolation that overshoot is
            # unbounded - enough to overflow the int coordinates drawRect() needs
            # below. Linear interpolation/extrapolation is bounded by construction.
            z_to_y = _LinearInterpExtrapolate(z_sorted, y_sorted)

            # next we draw the vias, based on the screen position of metals that we have stored
            # via position alternates between 3 positions along x axis
            pos = 1
            w = (xmax - xmin) / 10

            painter.setBrush(QColor(136, 192, 200, 80))
            for metal in self.metals_list.metals:
                if metal.is_via or metal.is_dielectric:

                    material = self.materials_list.get_by_name(metal.material)
                    label_suffix = self.via_label_suffix_fn(metal, material)

                    y1 = z_to_y(metal.zmin)
                    y2 = z_to_y(metal.zmax)
                    h = abs(y2 - y1)

                    if pos == 1:
                        xvia = (xmax + xmin) / 2 - 4 * w / 2
                        pos = 2
                    elif pos == 2:
                        xvia = (xmax + xmin) / 2 - w / 2
                        pos = 3
                    else:
                        xvia = (xmax + xmin) / 2 + w
                        pos = 1

                    painter.setPen(penBlack)
                    painter.drawRect(xvia, flipy(y1), w, -h)
                    painter.drawText(xvia + 5, flipy(y1 + 5), f"{metal.name} ({metal.layernum})" + label_suffix)

        painter.end()
