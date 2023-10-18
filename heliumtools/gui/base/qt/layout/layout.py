#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-08-01 15:59:24 srlab"

#  file       layout.py
#  copyright  (c) Sebastian Blatt 2017, 2018

from heliumtools.gui.base.qt import qt

DEFAULT_SPACING = 5
DEFAULT_MARGINS = 0


def apply_layout_spacings(layout, spacing, margins):
    layout.setSpacing(spacing)
    if isinstance(margins, int):
        layout.setContentsMargins(margins, margins, margins, margins)
    else:
        layout.setContentsMargins(margins[0], margins[1], margins[2], margins[3])
    return layout


class GridLayout(qt.QtWidgets.QGridLayout):
    def __init__(self, spacing=DEFAULT_SPACING, margins=DEFAULT_MARGINS):
        super().__init__()
        apply_layout_spacings(self, spacing, margins)


class HBoxLayout(qt.QtWidgets.QHBoxLayout):
    def __init__(self, spacing=DEFAULT_SPACING, margins=DEFAULT_MARGINS):
        super().__init__()
        apply_layout_spacings(self, spacing, margins)


class VBoxLayout(qt.QtWidgets.QVBoxLayout):
    def __init__(self, spacing=DEFAULT_SPACING, margins=DEFAULT_MARGINS):
        super().__init__()
        apply_layout_spacings(self, spacing, margins)


class FormLayout(qt.QtWidgets.QFormLayout):
    def __init__(self, spacing=DEFAULT_SPACING, margins=DEFAULT_MARGINS):
        super().__init__()
        apply_layout_spacings(self, spacing, margins)


def apply_fixed_horizontal_size_policy(widget):
    sp = qt.QtWidgets.QSizePolicy(
        qt.QtWidgets.QSizePolicy.Fixed, qt.QtWidgets.QSizePolicy.Preferred
    )
    widget.setSizePolicy(sp)


def apply_fixed_vertical_size_policy(widget):
    sp = qt.QtWidgets.QSizePolicy(
        qt.QtWidgets.QSizePolicy.Preferred, qt.QtWidgets.QSizePolicy.Fixed
    )
    widget.setSizePolicy(sp)


def apply_stretch_vertical_size_policy(widget, stretch):
    sp = qt.QtWidgets.QSizePolicy(
        qt.QtWidgets.QSizePolicy.Preferred, qt.QtWidgets.QSizePolicy.Preferred
    )
    sp.setVerticalStretch(stretch)
    widget.setSizePolicy(sp)


# layout.py ends here
