#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-03-01 13:42:09 sb"

#  file       about.py
#  copyright  (c) Sebastian Blatt 2017, 2018

import os
import sys

from heliumtools.gui.base.qt import qt
from heliumtools.gui.base.qt.icons import icons
from heliumtools.gui.base.qt.layout import layout


class AboutDialog(qt.QtWidgets.QDialog):
    def __init__(self, parent, model):
        super().__init__(parent)

        self._model = model
        self._name = model.get_config("app", "name")
        self._version = model.get_config("app", "version")
        self._description = model.get_config("app", "description")
        self._authors = model.get_config("app", "authors")
        self._emails = model.get_config("app", "emails")

        self._initialize_ui()

    def _create_label(
        self, text, font_type="Sans", font_size=None, bold=False, italic=False
    ):
        label = qt.QtWidgets.QLabel(text, parent=self)

        font = qt.QtGui.QFont(font_type)
        font.setBold(bold)
        font.setItalic(italic)
        if font_size is not None:
            font.setPointSize(font_size)
        label.setFont(font)
        return label

    def _initialize_ui(self):
        self.setWindowTitle("About {}".format(self._name))
        self.setMinimumSize(480, 0)

        h_grid = layout.HBoxLayout(20, 20)
        v_grid = layout.VBoxLayout(10, (0, 20, 0, 0))

        icon_widget = icons.PixmapLabel(self, "qcontrol3_logo", "128.png")
        h_grid.addWidget(icon_widget)

        v_grid.addWidget(self._create_label(self._name, font_size=24, bold=True))

        v_grid.addWidget(self._create_label(self._version))
        v_grid.addWidget(self._create_label(self._description))
        for author, email in zip(self._authors, self._emails):
            text = "{} <{}>".format(author, email)
            v_grid.addWidget(self._create_label(text))

        v_grid.addStretch(1)
        h_grid.addLayout(v_grid)

        self.setLayout(h_grid)


# about.py ends here
