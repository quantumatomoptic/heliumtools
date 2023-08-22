#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-03-01 16:41:57 srlab"

#  file       page.py
#  author     Sebastian Blatt
#  created    2018-03-01 16:15:34
#
#  Copyright (C) 2011 -- 2016 Christoph Gohle, Christian Gross
#  Copyright (C) 2016 -- 2018 Christoph Gohle, Christian Gross, Sebastian Blatt
#
#  This file is part of qcontrol -- The MPQ-developed, python-based
#  experiment control program.
#
#  License:
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Commentary:
#

from heliumtools.gui.base.qt import qt
from heliumtools.gui.base.qt.layout import layout
from heliumtools.gui.base.qt.widget.logging import logtable

from . import userwidget


class MainWidgetPage(qt.QtWidgets.QWidget):
    aboutToCloseEvent = qt.QtCore.pyqtSignal()
    NAME = "Main"

    def __init__(self, parent, model):
        super().__init__(parent)
        self._model = model
        self._initialize_ui()
        parent.aboutToCloseEvent.connect(self._on_about_to_close)

    @property
    def name(self):
        return self.NAME

    def _on_about_to_close(self):
        self._save_geometry()
        self.aboutToCloseEvent.emit()

    def _initialize_ui(self):
        pass

    def _save_geometry(self):
        pass


class VerticalSplitPage(MainWidgetPage):
    NAME = "VSplit"
    CONFIG_SPLITTER_SIZES = ("gui", "vsplit_page_splitter_sizes")

    def __init__(self, parent, model):
        self._splitter = None
        super().__init__(parent, model)

    def _initialize_ui(self):
        top = self._create_top_widget()
        bottom = self._create_bottom_widget()

        splitter = qt.QtWidgets.QSplitter(qt.QtCore.Qt.Vertical)
        splitter.addWidget(top)
        splitter.addWidget(bottom)
        splitter.setSizes(self._model.get_config(*self.CONFIG_SPLITTER_SIZES))

        self._splitter = splitter
        self._top = top
        self._bottom = bottom

        grid = layout.VBoxLayout()
        grid.addWidget(splitter)
        self.setLayout(grid)

    def _create_top_widget(self):
        return qt.QtWidgets.QWidget(self)

    def _create_bottom_widget(self):
        return qt.QtWidgets.QWidget(self)

    def _save_geometry(self):
        self._model.set_config(
            self.CONFIG_SPLITTER_SIZES[0],
            self.CONFIG_SPLITTER_SIZES[1],
            self._splitter.sizes(),
        )


class HorizontalSplitPage(MainWidgetPage):
    NAME = "HSplit"
    CONFIG_SPLITTER_SIZES = ("gui", "hsplit_page_splitter_sizes")

    def __init__(self, parent, model):
        self._splitter = None
        super().__init__(parent, model)

    def _initialize_ui(self):
        left = self._create_left_widget()
        right = self._create_right_widget()

        splitter = qt.QtWidgets.QSplitter()
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes(self._model.get_config(*self.CONFIG_SPLITTER_SIZES))

        self._splitter = splitter
        self._left = left
        self._right = right

        grid = layout.HBoxLayout()
        grid.addWidget(splitter)
        self.setLayout(grid)

    def _create_left_widget(self):
        return qt.QtWidgets.QWidget(self)

    def _create_right_widget(self):
        return qt.QtWidgets.QWidget(self)

    def _save_geometry(self):
        self._model.set_config(
            self.CONFIG_SPLITTER_SIZES[0],
            self.CONFIG_SPLITTER_SIZES[1],
            self._splitter.sizes(),
        )


class MainPage(VerticalSplitPage):
    NAME = "&Main"
    CONFIG_SPLITTER_SIZES = ("gui", "main_page_splitter_sizes")

    def __init__(self, parent, model):
        super().__init__(parent, model)

    def _create_top_widget(self):
        return userwidget.UserWidget(self, self._model)

    def _create_bottom_widget(self):
        storage = self._model._log_storage
        formatter = logtable.LogFormatterProcessed()
        logfilter = logtable.LogFilter()
        table_model = logtable.LogTableModel(
            storage, formatter, logfilter, storage.max_records()
        )
        w = logtable.LogTable(self, table_model)
        return w


class DebugPage(MainWidgetPage):
    NAME = "&Debug console"

    def _initialize_ui(self):
        storage = self._model._log_storage
        formatter = logtable.LogFormatterProcessed()
        logfilter = logtable.LogFilter()
        table_model = logtable.LogTableModel(
            storage, formatter, logfilter, storage.max_records()
        )
        w = logtable.LogTable(self, table_model)
        w.doubleClicked.connect(self._on_click)
        grid = layout.VBoxLayout()
        grid.addWidget(w)
        self.setLayout(grid)

    def _on_click(self, item):
        message = self._model._log_storage.get(item.row())["msg"]
        qt.QtWidgets.QMessageBox.about(self, "Log message", message)


class PreferencesPage(MainWidgetPage):
    NAME = "&Preferences"

    def __init__(self, parent, model):
        super().__init__(parent, model)

    def _initialize_ui(self):
        pass

    def _save_geometry(self):
        pass


# page.py ends here
