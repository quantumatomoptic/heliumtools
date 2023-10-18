#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2018-08-27 18:39:26 srlab"

#  file       logtable.py
#  copyright  (c) Sebastian Blatt 2017, 2018


import sys
import copy

from heliumtools.gui.base.qt import qt

from . import signalling_queue


class LogFilter(object):
    def __init__(self):
        pass

    def match(self, record):
        return True


class LogFilterLevelname(object):
    def __init__(self):
        self.levels = {
            "NOTSET": True,
            "DEBUG": True,
            "INFO": True,
            "WARNING": True,
            "WARN": True,
            "ERROR": True,
            "CRITICAL": True,
            "FATAL": True,
        }

    def match(self, record):
        if "levelname" in record and record["levelname"] in self.levels:
            return self.levels[record["levelname"]]
        return False


class LogFormatter(object):
    _record_labels = []
    _column_labels = []

    def __init__(self):
        pass

    def columns(self):
        return len(self._column_labels)

    def record_labels(self):
        return self._record_labels

    def column_labels(self):
        return self._column_labels

    def background_color(self, record, row, column):
        return qt.QtCore.QVariant()

    def foreground_color(self, record, row, column):
        return qt.QtCore.QVariant()


class LogFormatterBare(LogFormatter):
    _record_labels = ["asctime", "bare_text"]
    _column_labels = ["Time", "Message"]

    def __init__(self):
        super().__init__()


class LogFormatterProcessed(LogFormatter):
    _record_labels = [
        "levelname",
        "asctime",
        # 'pathname',
        # 'funcName',
        # 'lineno',
        # 'process',
        # 'threadName',
        "msg",
    ]
    _column_labels = [
        "Level",
        "Time",
        # 'Path',
        # 'Function',
        # 'Line',
        # 'Process',
        # 'Thread',
        "Message",
    ]

    _levelname_colors = {
        # white
        "NOTSET": ((0xff, 0xff, 0xff), (0xee, 0xee, 0xee)),
        "DEBUG": ((0xff, 0xff, 0xff), (0xee, 0xee, 0xee)),
        # green
        "INFO": ((0xa6, 0xd9, 0x6a), (0x66, 0xbd, 0x63)),
        # yellow
        "WARNING": ((0xfd, 0xae, 0x61), (0xfe, 0xe0, 0x8b)),
        "WARN": ((0xfd, 0xae, 0x61), (0xfe, 0xe0, 0x8b)),
        # red
        "ERROR": ((0xf4, 0x6d, 0x43), (0xd7, 0x30, 0x27)),
        "CRITICAL": ((0xf4, 0x6d, 0x43), (0xd7, 0x30, 0x27)),
        "FATAL": ((0xf4, 0x6d, 0x43), (0xd7, 0x30, 0x27)),
    }

    def __init__(self):
        super().__init__()
        self._colors = {}
        for k, (v, w) in self._levelname_colors.items():
            c_even = qt.QtGui.QColor(v[0], v[1], v[2])
            c_odd = qt.QtGui.QColor(w[0], w[1], w[2])
            self._colors[k] = (c_even, c_odd)

    def background_color(self, record, row, column):
        colors = None
        if "levelname" not in record:
            colors = self._colors["NOTSET"]
        else:
            colors = self._colors[record["levelname"]]

        if row % 2 == 0:
            return colors[0]
        return colors[1]

    def foreground_color(self, record, row, column):
        return qt.QtGui.QColor(0x00, 0x00, 0x00)


class LogStorage(object):
    def __init__(self, max_records):
        self._queue = signalling_queue.SignallingQueue(max_records, self._on_update)
        self._viewers = []

    def register_viewer(self, viewer):
        self._viewers.append(viewer)

    def _on_update(self, index_map):
        for viewer in self._viewers:
            viewer.apply_changes(self._queue, index_map)

    def append(self, record):
        self._queue.append(record)

    def extend(self, records):
        self._queue.extend(records)

    def max_records(self):
        return self._queue.max_length()

    def length(self):
        return self._queue.length()

    def apply_filter(self, log_filter):
        rc = []
        for j in range(self._queue.length()):
            record = self._queue.get(j)
            if log_filter.match(record):
                rc.append(j)
        return rc

    def get(self, index):
        return self._queue.get(index)

    def clear(self):
        self._queue.clear()


class LogTableModel(qt.QtCore.QAbstractTableModel):
    def __init__(self, log_storage, log_formatter, log_filter, max_records):
        super().__init__()
        self._log_storage = log_storage
        self._log_formatter = log_formatter
        self._log_filter = log_filter
        self._max_records = max_records
        self._viewed_indices = []

        log_storage.register_viewer(self)

    def set_filter(self, log_filter):
        self.layoutAboutToBeChanged.emit()

        self.beginResetModel()
        self._log_filter = log_filter
        self._viewed_indices = self._log_storage.apply_filter(self._log_filter)
        self.endResetModel()

        index_tl = self.createIndex(0, 0)
        index_br = self.createIndex(self.rows() - 1, self.columns() - 1)
        self.dataChanged.emit(index_tl, index_br)

        self.layoutChanged.emit()

    def rows(self):
        return len(self._viewed_indices)

    def columns(self):
        return self._log_formatter.columns()

    def rowCount(self, index):
        return self.rows()

    def columnCount(self, index):
        return self.columns()

    def data(self, index, role):
        row = index.row()
        col = index.column()

        if role == qt.QtCore.Qt.DisplayRole:
            return self._handle_display(row, col)
        elif role == qt.QtCore.Qt.TextAlignmentRole:
            return self._handle_text_alignment(row, col)
        elif role == qt.QtCore.Qt.ForegroundRole:
            return self._handle_foreground_color(row, col)
        elif role == qt.QtCore.Qt.BackgroundRole:
            return self._handle_background_color(row, col)

        return qt.QtCore.QVariant()

    def _handle_display(self, row, col):
        key = self._log_formatter.record_labels()[col]
        row_index = self._viewed_indices[row]
        record = self._log_storage.get(row_index)
        if key in record:
            return record[key]
        return ""

    def _handle_text_alignment(self, row, col):
        if col == 0:
            return qt.QtCore.Qt.AlignTop | qt.QtCore.Qt.AlignRight
        return qt.QtCore.Qt.AlignTop | qt.QtCore.Qt.AlignLeft

    def _handle_background_color(self, row, column):
        record = self._log_storage.get(self._viewed_indices[row])
        return self._log_formatter.background_color(record, row, column)

    def _handle_foreground_color(self, row, column):
        record = self._log_storage.get(self._viewed_indices[row])
        return self._log_formatter.foreground_color(record, row, column)

    def headerData(self, section, orientation, role):
        if role == qt.QtCore.Qt.DisplayRole:
            if orientation == qt.QtCore.Qt.Horizontal:
                return self._log_formatter.column_labels()[section]
            elif orientation == qt.QtCore.Qt.Vertical:
                pass
        return qt.QtCore.QVariant()

    def _diff_viewed_indices(self, indices_old, indices_new):
        n_old = len(indices_old)
        n_new = len(indices_new)
        if n_old > 0 and n_new > 0 and indices_old[0] != indices_new[0]:
            return 0, (n_new - 1)
        if n_new == 0:
            return 0, 0

        for i, j in zip(indices_old, indices_new):
            if i != j:
                return j, (n_new - 1)

        return 0, (n_new - 1)

    def _index_map_to_row_range(self, index_map):
        min_row = None
        max_row = None
        new_indices = copy.deepcopy(self._viewed_indices)

        for old_index, new_index in index_map:
            if old_index in new_indices:
                idx = new_indices.index(old_index)
                if new_index is None:
                    del new_indices[idx]
                else:
                    new_indices[idx] = new_index
            elif old_index is None:
                record = self._log_storage.get(new_index)
                if self._log_filter.match(record):
                    new_indices.append(new_index)

        # return min_row, max_row
        diff = self._diff_viewed_indices(self._viewed_indices, new_indices)
        self._viewed_indices = new_indices

        return diff

    def apply_changes(self, queue, index_map):
        min_row, max_row = self._index_map_to_row_range(index_map)

        if min_row is not None and max_row is not None:
            self.layoutAboutToBeChanged.emit()

            index_tl = self.createIndex(min_row, 0)
            index_br = self.createIndex(max_row, self.columns() - 1)
            self.dataChanged.emit(index_tl, index_br)

            self.layoutChanged.emit()

    def get_record(self, row):
        return self._log_storage.get(self._viewed_indices[row])


class TextDelegate(qt.QtWidgets.QStyledItemDelegate):
    def __init__(self):
        super().__init__()

    def _handle_option_index(self, option, index):
        opt = qt.QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        return opt

    def _create_doc(self, opt):
        doc = qt.Qt.QTextDocument()
        doc.setPlainText(opt.text)
        return doc

    def paint(self, painter, option, index):
        opt = self._handle_option_index(option, index)

        painter.save()
        doc = self._create_doc(opt)
        opt.text = ""

        style = opt.widget.style() if opt.widget else qt.QtWidgets.QApplication.style()
        style.drawControl(qt.Qt.QStyle.CE_ItemViewItem, opt, painter)
        painter.translate(opt.rect.left(), opt.rect.top())
        clip = qt.Qt.QRectF(0, 0, opt.rect.width(), opt.rect.height())
        doc.drawContents(painter, clip)
        painter.restore()

    def sizeHint(self, option, index):
        opt = self._handle_option_index(option, index)
        doc = self._create_doc(opt)
        doc.setTextWidth(opt.rect.width())
        return qt.Qt.QSize(int(doc.idealWidth()), int(doc.size().height()))


class LogTable(qt.QtWidgets.QTableView):
    def __init__(self, parent, model):
        super().__init__(parent)
        self._model = model
        self.setModel(model)
        self._font_metrics = qt.QtGui.QFontMetrics(self.font())
        self._default_row_height = self._font_metrics.height() + 2
        self._initialize_ui()

    def _initialize_ui(self):
        # self.setWordWrap(False)
        self.setWordWrap(True)
        self.setTextElideMode(qt.QtCore.Qt.ElideMiddle)

        delegate = TextDelegate()
        self.setItemDelegate(delegate)

        h = self.horizontalHeader()
        h.setVisible(True)
        h.setSectionResizeMode(
            self._model.columns() - 1, qt.QtWidgets.QHeaderView.Stretch
        )
        h.setSectionsClickable(False)
        h.setSectionsMovable(False)

        v = self.verticalHeader()
        v.setVisible(False)
        # v.setSectionResizeMode(qt.QtWidgets.QHeaderView.Fixed)
        # v.setSectionResizeMode(qt.QtWidgets.QHeaderView.Stretch)
        v.setSectionResizeMode(qt.QtWidgets.QHeaderView.ResizeToContents)
        v.setDefaultSectionSize(self._default_row_height)
        v.setSectionsClickable(False)
        v.setSectionsMovable(False)

        self.setVerticalScrollMode(qt.QtWidgets.QAbstractItemView.ScrollPerPixel)

        self.setSelectionBehavior(qt.QtWidgets.QTableView.SelectRows)
        self._model.dataChanged.connect(self.on_data_changed)
        self._model.layoutChanged.connect(self.on_layout_changed)

    def on_data_changed(self, top_left, bottom_right, roles):
        # top_row = top_left.row()
        # bottom_row = bottom_right.row()
        # for j in range(top_row, bottom_row + 1):
        #     self.resizeRowToContents(j)
        pass

    def on_layout_changed(self, parents, hint):
        self.scrollToBottom()

    def set_column_widths(self, widths):
        h = self.horizontalHeader()
        for j, w in enumerate(widths):
            h.resizeSection(j, w)

    def get_column_widths(self):
        h = self.horizontalHeader()
        return [h.sectionSize(j) for j in range(self._model.columns() - 1)]


# logtable.py ends here
