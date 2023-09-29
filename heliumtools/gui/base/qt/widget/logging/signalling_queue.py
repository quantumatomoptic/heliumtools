#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2017-12-19 12:01:13 srlab"

#  file       signalling_queue.py
#  author     Sebastian Blatt
#  created    2017-12-15 10:45:49
#
#  Copyright (C) 2011 -- 2016 Christoph Gohle, Christian Gross
#  Copyright (C) 2016 -- 2017 Christoph Gohle, Christian Gross, Sebastian Blatt
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

import collections
import math


class SignallingQueue(object):
    """A queue with a maximum length based on collections.deque that can
    be appended, extended, and cleared. On each of these iterations,
    the queue calls `on_update_callback` with an index map consisting
    of pairs of indices (old_index, new_index). The index map is also
    the return value of each of these functions.

    Here, old_index == None indicates that the item at old_index was
    deleted, and new_index == None indicates that the item now at
    new_index is new.

    The idea of this object is to let a subscriber know which items in
    its list of items have changed. For instance, we would like to
    create filtered views of SignallingQueue that can be updated
    automatically with the minimum number of resources instead of
    updating the whole view on every change.

    """

    def __init__(self, max_length=None, on_update_callback=None):
        self._deque = collections.deque(maxlen=max_length)
        self._on_update_callback = on_update_callback

    def __repr__(self):
        rc = (
            "SignallingQueue(length={}, max_length={}, " + "on_update_callback={})"
        ).format(self.length(), self.max_length(), self._on_update_callback)

        if self.length() > 0:
            n_zero = int(math.ceil(math.log10(self.length())))
            for j, x in enumerate(self._deque):
                rc += ("\n  {:0" + str(n_zero) + "d}: {}").format(j, x)
        return rc

    def length(self):
        return len(self._deque)

    def max_length(self):
        return self._deque.maxlen

    def append(self, x):
        if x is None:
            return []

        index_map = []
        n_old = self.length()
        if n_old == self.max_length():
            index_map = [(0, None)]
            index_map.extend([(j, j - 1) for j in range(1, n_old)])
            index_map.append((None, n_old - 1))
        else:
            index_map = [(None, n_old)]

        self._deque.append(x)
        if self._on_update_callback is not None:
            self._on_update_callback(index_map)
        return index_map

    def extend(self, iterable):
        if iterable is None:
            return []

        index_map = []

        N = self.max_length()
        n = self.length()
        m = len(iterable)

        if n + m <= N:
            index_map = [(None, j) for j in range(n, n + m)]
        else:
            if m >= N:
                index_map = [(j, None) for j in range(N)] + [
                    (None, j) for j in range(N)
                ]
            else:
                n_drop = n + m - N
                index_map = (
                    [(j, None) for j in range(n_drop)]
                    + [(j, j - n_drop) for j in range(n_drop, n_drop + N - m)]
                    + [(None, j) for j in range(N - n_drop, N)]
                )

        self._deque.extend(iterable)
        if self._on_update_callback is not None:
            self._on_update_callback(index_map)
        return index_map

    def clear(self):
        index_map = [(j, None) for j in range(self.length())]
        self._deque.clear()
        if self._on_update_callback is not None:
            self._on_update_callback(index_map)
        return index_map

    def filter_indices(self, selector, max_records, queue_indices=None):
        idx = (
            queue_indices
            if queue_indices is not None
            else reversed(range(self._queue.length()))
        )

        rc = []
        for j in idx:
            record = self._queue[idx]
            if self._selector(record):
                rc.append(j)
            if len(rc) > max_records:
                break
        return rc

    def get(self, index):
        if index < 0 or index >= self.length():
            raise Exception
        return self._deque[index]


if __name__ == "__main__":
    pass

# signalling_queue.py ends here
