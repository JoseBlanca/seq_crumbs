# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of seq_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with seq_crumbs. If not, see <http://www.gnu.org/licenses/>.
from bisect import bisect


class OrderedSet(object):

    def __init__(self, ordered_set=None):
        if ordered_set is None:
            ordered_set = []
        self.list = ordered_set

    def check_add(self, item):
        index = bisect(self.list, item)
        if index == 0:
            self.list.insert(0, item)
            return True
        else:
            present = True if item == self.list[index - 1] else False
            if present:
                return False
            else:
                self.list.insert(index, item)
                return True

    def __contains__(self, item):
        index = bisect(self.list, item)
        if index == 0:
            return False
        return True if item == self.list[index - 1] else False

    def __len__(self):
        return len(self.list)
