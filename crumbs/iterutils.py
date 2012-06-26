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

'''
Created on 2012 eka 26

@author: peio
'''
import random


def sample(iterator, length, sample_size):
    'It makes a sample for the given iterator'

    if not 0 <= sample_size <= length:
        raise ValueError("sample larger than population")

    selected = set(random.randint(0, length - 1) for num in range(sample_size))
    selected_add = selected.add
    while len(selected) < sample_size:
        selected_add(random.randint(0, length - 1))

    for index, item in enumerate(iterator):
        if index in selected:
            yield item


def length(iterator):
    'it counts the items in an iterator. It consumes the iterator'
    count = 0
    for item in iterator:
        count += 1
    return count
