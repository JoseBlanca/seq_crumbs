#!/usr/bin/env python

# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of bam_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with bam_crumbs. If not, see <http://www.gnu.org/licenses/>.

import pysam


def sam_to_bam(sam_fhand, bam_fhand):
    'It converts the file from sam to bam'
    for line in  pysam.view('-hSb', sam_fhand.name):
        bam_fhand.write(line)
    bam_fhand.flush()


def bam_to_sam(bam_fhand, sam_fhand):
    'It converts the file from sam to bam'
    for line in pysam.view('-h', bam_fhand.name):
        sam_fhand.write(line)
    sam_fhand.flush()







