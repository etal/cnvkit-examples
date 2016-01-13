#!/usr/bin/env python

"""Append the "gc" and "rmask" columns from one reference onto another.

This lets us skip recalculating GC and RepeatMasker values from the reference
genome sequence when creating another CNVkit reference.
Both CNVkit references must have the same number of rows (corresponding to the
same positions).
"""

import argparse
import sys

import cnvlib

AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument("orig")
AP.add_argument("other")
AP.add_argument("-o", "--output", type=argparse.FileType('w'),
                default=sys.stdout)
args = AP.parse_args()

orig_arr = cnvlib.read(args.orig)
other_arr = cnvlib.read(args.other)
assert len(other_arr) == len(orig_arr)

other_arr["gc"] = orig_arr["gc"]
other_arr["rmask"] = orig_arr["rmask"]
other_arr.sort()
other_arr.sort_columns()
other_arr.write(args.output)
