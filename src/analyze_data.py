#!/usr/bin/env python

import csv

reader = csv.DictReader(open("sim_kazan_05_12_2013_fec_3.csv"),delimiter=";")

first = True

for row in reader:

    if first:
        print ";".join(row.keys())
        first = False

    if row["type"] == "1":
        print ";".join(row.values())



