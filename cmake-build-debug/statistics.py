#!/usr/bin/python

import os
import sys
import pandas as pd

fn = sys.argv[1]
if os.path.exists(fn):
    with open(fn) as fp:
        Lines = fp.readlines()

    values = [float(i) for i in Lines]
    s = pd.Series(values)
    print(s.describe())

