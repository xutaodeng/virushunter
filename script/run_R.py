#!/usr/bin/env python

import os
import sys

f=open(sys.argv[1]) #command sh

for line in f:
	os.system(line.strip())
f.close()