#!/usr/bin/env python

import os
import sys
import multiprocessing

cpu=multiprocessing.cpu_count()

mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
mem_gib = int(mem_bytes/(1024.**3) )

print sys.argv[1], cpu, mem_gib #server and number of CPU