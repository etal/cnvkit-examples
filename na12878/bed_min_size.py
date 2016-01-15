#!/usr/bin/env python
import sys
import pandas as pd
d = pd.read_table(sys.argv[1], names=['chrom', 'start', 'end'], header=None)
d = d[(d.end - d.start) > int(sys.argv[2])]
d.to_csv(sep='\t', header=False, index=False)
