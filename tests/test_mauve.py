
"""
test by building with: 
python setup.py build_ext --inplace
"""

import os
import sys
LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

mauvepy_path = os.path.dirname(LOCAL_DIR)
sys.path.append(mauvepy_path)

import mauve

GENOME = os.path.join(LOCAL_DIR, 'mds42_recoded.fa')
REF_GENOME = os.path.join(LOCAL_DIR, 'mds42_full.fa')

print("Building lookup table")
index_lookup_table = mauve.buildIndex(GENOME, REF_GENOME)
print("Finished")