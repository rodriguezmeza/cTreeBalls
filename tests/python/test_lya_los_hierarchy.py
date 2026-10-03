#!/usr/bin/env python3
"""LOS hierarchy oracles, MPI failures and recovery; checks remain active under -O."""
import sys
from test_lya_hierarchy import main

if __name__ == "__main__":
    if "--los-tree" not in sys.argv:
        sys.argv.append("--los-tree")
    main()
