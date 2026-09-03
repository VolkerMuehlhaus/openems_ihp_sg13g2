#!/usr/bin/env python3
# Keep more_examples/local_modules_copy/modules/ in sync with workflow/modules/, the canonical
# source. That folder is a real (not symlinked) local copy - GitHub's "Download ZIP" doesn't
# reliably preserve symlinks across platforms - demonstrating the local "modules" folder
# distribution method for users who don't want the gds2openEMS PyPI package. Run this whenever
# workflow/modules/ changes.

import os
import shutil

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO_ROOT, 'workflow', 'modules')
DST = os.path.join(REPO_ROOT, 'more_examples', 'local_modules_copy', 'modules')

def main():
    if os.path.isdir(DST):
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST, ignore=shutil.ignore_patterns('__pycache__'))
    print(f'Synced {DST} from {SRC}')

if __name__ == '__main__':
    main()
