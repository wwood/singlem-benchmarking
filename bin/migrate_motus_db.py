#!/usr/bin/env python3

import os
import shutil
import motus.downloadDB
import sys

db_path = motus.downloadDB.relative_path + 'db_mOTU'
destination = sys.argv[1]
if os.path.exists(destination):
    print("Destination already exists. Deleting...")
    # delete it
    shutil.rmtree(destination)
    print("Deleted", destination)
print("Copying", db_path, "to", destination)
shutil.copytree(db_path, destination)
print("done")
