#!/bin/bash
# Build the AWS cell package. Run from repo root. Flattens core+cell into one tarball,
# because the userdata script extracts everything into a single /work directory.
set -e
cd "$(dirname "$0")/.."
tar -czf /tmp/dmlab_pkg.tgz dmlab.py nbody_3d.py nbody_dm_ic.py nfw_anisotropic_ic.py sidm.py feedback.py -C fleet cal_cell.py
echo "built /tmp/dmlab_pkg.tgz ($(tar -tzf /tmp/dmlab_pkg.tgz | wc -l | tr -d ' ') files)"
tar -tzf /tmp/dmlab_pkg.tgz | sed 's/^/  /'
