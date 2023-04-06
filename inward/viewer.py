#!/usr/bin/env python

import argparse

import polyscope as ps
import pymem3dg.visual as dg_vis
from pathlib import Path

import driver


parser = argparse.ArgumentParser(description="Visualize.")
parser.add_argument("filename", action="store", help="File to visualize")

args = parser.parse_args()

file = Path(args.filename)

if file.suffix == ".nc":
    dg_vis.animate(
        str(file), parameters=driver.p, geodesicDistance=True, externalForce=True
    )
elif file.suffix == ".ply":
    dg_vis.visualizePly(str(file), "lineCapillaryForce", "externalForce")
else:
    raise RuntimeError("Unknown file type")
ps.show()
