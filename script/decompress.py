#!/usr/bin/env python3

import os
import pandas as pd
import tempfile
import shutil
import sys

directory = sys.argv[1]

for filename in os.listdir(directory):
    if not filename.endswith(".csv"):
        continue

    file_path = os.path.join(directory, filename)

    with open(file_path, "r") as f:
        pmin, pmax, dcompress = map(float, f.readline().split())

    df = pd.read_csv(
        file_path,
        skiprows=1,
        sep=r"\s+",
        dtype={"pressure": float}
    )

    df["pressure"] = df["pressure"] * dcompress + pmin

    with tempfile.NamedTemporaryFile("w", delete=False, dir=directory) as tmp:
        tmp_path = tmp.name
        df.to_csv(
            tmp_path,
            index=False,
            sep= " ",
            float_format="%.6f"
        )

    shutil.move(tmp_path, file_path)