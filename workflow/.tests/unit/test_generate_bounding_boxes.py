import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_generate_bounding_boxes():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/generate_bounding_boxes/data")
        expected_path = PurePosixPath(".tests/unit/generate_bounding_boxes/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("../workflow/results/bounding_boxes/ISPY2-115638_0.npy ../workflow/results/bounding_boxes/ISPY2-115638_1.npy ../workflow/results/bounding_boxes/ISPY2-115638_2.npy ../workflow/results/bounding_boxes/ISPY2-115638_3.npy", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "generate_bounding_boxes",
            "-f", 
            "-j1",
            "--keep-target-files",
            "--use-conda",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
