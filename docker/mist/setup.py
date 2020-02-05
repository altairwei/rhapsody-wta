from __future__ import absolute_import
from setuptools import find_packages
from codecs import open
from os import path
from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import os

import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))
import mist.apps._version as v

cython_flag = "CYTHON_ENABLED_FOR_MIST_PIPELINE"
cython_enabled = (
    False if cython_flag not in os.environ else bool(int(os.environ[cython_flag]))
)

with open("README.md", encoding="utf-8") as readme, \
        open(path.join(path.dirname(__file__), "dependency_resolution", "requirements.txt")) as req, \
        open(path.join(path.dirname(__file__), "dependency_resolution", "requirements-test.txt")) as req_test:
    long_description = readme.read()
    required_packages = req.read()
    required_packages_test = req_test.read()


setup(
    name="mist",
    version=v.__version__,
    install_requires=required_packages,
    tests_require=required_packages_test,
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    entry_points={
        "console_scripts": [
            "mist_add_to_bam.py = mist.apps.AddtoBam:main",
            "mist_quality_filter.py = mist.apps.QualityFilter:main",
            "mist_check_references.py = mist.apps.CheckReferences:main",
            "mist_annotate_R1.py = mist.apps.AnnotateR1:package_main",
            "mist_annotate_R2.py = mist.apps.AnnotateR2:main",
            "mist_annotate_reads.py = mist.apps.AnnotateReads:main",
            "mist_annotate_molecules.py = mist.apps.AnnotateMolecules:main",
            "mist_get_datatables.py = mist.apps.GetDataTables:main",
            "mist_metrics.py = mist.apps.Metrics:main",
            "mist_cluster_analysis.py = mist.apps.ResolveClusterAnalysis:main",
            "mist_check_fastqs.py = mist.apps.CheckFastqs:main",
            "mist_split_fastq.py = mist.apps.SplitFastq:main",
            "mist_sparse_to_dense.py = mist.apps.sparse_to_dense_datatable:main",
        ]
    },
    ext_modules=cythonize(["src/mist/**/*.py"]) if cython_enabled else None,
    cmdclass={"build_ext": build_ext} if cython_enabled else {},
)
