import os
import glob
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

pwd = os.getcwd()

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
gtf_file = os.path.join(pwd, snakemake.input.gtf)
bam_file = os.path.join(pwd, snakemake.input.bam)

with TemporaryDirectory() as tempdir:
    # Output of TPMCalculator can't be changed, so we change working directory
    shell(
        "(cd {tempdir} &&"
        " TPMCalculator -v -p"
        " -g {gtf_file} -b {bam_file} -k gene_id)"
        " {log}"
    )
    shutil.move(os.path.join(tempdir, Path(bam_file).stem + "_genes.out"), snakemake.output.out)
    shutil.move(os.path.join(tempdir, Path(bam_file).stem + "_genes.ent"), snakemake.output.ent)
    shutil.move(os.path.join(tempdir, Path(bam_file).stem + "_genes.uni"), snakemake.output.uni)
