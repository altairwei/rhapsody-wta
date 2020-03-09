from toil.job import Job
from toil.common import Toil
import subprocess
import os


def initialize_jobs(job):
    job.fileStore.logToMaster('initialize_jobs')


def runQC(job, cwl_file, cwl_filename, yml_file, yml_filename, outputs_dir, output_num):
    job.fileStore.logToMaster("runQC")
    tempDir = job.fileStore.getLocalTempDir()

    cwl = job.fileStore.readGlobalFile(cwl_file, userPath=os.path.join(tempDir, cwl_filename))
    yml = job.fileStore.readGlobalFile(yml_file, userPath=os.path.join(tempDir, yml_filename))

    subprocess.check_call(["cwltoil", cwl, yml])

    output_filename = "output.txt"
    output_file = job.fileStore.writeGlobalFile(output_filename)
    job.fileStore.readGlobalFile(output_file, userPath=os.path.join(outputs_dir, "sample_" + output_num + "_" + output_filename))
    return output_file


if __name__ == "__main__":
    options = Job.Runner.getDefaultOptions("results/rhapsody-wta-job-store")
    options.clean = "never"
    options.workDir = "tmp"
    options.outdir = "results"
    options.writeLogs = "logs"
    options.logFile = "cwltoil.log"
    options.logLevel = "INFO"
    options.retryCount = "2"
    options.maxLogFileSize = "20000000000"

    
    with Toil(options) as toil:
        # specify the folder where the cwl and yml files live
        inputs_dir = os.path.dirname(os.path.abspath(__file__))
        # specify where you wish the outputs to be written
        outputs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

        job0 = Job.wrapJobFn(initialize_jobs)

        cwl_filename = "rhapsody-wta-yaml.cwl"
        cwl_file = toil.importFile("file://" + os.path.abspath(os.path.join(inputs_dir, cwl_filename)))

        # add list of yml config inputs here or import and construct from file
        yml_filename =["template_wta.yml"
        yml_file = toil.importFile("file://" + os.path.abspath(os.path.join(inputs_dir, yml_filename)))

        job = Job.wrapJobFn(runQC, cwl_file, cwl_filename, yml_file, yml_filename, outputs_dir, output_num=str(i))
        job0.addChild(job)

        toil.start(job0)