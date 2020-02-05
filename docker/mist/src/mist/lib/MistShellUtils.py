import subprocess
import logging


def shell_command_log_stderr(cmd, **kwargs):
    """run a shell command using the subprocess.call api, but redirect output to log"""
    logging.debug('Running the following: `{}`'.format(cmd))
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        **kwargs
    )
    for line in iter(p.stdout.readline, ''):
        logging.debug(line.strip())
    p.stdout.close()
    res = p.wait()
    if res != 0:
        raise subprocess.CalledProcessError(res, cmd)