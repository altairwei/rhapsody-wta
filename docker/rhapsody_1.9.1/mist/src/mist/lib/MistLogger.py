import os
import sys
import traceback
import functools
import contextlib
# merge the logging module namespace with our own to enable as a drop in replacement
from logging import *


node_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
logger = getLogger('')  # grab the root logger
logger.setLevel(DEBUG)

# Log to standard error
logger.addHandler(StreamHandler(stream=sys.stderr))

# Log to file
try:
    log_file_handler = FileHandler('{}.log'.format(node_name))
except IOError:  # jenkins doesn't deal well with logging to a file
    pass
else:
    ff = Formatter('%(asctime)s - {node_name} - %(name)s - %(levelname)s - %(message)s'.format(node_name=node_name))
    log_file_handler.setFormatter(ff)
    logger.addHandler(log_file_handler)

# Quiet matplotlib chatter
getLogger("matplotlib").setLevel(WARNING)


def log_death(f):
    """capture and log exceptions, then die"""
    @functools.wraps(f)
    def log_death_inner(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except:
            e = traceback.format_exc()
            error(e)
            sys.exit(1)
    return log_death_inner


@contextlib.contextmanager
def log_bookend(message):
    import time
    info("{}...".format(message))
    start = time.time()
    yield
    end = time.time()
    info("...done. (That took {:.2f} seconds.)".format(end - start))
