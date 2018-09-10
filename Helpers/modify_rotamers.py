import sys
import re
import logging

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def change_angle(filename, new_angle):
    with open(filename) as rotamers_file:
        rotamers = rotamers_file.read()
    rotamers_replaced = re.sub("FREE\d+", "FREE{}".format(new_angle), rotamers)
    with open(filename, "w") as rot_to_write:
        rot_to_write.write(rotamers_replaced)
