import sys
import pandas as pd
import sklearn as sk
# error handling code from https://github.com/gymreklab/cse185-demo-project/blob/main/mypileup/myutils.py

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg):
        """
        Print an error message and die

        Parameters
        ----------
        msg : str
           Error message to print
        """
        sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
        sys.exit(1)