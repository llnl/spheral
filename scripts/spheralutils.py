#!/usr/bin/env python3

import sys

if sys.version_info < (3, 7):
    raise Exception("Must have Python 3.7 or higher")

import subprocess

# Helper function for executing commands stolen from uberenv
def sexe(cmd,ret_output=False,echo=True):
    """ Helper for executing shell commands. """
    if echo:
        print("[exe: {0}]".format(cmd))

    # If we want to return the output as string a print to stdout
    # in real-time we need to let subprocess print as normal to
    # PIPE and STDOUT. We then need to read it back ourselves and
    # append to an ouput string of our own making. There is no way
    # to do this with subprocess currently.
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             encoding='utf8')
        out = "";
        while True:
            realtime_output = p.stdout.readline()

            if realtime_output == '' and p.poll() is not None:
                break

            if realtime_output:
                print(realtime_output.strip(), flush=True)
                out += realtime_output

        if echo:
            print(out)
        return out
    else:
        # If we do not need to return the output as a string, run()
        # will suffice.
        p = subprocess.run(cmd, shell=True,
                           check=True, text=True)
