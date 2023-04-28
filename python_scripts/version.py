#/usr/bin/python3

import os
import subprocess

coding='utf-8'

def analysis_version_string():
    """Return a string representing the version of the analysis suite."""
    # change to the directory of this source file,
    # assuming it is in the git repository
    old_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    try:
        version_string = subprocess.check_output(["git", "describe"],encoding=coding).rstrip('\n')
    except subprocess.CalledProcessError:
        version_string = 'SMASH-analysis-unknown version'
    os.chdir(old_dir)
    return version_string

def analysis_version_number():
    """Return a string representing the version of the analysis suite."""
    # change to the directory of this source file,
    # assuming it is in the git repository
    old_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    try:
        version_string = subprocess.check_output(["git", "describe"],encoding=coding).split("-")[2].rstrip('\n')
    except subprocess.CalledProcessError:
        version_string = '?.?'
    os.chdir(old_dir)
    return version_string

print(analysis_version_number())
