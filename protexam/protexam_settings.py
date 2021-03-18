#!/usr/bin/python
#protexam_settings.py

'''
This is the ProtExAM settings file.
It sets default environment variables.
Because these are environment variables, they may vary by environment.
It's all right there in the name. Fun, right?
'''

from environs import Env
from pathlib import Path

env = Env()

PERSONAL_EMAIL = "jcaufield@mednet.ucla.edu"

QUERY_PATH = Path('../queries')

