#!/usr/bin/python
#protexam.py
'''
This is the primary file for ProtExAM.
It is intended to call all submethods as modules.

ProtExAM is a system for extracting a set of protein mentions and
contexts from biomedical literature, such that the mentions correspond
to a given topic (e.g., the heart, lung cancer), are as comprehensive
as the literature will allow, and are context-sentitive. The last of
these properties allows for protein mentions to be filtered by type,
e.g., those observed to be expressed vs. negative mentions.

'''
__author__= "Harry Caufield"
__email__ = "jcaufield@mednet.ucla.edu"

import sys, argparse

import protexam_helpers as p_help
import protexam_query as p_query
import protexam_input as p_input
import protexam_process as p_proc
import protexam_output as p_output
import protexam_settings as p_settings

## Constants and Options
parser = argparse.ArgumentParser()
parser.add_argument("--query", help="Search for documents matching a query, in quotes."
                                     " This will be passed to PubMed so please use"
                                     " PubMed search options, including MeSH terms.", 
					action="append")
parser.add_argument("--query_file", help="Search for documents matching a query," 
                                         " starting with a text file containing one"
                                         " search term per line.",
					action="append")
args = parser.parse_args()

## Classes

## Functions

## Main
def main():
	
	print("Done.")

if __name__ == "__main__":
	sys.exit(main())

