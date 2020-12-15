#!/usr/bin/python
#protexam_start.py
'''
Functions for ProtExAM to perform queries on knowledgebases,
literature resources (e.g., PubMed), and annotation systems. 
'''

import os, random, sys, datetime

from bs4 import BeautifulSoup

from pathlib import Path
from tqdm import *

## Constants

## Functions

def run_pubmed_query(query):
 ''' 
 Runs a query on PubMed, with assistance from the BioPython Entrez
 module. The input may be a string or a list. If it's the latter,
 it will be combined into an appropriate search query with OR between
 all terms. The output is a list of matching PubMed IDs.
 '''
 
 pmid_list = []
 
 now = datetime.datetime.now()
 
 #If the query is a list we need to parse it
 if isinstance(query, list):
  flat_query = " OR ".join(query)
  query = flat_query
 
 print("Searching PubMed with the query: %s" % (query))
 print("Query date and time: %s" % (now.strftime("%Y-%m-%d %H:%M:%S")))
 
 return pmid_list
