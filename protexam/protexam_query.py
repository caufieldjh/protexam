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

from Bio import Entrez

import protexam_settings as pset

## Constants

QUERY_PATH = Path('../queries')

PERSONAL_EMAIL = pset.PERSONAL_EMAIL

## Functions

def run_pubmed_query(query):
 ''' 
 Runs a query on PubMed, with assistance from the BioPython Entrez
 module. The input may be a string or a list. If it's the latter,
 it will be combined into an appropriate search query with OR between
 all terms. The output is a list of matching PubMed IDs.
 This function also saves records (MEDLINE format with abstracts) and
 PubMed Central full text records - where available - to the queries
 folder.
 '''
 
 pmid_list = []
 
 Entrez.email = PERSONAL_EMAIL
 
 now = datetime.datetime.now()
 nowstring = now.strftime("%Y-%m-%d %H:%M:%S")
 
 #Create the queries folder if it does not yet exist
 QUERY_PATH.mkdir(exist_ok=True)
 
 #If the query is a list we need to parse it
 if isinstance(query, list):
  flat_query = " OR ".join(query)
  query = flat_query
 
 #Filename and directory setup
 query_dir_name = (query[0:20] + nowstring).replace(" ", "_").replace(":", "-")
 query_dir_path = QUERY_PATH / query_dir_name
 query_dir_path.mkdir()
 
 query_list_fn = "pmid_list_" + query_dir_name + ".txt"
 query_list_path = query_dir_path / query_list_fn
 query_absts_fn = "absts_" + query_dir_name
 query_absts_path = query_dir_path / query_absts_fn
 query_full_fn = "fulltexts_" + query_dir_name
 query_full_path = query_dir_path / query_full_fn
 
 #Now we're ready to search and retrieve
 print("Searching PubMed with the query: %s" % (query))
 print("Query date and time: %s" % (nowstring))
 
 #Use the history server by default
 
 retmax = 10000
 handle = Entrez.esearch(db="pubmed", term=query, retmax = retmax, usehistory="y")
 record = Entrez.read(handle)
 handle.close()
 
 for pmid in record['IdList']:
  pmid_list.append(pmid)
 webenv = record['WebEnv']
 count = int(record['Count'])
 
 #If we hit retmax, then it's time to iterate through the rest
 
 if count > retmax:
  pbar = tqdm(total = count, unit=" PMIDs retrieved")
  index = retmax - 1
  while index < count + retmax - 1: #Need to add retmax to get the last chunk
   handle = Entrez.esearch(db="pubmed", term=query, retmax = retmax, usehistory="y", webenv = webenv, retstart = index)
   record = Entrez.read(handle)
   handle.close()
   for pmid in record['IdList']:
    pmid_list.append(pmid)
    pbar.update(len(record['IdList']))
   index = index + retmax
  pbar.close()

 print("Query returned %s PubMed records." % (len(pmid_list)))
 
 with open(query_list_path, "w") as pmid_list_file:
  for pmid in pmid_list:
   pmid_list_file.write(pmid + "\n")
 
 print("Wrote PMID list to %s." % (query_list_path))
 
 #Next up - save records for each in MEDLINE format
 
 #Then get and save PMC entries
 
 return pmid_list
 
