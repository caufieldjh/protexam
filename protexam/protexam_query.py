#!/usr/bin/python
#protexam_start.py
'''
Functions for ProtExAM to perform queries on knowledgebases,
literature resources (e.g., PubMed), and annotation systems. 
'''

import os, random, sys, datetime
import urllib

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
 all terms. The output is a list of matching PubMed IDs (also saved
 as a file), the path to the query directory, and the NCBI WebEnv
 address, so we can continue to use the history server.
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
 query_full_fn = "fulltexts_" + query_dir_name
 query_full_path = query_dir_path / query_full_fn
 
 #Now we're ready to search and retrieve
 print("Searching PubMed with the query: %s" % (query))
 print("Query date and time: %s" % (nowstring))
 
 #Use the history server by default
 
 retmax = 1000
 handle = Entrez.esearch(db="pubmed", term=query, retmax = retmax, usehistory="y")
 record = Entrez.read(handle)
 handle.close()
 
 for pmid in record['IdList']:
  pmid_list.append(pmid)
 webenv = record['WebEnv']
 count = int(record['Count'])
 
 #If we hit retmax, then it's time to iterate through the rest
 
 try:
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
 except urllib.error.HTTPError as e:
  print(e)

 print("Query returned %s PubMed records." % (len(pmid_list)))
 
 with open(query_list_path, "w") as pmid_list_file:
  for pmid in pmid_list:
   pmid_list_file.write(pmid + "\n")
 
 print("Wrote PMID list to %s." % (query_list_path))
 
 return pmid_list, query_dir_path, webenv
 
def download_pubmed_entries(pmid_list, query_dir_path, webenv):
 '''
 Retrieve contents of PubMed entries with IDs in the provided list.
 Input is a list of PMIDs, the path to a previously created query
 directory, and a previously created WebEnv address.
 Output, the set of PubMed records, is returned and is saved to 
 the same directory created for the query.
 '''
 
 from Bio import Medline
 
 query_entries_fn = "entries.txt"
 query_entries_path = query_dir_path / query_entries_fn
 
 count = len(pmid_list)
 
 print("Retrieving contents for %s PubMed entries." % (count))
 
 index = 0
 retmax = 1000
 pm_recs = []
 pbar = tqdm(total = count, unit=" entries retrieved")
 
 try:
  if count > retmax:
   while index < count + retmax - 1: #Need to add retmax to get the last chunk
    handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline", retmode = "text", retmax = retmax, usehistory="y", webenv = webenv, retstart = index)
    these_pm_recs = list(Medline.parse(handle))
    handle.close()
    new_pm_recs = pm_recs + these_pm_recs
    pm_recs = new_pm_recs
    pbar.update(len(these_pm_recs))
    index = index + retmax
  else:
   handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline", retmode = "text", retmax = retmax, usehistory="y", webenv = webenv, retstart = index)
   pm_recs = list(Medline.parse(handle))
   handle.close()
   pbar.update(len(pm_recs))
  pbar.close()
 except urllib.error.HTTPError as e:
  print(e)
  
 #ERROR - the download process throws a 400 error when finishing,
 #so I guess it's running one more time than necessary.
 #The one below does it too.
  
 print("Retrieved %s entries." % (len(pm_recs)))
 
 with open(query_entries_path, "w", encoding="utf-8") as entryfile:
  for rec in pm_recs:
   entryfile.write(str(rec) + "\n")
 
 print("Wrote entries to %s." % (query_entries_path))
 
 return pm_recs

def download_pmc_entries(pm_recs, query_dir_path, webenv):
 '''
 Retrieve contents of PubMed Central full text entries with IDs in the 
 provided list.
 Input is a set of PubMed records as returned from the 
 download_pubmed_entries function, the path to a previously created 
 query directory, and a previously created WebEnv address.
 Output is saved to the same directory created for the query.
 '''
 
 from Bio import Medline
 import xml.dom.minidom
 import xml.etree.ElementTree as ET
 
 query_entries_fn = "pmc_fulltexts.txt"
 query_entries_path = query_dir_path / query_entries_fn
 
 pmc_ids = []
 for rec in pm_recs:
  if 'PMC' in rec.keys():
   pmc_ids.append(rec['PMC'])
  
 count = len(pmc_ids)
 
 print("Retrieving contents for %s PubMed Central texts." % (count))
 
 #Need to trim off the "PMC" from each ID
 pmc_ids_trim = [id[3:] for id in pmc_ids]
 
 index = 0
 retmax = 1000
 pm_recs = []
 pbar = tqdm(total = count, unit=" entries retrieved")
 
 try:
  if count > retmax:
   while index < count + retmax - 1: #Need to add retmax to get the last chunk
    handle = Entrez.efetch(db="pmc", id=pmc_ids_trim, rettype="full", retmode = "xml", retmax = retmax, usehistory="y", webenv = webenv, retstart = index)
    hxml = xml.dom.minidom.parseString(handle.read())
    hxml_pretty = hxml.toprettyxml()
    with open(query_entries_path, "a", encoding="utf-8") as outfile:
     for line in hxml_pretty:
      outfile.write(line)
    handle.close()
    pbar.update(retmax) #Not quite right 
    index = index + retmax
  else:
   handle = Entrez.efetch(db="pmc", id=pmc_ids_trim, rettype="full", retmode = "xml", retmax = retmax, usehistory="y", webenv = webenv, retstart = index)
   hxml = xml.dom.minidom.parseString(handle.read())
   hxml_pretty = hxml.toprettyxml()
   with open(query_entries_path, "w", encoding="utf-8") as outfile:
     for line in hxml_pretty:
      outfile.write(line)
   handle.close()
   pbar.update(count)
  pbar.close()
 except urllib.error.HTTPError as e:
  print(e)
 
 #ERROR - Writing multiple XML docs to the same file means they need merging.
 #Still need to do this.
 
 print("Wrote entries to %s." % (query_entries_path))
 
 #Need to check PMC results as some may not have been available
 print("Checking on documents...")
 
 pub_ids = []
 fulldoc_ids =[]
 
 try:
  tree = ET.parse(query_entries_path)
  root = tree.getroot()
  for article in root.findall("./article"):
   for article_id in article.findall("front/article-meta/article-id"):
     if article_id.attrib["pub-id-type"] == "pmid":
      pmid = article_id.text
      pub_ids.append(pmid)
   body = article.find("body")
   if body is not None:
    fulldoc_ids.append(pmid)
 except xml.etree.ElementTree.ParseError as e:
  print("Encountered this error in parsing the PMC output: %s" % (e))
  print("There may be a problem with the tree structure (e.g., two roots).")
   
 print("Retrieved %s PMC entries." % (len(pub_ids)))
 print("PMC entries with full body text: %s" % (len(fulldoc_ids)))
