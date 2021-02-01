#!/usr/bin/python
#protexam_start.py
'''
Functions for ProtExAM to perform queries on knowledgebases,
literature resources (e.g., PubMed), and annotation systems. 
'''

import os, random, requests, sys, datetime
import urllib
import xml.dom.minidom

from bs4 import BeautifulSoup

from pathlib import Path
from tqdm import *

from Bio import Entrez

import protexam_settings as pset
import protexam_helpers as phlp

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
 randnum = random.randint(0,10000000)
 
 #Create the queries folder if it does not yet exist
 QUERY_PATH.mkdir(exist_ok=True)
 
 #If the query is a list we need to parse it
 if isinstance(query, list):
  flat_query = " OR ".join(query)
  query = flat_query
 
 #Filename and directory setup
 repchars = ":[]()\\\/"
 query_dir_name = (query[0:40] + "_" + str(randnum)).replace(" ", "_")
 for char in repchars:
  query_dir_name = query_dir_name.replace(char,"-")
 query_dir_name = query_dir_name.lstrip("-")
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
   multiple_xml = True
   while index < count + retmax - 1: #Need to add retmax to get the last chunk
    handle = Entrez.efetch(db="pmc", id=pmc_ids_trim, rettype="full", retmode = "xml", retmax = retmax, usehistory="y", webenv = webenv, retstart = index)
    hxml = xml.dom.minidom.parseString(handle.read())
    hxml_pretty = hxml.toprettyxml()
    #Don't include multiple sets of headers/footers or we'll end up with multiple roots
    with open(query_entries_path, "a", encoding="utf-8") as outfile:
     outfile.write(hxml_pretty)
    handle.close()
    pbar.update(retmax) #Not quite right 
    index = index + retmax
  else:
   multiple_xml = False
   handle = Entrez.efetch(db="pmc", id=pmc_ids_trim, rettype="full", retmode = "xml", retmax = retmax, usehistory="y", webenv = webenv, retstart = index)
   hxml = xml.dom.minidom.parseString(handle.read())
   hxml_pretty = hxml.toprettyxml()
   with open(query_entries_path, "w", encoding="utf-8") as outfile:
    outfile.write(hxml_pretty)
   handle.close()
   pbar.update(count)
  pbar.close()
 except urllib.error.HTTPError as e:
  print(e)
 
 #Merge xml with multiple roots
 if multiple_xml:
  print("Merging entries...")
  with open(query_entries_path, "r+", encoding="utf-8") as outfile:
   new_outfile = outfile.readlines()
   outfile.seek(0)
   have_header = False
   for line in new_outfile:
    if "<pmc-articleset>" in line and not have_header:
     have_header = True
     outfile.write(line)
    elif "<pmc-articleset>" in line and have_header:
     pass
    elif any(noline in line for noline in ["</pmc-articleset>",
                                           "<?xml version=\"1.0\" ?>",
                                           "<!DOCTYPE pmc-articleset",
                                           "  PUBLIC '-//NLM//DTD ARTICLE SET 2.0//EN'",
                                           "  'https://dtd.nlm.nih.gov/ncbi/pmc/articleset/nlm-articleset-2.0.dtd'>"]):
     pass
    else:
      outfile.write(line)
   outfile.write("</pmc-articleset>")
   outfile.truncate()
     
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
 
def download_ptc_gene_annotations(idlist, query_dir_path):
 '''
 Retrieve annotations from PubTator Central for genes, in BioC XML 
 format, given a list of PMIDs.
 Also requires a previously created path to save annotations to.
 Annotations will include both those on title/abstract and those
 covering full texts (i.e., from PMC docs).
 Note that details of BioC format are here:
 https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/tmTools/Format.html.
 '''
 
 annotations_fn = "gene_annotations.txt"
 annotations_path = query_dir_path / annotations_fn
 
 count = len(idlist)
 batch_size = 100
 print("Retrieving PubTator Central gene annotations for %s PMIDs." % (count))
 
 pbar = tqdm(total = count, unit=" annotation sets retrieved")
 for batch in phlp.batch_this(idlist,batch_size): 
  pmids_joined = ",".join(batch)
  r = requests.get("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?&concepts=gene&pmids=" + pmids_joined)
  rxml_docs = xml.dom.minidom.parseString(r.text)
  rxml_pretty = rxml_docs.toprettyxml()
  
  with open(annotations_path, "a", encoding="utf-8") as outfile:
   outfile.write(rxml_pretty)
   
  pbar.update(len(batch))
 pbar.close()
 
 multiple_xml = False
 if batch_size < count:
  multiple_xml = True
  
 if multiple_xml:
  print("Merging annotation sets...")
  with open(annotations_path, "r+", encoding="utf-8") as outfile:
   new_outfile = outfile.readlines()
   outfile.seek(0)
   have_header = False
   for line in new_outfile:
    if "<collection>" in line and not have_header:
     have_header = True
     outfile.write(line)
    elif "<collection>" in line and have_header:
     pass
    elif any(noline in line for noline in ["</collection>",
                                           "<?xml version=\"1.0\" ?>",
                                           "<!DOCTYPE collection",
                                           "  SYSTEM 'BioC.dtd'>",
                                           "	<source>PubTator</source>",
                                           "	<date/>",
                                           "	<key>BioC.key</key>"]):
     pass
    else:
      outfile.write(line)
   outfile.write("</collection>")
   outfile.truncate()
 
 print("Wrote entries to %s." % (annotations_path))
 
def download_uniprot_entries(idlist, mode):
 '''
 Retrieve full entries for a list of UniProtKB protein accession codes.
 Also requires a mode (one of "full" or "alias", where the latter
 limits output to alternate gene names and IDs).
 Queries are written to files in a newly created directory under
 the queries folder.
 '''
 
 #Just XML output for now
 #Alias mode not functional yet
 
 import json
 import xmlschema
 import xml.dom.minidom
 import xml.etree.ElementTree as ET
 
 now = datetime.datetime.now()
 nowstring = now.strftime("%Y-%m-%d_%H_%M_%S")
 
 query_dir_name = "ProteinQuery_" + nowstring
 query_dir_path = QUERY_PATH / query_dir_name
 
 if mode == "full":
  proteins_fn = "prot_entries.txt"
 elif mode == "alias":
  proteins_fn = "aliases.txt"
 proteins_path = query_dir_path / proteins_fn
 proteins_xml_path = query_dir_path / "prot_entries.xml"
 query_dir_path.mkdir()
 
 schema = xmlschema.XMLSchema('https://www.uniprot.org/docs/uniprot.xsd')
 
 count = len(idlist)
 print("Retrieving UniProtKB entries for %s accessions." % (count))
 
 pbar = tqdm(total = count, unit=" protein entries retrieved")
 for entry in idlist:
  url = 'https://www.uniprot.org/uniprot/'+ entry + ".xml"
  with urllib.request.urlopen(url) as r:
   raw_data = r.read().strip()
   rxml_prot = xml.dom.minidom.parseString(raw_data)
   rxml_pretty = rxml_prot.toprettyxml(newl='')
  with open(proteins_xml_path, "a", encoding="utf-8") as outfile:
   outfile.write(rxml_pretty)
  pbar.update(1)
  
 pbar.close()
 
 print("Merging entries...")
 up_head = "<uniprot xmlns=\"http://uniprot.org/uniprot\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd\">"
 with open(proteins_xml_path, "r+", encoding="utf-8") as outfile:
  new_outfile = outfile.readlines()
  outfile.seek(0)
  have_header = False
  outfile.write(up_head + "\n")
  have_header = True
  for line in new_outfile:
   if up_head in line and have_header:
    pass
   elif any(noline in line for noline in ["</uniprot>"]):
    pass
   elif "copyright>" in line:
    pass
   elif "Copyrighted by the UniProt Consortium" in line:
    pass
   elif "Distributed under the Creative Commons Attribution" in line:
    pass
   else:
    outfile.write(line)
  outfile.write("\t<copyright>\n"
                "Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms\n"
                "Distributed under the Creative Commons Attribution (CC BY 4.0) License\n"
                "\t</copyright>\n"
                "</uniprot>")
  outfile.truncate()
 print("Wrote XML entries to %s." % (proteins_xml_path))

 print("Parsing XML entries...")
 tree = ET.parse(proteins_xml_path)
 entry_dict = schema.to_dict(tree)
 content = entry_dict['{http://uniprot.org/uniprot}entry']
 with open(proteins_path, "w", encoding="utf-8") as outfile:
  for entry in content:
   outfile.write(str(entry) + "\n")
 
 print("Wrote entries to %s." % (proteins_path))
