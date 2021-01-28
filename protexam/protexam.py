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

import protexam_helpers as phlp
import protexam_query as pqry
import protexam_process as ppro
import protexam_output as pout
import protexam_settings as pset

## Constants and Options
parser = argparse.ArgumentParser()
query_group = parser.add_mutually_exclusive_group()
prot_group = parser.add_mutually_exclusive_group()
query_group.add_argument("--query", help="Search for documents matching a query, in quotes."
                                     " This will be passed to PubMed so please use"
                                     " PubMed search options, including MeSH terms.", 
					action="append")
query_group.add_argument("--query_file", help="Search for documents matching a query," 
                                         " starting with the name of a text file"
                                         " containing one search term per line."
                                         " By default, this assumes an OR"
                                         " relationship between all terms.",
					action="append")
prot_group.add_argument("--get_protein_entries",
                    help = "Retrieve full UniProtKB entries for a list of"
                           " UniProt accession codes, provided in a specified"
                           " filename with one code per line."
                           " Output is written to a file in its own folder,"
                           " bypassing all other steps.",
 )
prot_group.add_argument("--get_protein_aliases", 
                    help = "Retrieve UniProtKB entries for a list of"
                           " UniProt accession codes, provided in a specified"
                           " filename with one code per line."
                           " Output is written to a file in its own folder,"
                           " bypassing all other steps and limited to"
                           " alternate protein names and identifiers.",
 )
parser.add_argument("--auto", help="Run in automatic mode, accepting all options"
                                   " with a Yes.", 
					action="store_true")
args = parser.parse_args()

## Classes

## Functions

## Main
def main():
 
 print("** ProtExAM **")
 
 #A quick version check
 if sys.version_info[0] < 3:
		sys.exit("Not compatible with Python2 -- sorry!\n"
					"Exiting...")
 
 have_query = False
 have_prot_query = False
 
 if args.get_protein_entries:
  have_prot_query = True
  query_list = []
  with open(args.get_protein_entries) as query_file:
			for query_item in query_file:
				query_list.append(query_item.rstrip())
  pqry.download_uniprot_entries(query_list, "full")
  sys.exit("Exiting...")
  
 if args.get_protein_aliases:
  have_prot_query = True
  query_list = []
  with open(args.get_protein_aliases) as query_file:
			for query_item in query_file:
				query_list.append(query_item.rstrip())
  pqry.download_uniprot_entries(query_list, "alias")
  sys.exit("Exiting...")
 
 if args.query:
  have_query = True
  query = (args.query)[0]
  pmid_list, query_dir_path, webenv = pqry.run_pubmed_query(query)
  log_file_path = phlp.write_log(query, query_dir_path)
  
 if args.query_file:
  have_query = True
  query_list = []
  with open(args.query_file) as query_file:
			for query_item in query_file:
				query_list.append(query_item.rstrip())
  pmid_list, query_dir_path, webenv = pqry.run_pubmed_query(query_list)
  log_file_path = phlp.write_log("\t".join(query_list), query_dir_path)
  
 if not have_query:
  sys.exit("No query provided. Please use the --query or --query_file "
           "options.\n"
  "Exiting...")
 
 print("Saving records and logs to %s" % (query_dir_path))
 print("See %s for query details." % (log_file_path))
 
 #Prompt to continue - the function here is just an input parser
 question = "Continue with document download? (Y/N) "
 if not args.auto:
  response = phlp.get_input(question, "truefalse")
  if not response:
   sys.exit("OK, exiting...")
 else:
  print("%s Y" % (question))
 
 recs = pqry.download_pubmed_entries(pmid_list, query_dir_path, webenv)
 pqry.download_pmc_entries(recs, query_dir_path, webenv)
 
 print("Document download complete.")
 
 #Prompt to continue again
 question = "Continue with annotation download? (Y/N) "
 if not args.auto:
  response = phlp.get_input(question, "truefalse")
  if not response:
   sys.exit("OK, exiting...")
 else:
  print("%s Y" % (question))
 
 pqry.download_ptc_gene_annotations(pmid_list, query_dir_path)

 print("Annotation download complete.")
 
if __name__ == "__main__":
	sys.exit(main())

