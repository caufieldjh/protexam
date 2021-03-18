#!/usr/bin/python
#protexam_process.py
'''
Processing functions for ProtExAM.
'''

import os, json, re
from pathlib import Path

import protexam_settings as pset

## Constants

QUERY_PATH = pset.QUERY_PATH

## Functions

def spacy_gene_ner(text):
	'''Takes a string and returns a list of tuples denoting genes 
	found. Each tuple in the list has two elements. The first is the
	substring identified as a gene. The second is a  tuple with the
	starting index of the substring and the substring's length, in that
	order.'''
 
	import scispacy
	import spacy
	
	__spacy_nlp = spacy.load("en_ner_bionlp13cg_md")
	
	genes_found = []
	doc = __spacy_nlp(text)
	for ent in doc.ents:
		if ent.label_ == "GENE_OR_GENE_PRODUCT":
			gene = ent.text
			str_start = ent.start_char
			str_length = len(gene)
			entry = (gene, (str_start, str_length))
			genes_found.append(entry)
	return genes_found

def extract_full_text_json():
	'''If entries for full texts have been downloaded from PubMed 
	Central and now reside in query folders, this function will produce
	a JSON file containing available document full-texts for each
	PubMed ID in the pmid_list for each query.
	The full text entries may then be used to enrich structured PubMed
	entries with text beyond abstracts.'''
	
	import xml.etree.ElementTree as ET
	
	fulltext_output_fn = "pmid_fulltext.json"
	fulltext_output_path = QUERY_PATH / fulltext_output_fn
	
	def PubmedXML(file_name):
		with open(file_name,"rb") as f:
			context = ET.iterparse(f, events=("start","end"))
			for event, elem in context:
				yield event,elem
				elem.clear()
                        
	def FullTextFromXML(file_name):
		current_vals = {"pmid":None, "pmcid": None, "text":None}
		xml = PubmedXML(file_name)
		flag = False
		for event, elem in xml:
		#found new article
			if elem.tag == "article" and event == "start":
			#yield only if its not the first article tag found
				if flag:
					yield current_vals
				flag = True
				current_vals = {"pmid":None, "pmcid": None, "text":None}
			if elem.tag == "article-id" and event == "start":
				if elem.get('pub-id-type') == 'pmid':
					pmid = elem.text
					current_vals["pmid"] = pmid
				if elem.get('pub-id-type') == 'pmc':
					pmcid = elem.text
					current_vals["pmcid"] = pmcid
			if elem.tag == "body" and event == "start":
				current_vals["text"] = ''.join(elem.itertext())
				replaced_vals = current_vals["text"].replace("\n", " ")
				replaced_vals = re.sub("\t+", "\t", replaced_vals)
				replaced_vals = re.sub(" +", " ", replaced_vals)
				current_vals["text"] = replaced_vals
				print(current_vals["text"])
            
	def pmids_from_file(file_name):
		with open(file_name,"r") as f:
			pmids = set(line.strip() for line in f.readlines())
			return pmids
	
	print("Will output PMC full texts for each PMID across all lists in the queries folders.")
	
	folders_to_ignore = set(["__pycache__",".",".."])
	folders = [folder for folder in os.scandir(QUERY_PATH) if folder.is_dir() and folder.name not in folders_to_ignore]
	pmid_fulltext_dict = {}
	for folder in folders:
		files = os.scandir(folder)
		pmid_list_path, xml_path = None, None 
		for file_ in files:
			if "fulltexts.txt" in file_.path:
				xml_path = file_.path
			if "pmid_list" in file_.path:
				pmid_list_path = file_.path
		if pmid_list_path is None:
			print(f"{folder} does not have a PMID list. Skipping this folder...")
			continue
		if xml_path is None:
			print(f"{folder} does not have any PMC xml files. Skipping this folder...")
			continue
		print(f"Working on dir:{folder}")
		full_text_dicts = FullTextFromXML(xml_path)
		pmids = pmids_from_file(pmid_list_path)
		for full_text_dict in full_text_dicts:
			if full_text_dict["pmid"] in pmids:
				pmid_fulltext_dict[full_text_dict["pmid"]] = full_text_dict["text"]
	json.dump(pmid_fulltext_dict, open(fulltext_output_path,"w"))
	print("Wrote output to " + str(fulltext_output_path))
    
if __name__ == "__main__":
	print("Starting tests for functions in protexam_process")
	text = "Tafazzin is a protein that in humans is encoded by the TAZ gene. It is important to note that the TAZ gene is frequently confused with a protein called TAZ (transcriptional coactivator with PDZ-binding motif, a 50kDA protein). "
	genes = spacy_gene_ner(text)
	print("Starting tests for spacy_gene_ner")
	for result in genes:
		assert(result[0] == text[result[1][0]:result[1][0]+result[1][1]])
	print("Tests for spacy_gener_ner passed")
	print("All tests passed")
