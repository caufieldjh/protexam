#!/usr/bin/python
#protexam_process.py
'''
Processing functions for ProtExAM.
'''

import scispacy
import spacy

__spacy_nlp = spacy.load("en_ner_bionlp13cg_md")

def spacy_gene_ner(text):
    '''Takes a string and returns a list of tuples denoting genes 
found. Each tuple in the list has two elements. The first is the
 substring identified as a gene. The second is a  tuple with the
 starting index of the substring and the substring's length, in that
 order.'''
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

if __name__ == "__main__":
    print("Starting tests for functions in protexam_process")
    text = "Tafazzin is a protein that in humans is encoded by the TAZ gene. It is important to note that the TAZ gene is frequently confused with a protein called TAZ (transcriptional coactivator with PDZ-binding motif, a 50kDA protein). "
    genes = spacy_gene_ner(text)
    print("Starting tests for spacy_gene_ner")
    for result in genes:
        assert(result[0] == text[result[1][0]:result[1][0]+result[1][1]])
    print("Tests for spacy_gener_ner passed")
    print("All tests passed")
