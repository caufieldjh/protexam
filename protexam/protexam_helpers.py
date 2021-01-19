#!/usr/bin/python
#protexam_helpers.py
'''
General purpose helper functions for ProtExAM.
'''

from zipfile import ZipFile
import xlrd
import csv, datetime

def batch_this(mylist, n=1):
 '''Prepares batches of maximum size n from a list.'''
 length = len(mylist)
 for ndx in range(0, length, n):
  yield mylist[ndx:min(ndx + n, length)]

def decompress(filepath, outpath):
	'''Takes a Path filename of a compressed file
		and the intended output path as input.
		Decompresses to same path.
		Doesn't have a return.
		Just for ZIP compression for now but will handle all
		necessary formats.'''
	
	with ZipFile(filepath, 'r') as zip_ref:
		zip_ref.extractall(outpath)
		
def convert_xlsx_to_tsv(filepath, tabfilepath):
	'''Converts an Excel spreadsheet file (XLSX) to TSV.
	No return here.
	Takes a Path filename of the input file
		and the intended output file path as input.'''
	
	with open(tabfilepath, 'w') as outfile:
		writer = csv.writer(outfile, delimiter="\t")
		xlsfile = xlrd.open_workbook(filepath)
		sheet = xlsfile.sheet_by_index(0)
		for rownum in range(sheet.nrows):
			writer.writerow(sheet.row_values(rownum))
   
def get_input(context, desired_output):
 '''Function for handling user input.'''
 text_ok = False
 while not text_ok:
  text = input(context)
  if text.lower() not in ["y","n"]:
   print("Y for Yes or N for No, please.")
  else:
   text_ok = True
  
 if desired_output == "truefalse" and text.lower() == "y":
  return True
 if desired_output == "truefalse" and text.lower() == "n":
  return False
  
def write_log(query, query_dir_path):
 '''Saves metadata for a query to a log file.'''
 
 log_fn = "query.log"
 log_file_path = query_dir_path / log_fn
 
 now = datetime.datetime.now()
 nowstring = now.strftime("%Y-%m-%d %H:%M:%S")
 
 with open(log_file_path, "w", encoding="utf-8") as outfile:
  outfile.write(query + "\n")
  outfile.write(nowstring)
 
 return log_file_path
  
	
