echo "Testing ProtExAM."
(cd protexam/ && python protexam.py --auto --query "(Barth Syndrome[MeSH Terms]) AND (("2000/01/01"[Date - Publication] : "3000"[Date - Publication]))") && echo "Ran without errors." || echo "Failed!"
