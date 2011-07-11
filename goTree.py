import gene_ontology as go
from gzip import GzipFile
from collections import defaultdict

def annotations(fileObj, evidenceCodes = []):
	'''
	Takes only 3 things from the annotation file, db_object_id and go_id, while evidence_code are used for filterDataing
	Returns a dict mapping go_id to a set of db_object_id
	
	fileObj - file containing annotations, can be filename, or file object
	evidenceCodes - the evidence codes to select, if the list is left empty, then all are allowed
	'''
	
	ann = go.AnnotationFile(fileObj)
	
	
	if len(evidenceCodes) == 0:
		filterData = False
	else:
		filterData = True
		evidenceCodes = set(evidenceCodes)
		
	associations = defaultdict(set)	# Maps go_id to a set of genes id
	for temp in ann:
		#filterDataing based on evidence codes
		if filterData  and temp.evidence_code not in evidenceCodes:
			continue
		
		associations[temp.go_id].add(temp.db_object_symbol)		#IMPORTANT: change to db_object_id
	
	return associations
		
		
	
	
	
