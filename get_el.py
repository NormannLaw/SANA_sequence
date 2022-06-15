from typing import List
import sys
from Bio import Entrez

# *Always* tell NCBI who you are
Entrez.email = "mluo@iwu.edu"


def get_gene_ids(line: str) -> List:
    # list for later return
    return_me = []
    # create a list contains strs after line splitted by tabs
    by_tab = line.split('\t')

    # get gene id of interactor A and B from col 0 and 1
    for i in range(0, 2):
        by_colon = by_tab[i].split(':')
        gene_id = by_colon[1]
        return_me.append(gene_id)

    return return_me



def retrieve_annotation(gene_id):

    
    handle = Entrez.esummary(db = "gene", id = gene_id) 
    annotations = Entrez.read(handle)

    return annotations['DocumentSummarySet']['DocumentSummary'][0]['Name']



read_me = sys.argv[1] #  <-- input file
el_file = sys.argv[2] #  <-- output .el
el_n_file = sys.argv[3] # <-- output .el file with only number ids


id_set = set()
id_dict = dict()


with open(read_me) as r:
    with open(el_file, 'w') as e:
        with open(el_n_file, 'w') as k:
            r.readline()
            for line in r:
                ids = get_gene_ids(line)
                
                if ids[0] not in id_set:
                    name0 = retrieve_annotation(ids[0])
                    id_set.add(ids[0])
                    id_dict[ids[0]] = name0
                else:
                    name0 = id_dict[ids[0]]
                
                if ids[1] not in id_set:
                    name1 = retrieve_annotation(ids[1])
                    id_set.add(ids[1])
                    id_dict[ids[1]] = name1
                else:
                    name1 = id_dict[ids[1]]

                k.write(ids[0] + '\t' + ids[1] + '\n')
                #e.write(ids[0] + '_' + name0 + '\t' + ids[1] + '_' + name1 + '\n')
                e.write(name0 + '_' + ids[0] + '\t' + name1 + '_' + ids[1] + '\n')
