#!/python
import sys
import operator
import Bio.Entrez as BE
import Bio.SeqIO as SeqIO


'''
reads a .el file and return a set of all nodes

@param: filename of .el file

@return: a set (complete and non-repeating) of all nodes (gene ids)
'''

'''
def readPPI(fileName):
    ofile = open(fileName, "r")

    vertices = set()
    
    for line in ofile:
        nodes = line.split('\t')
        for node in nodes:
            node_n = node.split('_')[1]
            vertices.add(node_n)
    
    ofile.close()
    return vertices
'''


def readPPI(fileName):
    ofile = open(fileName, "r")

    vertices = set()
    id_dict = dict()
    return_me = []
    
    for line in ofile:
        nodes = line.strip().split('\t')
        for node in nodes:
            node_n = node.split('_')[1]
            vertices.add(node_n)
            id_dict[node_n] = node.split('_')[0]

    return_me.append(vertices)
    return_me.append(id_dict)
    
    ofile.close()
    return return_me



'''
Take all the gene ids read from a .el file and output proper protein seq corresponds
to each gene. Longer, ref seq has higher priority to be chosen

@param: filename of the output fasta file
@param: set of all entrez_ids (aka gene ids)

@output: prints
'''
def fasta_from_entrezid(fileName, entrez_ids):
    ofile = open(fileName, "w")
    count = 1
    for eid in entrez_ids[0]:
        try:
            
            print("processing ", count, " / ", str(len(entrez_ids[0])) + ':', eid)

            # find all the protein ids of the current gene
            res = BE.elink(db='protein', dbfrom='gene', id=eid)
            record = BE.read(res)
            res.close()
            
            # Get a list of dict that contains 'Link', 'DbTo', 'LinkName' as keys
            linkset = record[0]['LinkSetDb']
            
            # locates the dict with LinkName gene_protein, initialized to -1
            normal_pos = -1
            # locates the dict with LinkName gene_protein_ref_seq, initialized to -1
            ref_pos = -1
            
            # find the dict of gene proteins seqs and refseq(s)
            for i in range(0, len(linkset)):
                if linkset[i]['LinkName'] == "gene_protein":
                    normal_pos = i
                elif linkset[i]['LinkName'] == "gene_protein_refseq":
                    ref_pos = i
            
            # holder for the seq we choose to output
            chosen_seq = ''
            
            # if the ref seq(s) exists
            if ref_pos != -1:
                # find and keep the longest ref seq
                for id_dict in linkset[ref_pos]['Link']:
                    pid = id_dict['Id']
                    handle = BE.efetch(db='protein', id=pid, rettype='fasta')
                    seq = str(SeqIO.read(handle, 'fasta').seq)
                    
                    if len(seq) > len(chosen_seq):
                        chosen_seq = seq
            # otherwise, find and keep the longest seq
            else:
                for id_dict in linkset[normal_pos]['Link']:
                    pid = id_dict['Id']
                    handle = BE.efetch(db='protein', id=pid, rettype='fasta')
                    seq = str(SeqIO.read(handle, 'fasta').seq)
                    
                    if len(seq) > len(chosen_seq):
                        chosen_seq = seq

            # output the chosen protein seq
            #ofile.write(">%s\n"%(eid))
            ofile.write('>' + eid + '_' + entrez_ids[1][eid] + '\n')
            ofile.write("%s\n"%(chosen_seq))
        
        except Exception as e:
            print(e)
            print("seq: ", eid, "not found")
            
        count = count + 1
        
    ofile.close()

BE.email = 'mluo@iwu.edu'

PPI_File = sys.argv[1] #  <-- input .el PPI network
Fasta_File = sys.argv[2] #  <-- output .fasta file


print("Loading entrez-ids from PPI file", PPI_File)

# get a set of gene ids read from a PPI file
entrez_list = readPPI(PPI_File)


print(len(entrez_list[0]))


print("Retrieving FASTA sequences", "\n")

# output all protein seqs correspond to the gene ids
nodes = fasta_from_entrezid(Fasta_File, entrez_list)





# TODO:
# ADD COMMENTS, RENAME VARIABLES MEANINGFULLY
# Run on all species and check that all get a fetched protein sequence
# Spot-check protein FASTA files, make sure they downloaded properly
# Maybe TODO:
# Always pull the protein ID from the gene_protein_refseq link if it exists
# If there is more than 1 protein link, use the LONGEST protein sequence
# Can this be sped up? I.e. fetch all the elinks at once rather than submitting as separate requests
# Similarly for fetching multiple FASTAs, if we need to, and/or determining which protein is longest without fetching all of them
# Next:
# Change script to submit jobs to run in parallel
# Add all network species from networks
# Run scripts
# 



# limited version
'''
def fasta_from_entrezid(fileName, entrez_ids):
    ofile = open(fileName, "w")
    count = 1
    for eid in entrez_ids:
        try:
            #print 'processing ',count, ' / ', len(entrez_ids)
            print("processing ", count, " / ", str(len(entrez_ids)) + ':', eid)

            res = BE.elink(db='protein', dbfrom='gene', id=eid)
            record = BE.read(res)
            res.close()
            
            pid = record[0]['LinkSetDb'][0]['Link'][0]['Id']
            #pid = mylst[1]
            handle = BE.efetch(db='protein', id=pid, rettype='fasta')
            seq = str(SeqIO.read(handle, 'fasta').seq)
            ofile.write(">%s\n"%(eid))
            ofile.write("%s\n"%(seq))
        except Exception as e:
            print(e)
            #print "seq: %s not found"%(eid)
            print("seq: ", eid, "not found")
            
        count = count + 1
    ofile.close()
'''
