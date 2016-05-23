#!usr/lib/python
#coding:utf8
#This is for Python 2.7
#This module is for general analysis of biomolecule like Protein, Gene, and so on.

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib
import requests
import sys
import os
import os.path

#instructions
#1. html5lib is required for pandas to read table from html.

#reference
#Amino acide table is from Wikipedia, Addgene and Invitrogen.

#data_base

AA_table = pd.read_pickle('AA_table.pickle')#collection of Amino acids info

complementory_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

seq_re = re.compile(r'[^ARNDCEQGHILKMFPSTWYV]+')#used for clean_seq function

gene_link_re = re.compile(r'<a href="/gene/(\d+)')#used in gene_search function

protein_link_re = re.compile(r'<a href="/uniprot/\w+?">(\w+?)</a>.*?>(\w+?)</td>')#used in uniprot_search function


#general functions

def clean_seq(s):
    '''
    Remove number, whitespace chars in a DNA or protein sequence.
    Return the cleaned sequence.
    '''
    return ''.join(seq_re.split(s.upper()))

def reverse_complementory(s):
    '''
    Get the reversed complementory sequence of a DNA sequence.
    Return the reversed complementory sequence.
    '''
    try:
        return ''.join([complementory_dict[base] for base in list(s.upper())])[::-1]
    except KeyError:
        print('Only "ATGCatgc" are accepted. Please check your sequence.')
        
def GC_content(s):
    '''
    Calculate GC content of DNA.
    '''
    s = s.upper()
    G_content = s.count('G')
    C_content = s.count('C')
    ratio = (G_content + C_content) / float(len(s))
    return format(ratio, '.2%')

def primers(s, length = 21, f_addon = '', r_addon = ''):
    '''
    This function is to design primers for a DNA sequence.
    '''
    s = s.upper()
    for base in s:
        if not base in ['A', 'T', 'G', 'C']:
            break
        return 'Illegal character found. Please check your sequence.'
    result = [[],[],[]]
    f_primer = f_addon + s[0:length]
    r_primer = r_addon + reverse_complementory(s[-length:])
    attr = [(str(length) + '-' + 'mer'), {'f_primer_GC':GC_content(f_primer)}, {'r_primer_GC':GC_content(r_primer)}]
    result[0].extend(attr)
    result[1].append({'f_primer':f_primer})
    result[2].append({'r_primer':r_primer})
    return result

def read_fasta(path):#read a local fasta file downloaded from NCBI
    with open(path) as foo:
        return ''.join([item.strip() for item in foo.readlines()[1:]])
    

def list_files(dir_name, type=None, *args):
    '''
    This function is used to list all the files in the requested dir. One can also get specific file types
    by defining the type parameter.
    The type can a string or a collection of strings that declearing file types.
    e.g. type = '.py', type = ['.txt', '.csv', '.data']
    '''
    if type is None:
        return [path for path in os.listdir(dir_name) if os.path.isfile(path)]
    return [path for path in os.listdir(dir_name) if os.path.isfile(path) and os.path.splitext(path)[1] == type]
    
def join_sequecning_fragments(file_list):
    '''
    When you send your PCR products or plasmid for sequencing, you probably get several .seq files back. 
    Usually, one need to align these fragments with your reference sequence individually.
    join_sequencing_fragments function provide you an easy way to get things done.
    Prerequisits:
    1. change working directory to the folder that contains the sequencing results
    2. make a file list using the list_files function with type parameter set as '.seq'
    For proper use of this function, please make sure your sequencing results are aranged in order. 
    In other words, up stream fragments should be prior to down stream fragments.
    '''
    seqs = []
    for each_file in file_list:
        with open(each_file) as foo:
            seq = foo.read().replace('\n', '')[100:]
            first_N = seq.find('N')
            if first_N <= 800:
                seq = seq[:first_N]
            else:
                seq = seq[:800]
            seqs.append(seq)
                
    joined_seq = seqs[0]
    
    for i in range(len(seqs)-1):
        idx = seqs[i+1].find(seqs[i][-20:])
        if idx == -1:
            return 'Something wrong between %d and %d seq. Make sure your .seq files are in order.'%(i, i+1)
        addon = seqs[i+1][idx+20:]
        joined_seq += addon
        
    return joined_seq
    

def help_info():
    print(
        '''
        *clean_seq(input_seq)*: remove white space or numbers in the input gene or protein sequence.
        *reverse_complementory(input_seq)*: get the reversed complementory sequence of the input DNA sequence.
        *GC_content(input_seq)*: calculate the GC content of the input sequence.
        
        --help or -h: get help info of Mybiotools
    ''')
        
def gene_search(gene_name, Species = 'Human'):
    '''
    This function is used to search genes by gene name and species on NCBI website.
    gene name is a required parameter, while species is set default to human.
    One can declare species when calling the function.
    A list of gene id and species will returned if serach sucessfully.
    '''
    url = 'http://www.ncbi.nlm.nih.gov/gene/'
    query = {'term':gene_name + ' ' + Species}
    r = requests.get(url, query)
    if r.status_code != 200:
        return 'Search failed.'
    gene_IDs = gene_link_re.findall(r.text)#return a list of the gene id
    
    res = re.compile(r'<.*?>')
    
    reg1 = re.compile(r'ID:.*?</span></td><td>(.*?)\[<em>')#get the gene names
    gene_names = reg1.findall(r.text)
    gene_names = [''.join(res.split(item)) for item in gene_names]
    
    reg2 = re.compile(r'\[<em>(.*?)\]')
    gene_species = reg2.findall(r.text)#find pattern like '[<em>Mus musculus</em> (house mouse)]'
    gene_species = [''.join(res.split(item)) for item in gene_species]
    
    return list(zip(gene_IDs, gene_names, gene_species))
        
def uniprot_search(protein_name, Species = 'Human', Uniprot_ID = None, *args, **kwargs):
    '''
    Search Uniprot for protein of interest.
    Parameters: protein_name(essential), Species set Human as default.
    Return: a list of protein uniprot ID and species.
    '''
    url = 'http://www.uniprot.org/uniprot/'
    data = {'query':protein_name + ' ' + Species, 'sort':'score'}
    r = requests.get(url, data)
    if r.status_code != 200:
        return 'Search failed.'
    
    alist = protein_link_re.findall(r.text)
    protein_id = [item[0] for item in alist]
    entry_name = [item[1] for item in alist]
    
    species_re = re.compile(r'<a href="/taxonomy/.*?">(.+?)</a>')
    species_ = species_re.findall(r.text)
    
    return list(zip(protein_id, entry_name, species_))

def psipred(protein_seq, subject, email=None, passwd=None):
    url = 'http://bioinf.cs.ucl.ac.uk/psipred/submit'
    payload = {
        'utf8':'%E2%9C%93', 'program_psipred':'1', 'sequence':protein_seq, 'email':email, 'passwd':passwd, 
        'subject':subject, 'commit':'Predic', 'msa_control':'all', 'output':'opnone', 'seqalign':'yes', 
        'database':'PfamA', 'eval':'0.01', 'iterations':'5', 'domssea':'yes', 'secpro':'yes', 'pp':'yes'
    }
    r = requests.get(url, payload)
    if r.status_code != 200:
        return 'Can not fetch the site.'
    reg = re.compile(r'check the progress at <a href="(.*?)">')
    lnk = reg.findall(r.text)[0]
    if lnk is None:
        return 'No match returned. Check with Mybiotools Author for this issue.'
    print('Job has been submitted to the server. You can check it later at the link below. To be noted, it typically takes at least 30 minutes and may be as long as 2 hours.')
    return lnk
        
#run from shell
    
accepted_argvs = {'clean_seq':clean_seq, 'reverse_complementory':reverse_complementory, 'GC_content':GC_content, '--help':help_info, '-h':help_info}


#classes

class isoforms():
    def __init__(self, code):
        self.name = None
        self.description = None
    
class Gene:
    def __init__(self, gene_name, Species = 'Human', Gene_ID = None, *args, **kwargs):
        self.name = gene_name
        self.Species = Species
        init_list = gene_search(self.name, Species = self.Species)
        if init_list == 'Search failed.':
            return 'Connection to internet failed.'
        elif init_list is None:
            return 'Get nothing back from NCBI gene search. Please check your input.'
        print('Please look up in the table below and set your gene of interest.')
        print(pd.DataFrame(init_list, columns = ['Gene ID', 'Gene description', 'Species']))
        user_choice = input('Please select the index number(starts from 0) of the gene of your interest:  ')
        while True:
            if int(user_choice) < -1 or int(user_choice) >= len(init_list):
                print('Invalid number, please input again.')
            elif not user_choice.isnumeric():
                print('Only numbers are accepted.')
            elif int(user_choice) == -1:
                return 'You decide there is no gene of your choice. Please try other key words.'
            else:
                break
                
        self.description = init_list[int(user_choice)][1]
        self.ID = init_list[int(user_choice)][0]
        self.Species = init_list[int(user_choice)][2]
        print('\n\nYour gene object %s from %s has been created successfully.'%(self.description, self.Species))
        
        self.mRNA = None
        self.default_mRNA = None
        self.NG_num = None
        
    def get_Gene_ID(self):
        print(self.ID)
        
    def get_NG_num(self):
        url_base = 'http://www.ncbi.nlm.nih.gov/gene/'
        r = requests.get(url_base + str(self.ID))
        if r.status_code != 200:
            return 'Connection to NCBI failed.'
        NG_re = re.compile(r'<p>NG_(.*?)RefSeqGene</p>')
        self.NG_num = NG_re.findall(r.text)[0].strip()
        return self.NG_num
    
    def get_mRNA(self):
        url_base = 'http://www.ncbi.nlm.nih.gov/gene/'
        r = requests.get(url_base + str(self.ID))
        if r.status_code != 200:
            return 'Connection to NCBI failed.'
        mRNA_re = re.compile(r'<p><a href="/nuccore/(.*?)">NM_')
        mRNA_list = mRNA_re.findall(r.text)
        
        if mRNA_list:
            self.mRNA = list(range(len(mRNA_list) + 1))
            self.mRNA[0] = 'total: %d'%len(mRNA_list)
            for i in range(1, len(mRNA_list)+1):
                self.mRNA[i] = isoforms(self.mRNA[i])
        
        mRNA_name_re = re.compile(r'NP_.*?</a>(.*?)</p>')
        mRNA_names = mRNA_name_re.findall(r.text)
        mRNA_names = [item.strip() for item in mRNA_names]
        
        mRNA_des_re = re.compile(r'<dd>Transcript Variant:(.*?)</dd>')
        mRNA_descriptions = mRNA_des_re.findall(r.text)
        mRNA_descriptions = [item.strip() for item in mRNA_descriptions]
        
        for i in range(1, len(mRNA_list)+1):
            self.mRNA[i].name = mRNA_names[i-1]
            self.mRNA[i].description = mRNA_descriptions[i-1]
            self.mRNA[i].NM_id = mRNA_list[i-1]
            
        mRNA_df = pd.DataFrame([(item.NM_id, item.name, item.description) for item in self.mRNA[1:]], columns = ['NM_id', 'name', 'description'])
        print(mRNA_df)
        
        self.default_mRNA = self.mRNA[1]
        
    def get_mRNA_seq(self):
        if self.default_mRNA is None:
            return 'Please use get_mRNA method to set the default mRNA first.'
        url_base = 'http://www.ncbi.nlm.nih.gov/nuccore/'
        r = requests.get(url_base + str(self.default_mRNA.NM_id))
        if r.status_code != 200:
            return 'Connection to NCBI failed.'
        
        nuccore_re = re.compile(r'genbank_fasta" href="/nuccore/(.*?)\?report=fasta"')
        nuccore_id = nuccore_re.findall(r.text)[0]
        
        _url = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi'
        payload = {'val':nuccore_id, 'db':'nuccore', 'dopt':'fasta', 'extrafeat':'0', 'fmt_mask':'0', 'retmode':'html', 'withmarkup':'on', 'log$':'seqview', 'maxdownloadsize':'1000000'}
        
        fasta_seq = requests.get(_url, payload).text
        
        reg1 = re.compile(r'<s.*?">')
        reg2 = re.compile(r'if.*</script>')
        
        fasta1 = ''.join(reg1.split(fasta_seq))
        fasta2 = ''.join(reg2.split(fasta1))
        fasta3 = ''.join(re.split(r'</span>', fasta2))#mRNA in fasta format
        mRNA_seq = ''.join(re.split(r'\n+', fasta3)[1:])
        self.default_mRNA_seq = mRNA_seq
        
        payload1 = {'val':nuccore_id, 'db':'nuccore', 'dopt':'genbank', 'extrafeat':'976', 'fmt_mask':'0', 'retmode':'html', 'withmarkup':'on', 'log$':'seqview', 'maxplex':'3', 'maxdownloadsize':'1000000'}
        
        cds_loc = requests.get(_url, payload1).text
        cds_re = re.compile(r'features\["CDS"\].push\(\[\[(.*?)\]\]\);')
        cds_index = cds_re.findall(cds_loc)[0].strip().split(',')
        cds_seq = mRNA_seq[int(cds_index[0])-1:int(cds_index[1])]
        self.cds_seq = cds_seq
        
        return 'The default mRNA seq and CDS seq have been set as %s\'s attributes.'%self.name
            
        
    
class Protein:
    def __init__(self, protein_name, Species = 'Human', Uniprot_ID = None, seq = None, *args, **kwargs):
        self.name = protein_name
        self.Species = Species
        self.seq = seq
        if self.seq:
            return 'The input sequence will be used. No Uniprot Search will be performed.'
            
        candidates = uniprot_search(self.name, Species = self.Species)
        candi_df = pd.DataFrame(candidates, columns = ['UniProt_ID', 'Entry Name', 'Species'])
        print(candi_df)
        user_choice = input('Please look into the table above and make your choice by input the row index of the protein of your interest. To be noted, index starts from 0.\n\n Input a number between 0 and %d here:  '%(len(candidates)-1))
        while True:
            if not user_choice.isnumeric():
                print('Input a number between 0 and %d please. Or input -1 if there is no protein of your choice.'%(len(candidates)-1))
            elif int(user_choice) < -1 or int(user_choice) >= len(user_choice):
                print('Input a number between 0 and %d please. Or input -1 if there is no protein of your choice.'%(len(candidates)-1))
            elif int(user_choice) == -1:
                return 'You decide there is no your choice. Please try other key words.'
            else:
                break
        
        self.UniProt_id = candidates[int(user_choice)][0]
        self.entry_name = candidates[int(user_choice)][1]
        self.Species = candidates[int(user_choice)][2]
        
        print('Your protein object %s has been created successfully.'%(self.entry_name))

        
    def get_seq(self):
        if self.seq:
            print('The input sequence is being used.')
            return self.seq
        url = 'http://www.uniprot.org/uniprot/'
        response = requests.get(url + self.UniProt_id + '.fasta')
        self.fasta = response.text
        self.seq = ''.join(self.fasta.split('\n')[1:])
        
        return self.seq
    
    def get_isoforms(self):
        url = 'http://www.uniprot.org/uniprot/'
        response = requests.get(url + self.UniProt_id)
        iso_re = re.compile(r'Isoform \d')
        return set(iso_re.findall(response.text))
    
    def get_function(self):
        url = 'http://www.uniprot.org/uniprot/'
        response = requests.get(url + self.UniProt_id)
        Func_re = re.compile(r'<meta content="(.*?)" name="description"')
        return Func_re.findall(response.text)[0]

    def get_modifications(self):
        pass
    
    def pred_2d_structure(self):
        pass
    
    def get_3d_structure(self):
        pass