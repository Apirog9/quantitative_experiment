# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 15:30:27 2023

@author: APirog
"""

import os
import json
import pandas as pd
from .processing_utils import Processing
from .annotation_utils import Annotation
import itertools as it
import copy
import functools as func


column_default_dictionary = {'precursor_level' : {"precursor" : "Precursor.Id","unmodified_peptide":"Stripped.Sequence","modified_peptide":"Modified.Sequence","protein_group":"Protein.Group","protein_name":"First.Protein.Description","gene_name":"Genes"},
'peptide_level': {"unmodified_peptide":"Stripped.Sequence","modified_peptide":"Modified.Sequence","protein_group":"Protein.Group","protein_name":"First.Protein.Description","gene_name":"Genes"},
'protein_level': {"protein_group":"Protein.Group","protein_name":"First.Protein.Description","gene_name":"Genes"},
'GO_term' : {}}
multiple_item_delimiter = ';'


class Experiment:
    '''
    Contains a quantitative experiment class. The quantitation basis is the precursortable.
    datapath - path to store the experiment data
    mzml_dict - json derived object, dict with keys:
        "files" - mzml - quatitation column pairing
        "path"  - path to files
    conditions - currently supports conditions, biological replicates and technical replicates.
                json derived object, dict of dicts, keys:
        'raw' - conditions vs raw (not averaged between technical replicates)
        'technical_replicates' - keys are conditions, values are lists to average to obtain conditions
        'averaged' - conditions after averaging, empty dict before actual averaging
    metadata - any metadata in form of dict (also json-derived)
    precursortable - raw, unaveraged precursortable
    processed_tables - tables produced after averaging, annotating, and other transformations. dict of filenames or dataframe objects
    associatedlevels - peptide, protein and higher concatenation levels. All levels are defined using proper objects
    default_data - which data table are currently set for default use
    All tables are stored and operated as:
        {'precursortable' : [table,self.datapath+'/raw_precursortable.tsv']}
        {"name" : [DataFrame object,filepath]}
        saving writes dataframe object to file, and reading writes data from file to DataFrame
        
    Helps creation of annotated quantitation results in reproducible way.
        
    '''
    
    def __init__(self):
        self.datapath = ''
        self.mzml_dict = {}
        self.conditions = {}
        self.metadata = {}
        self.tables = {}
        self.associatedlevels = {}
        self.default_data = None
        self.rawquantcolumns = []
        self.processedquantcolumns = []
        processing = Processing()
        
        
    def set_datapath(self,pathlike):
        '''
        get datapath string from pathlike, if folder does not exist, create it
        pathlike -  path to store data , format d:/somecatalog/some
        
        sets self.datapath attribute
        '''
        datapath = pathlike
        if not os.path.exists(datapath):
            os.mkdir(datapath)
        self.datapath = datapath

        
        
    def set_mzmls(self,jsonlike):
        ''''
        set and check the existence of mzml files
        jsonlike - path to json containing mzml data
            
        set self.mzmls attribute
        '''
        jsonlike = open(jsonlike,'r')
        data = json.load(jsonlike)
        'totalpah is the directory of files'
        totalpath = data['mzml']['path']
        'files are individual file names'
        files = data['mzml']['files']
        filedict = {}
        'make full path file names'
        for key in files.keys():
            filedict[key] = totalpath +'/' + files[key]
        'check if they exist'    
        for key in filedict.keys():
            if os.path.exists(filedict[key]):
                pass
            else:
                print(key)
                print(filedict[key])
                print('Does not exist')
        self.mzml_dict = filedict
        

    def set_conditions(self,jsonlike):
        ''''
        read and set conditions
        jsonlike - path to json containing condition data
            
        set self.conditions attribute
        '''
        jsonlike = open(jsonlike,'r')
        data = json.load(jsonlike)
        self.conditions = data['rawconditions']
        self.rawquantcolumns = list(it.chain.from_iterable([self.conditions['raw'][key] for key in self.conditions['raw'].keys()]))
        

        
    def set_metadata(self,jsonlike):
        ''''
        read and set metadata
        jsonlike - path to json containing metadata
            
        set self.metadata attribute
        '''
        
        jsonlike = open(jsonlike,'r')
        data = json.load(jsonlike)
        
        self.metadata = data['metadata']
        
        

    def set_table(self,pathlike,datatype='precursor_level',name='precursortable',additional_required = []):
        '''
        read and check table. rename critical columns
        pathlike - path to table
        required_columns - columns required to be present in the table, useful to prevent importing table with wrong column names
        name -  key for self.tables dictionary
            
        adds a key 'name' : value [DataFrame,tsv file] pair to self.tables dictionary
        '''
        default_columns = column_default_dictionary[datatype]
        if self.metadata['columns'] == 'default':
            defaultrenamer = {}
        else:
            defaultrenamer = {}
            for key in default_columns.keys():
                defaultrenamer[self.metadata['columns'][key]] = default_columns[key]
           
        table = pd.read_csv(pathlike,sep='\t')
        columns = list(table.columns)
        required_columns = self.rawquantcolumns + list(defaultrenamer.keys())

        if self.metadata == {}:
            print('Define metadata first')
        if self.conditions == {}:
            print('Define conditions first')
        else:
            if any([x not in columns for x in required_columns]):
                print('Check column names')
            else:
                filename = self.datapath+'/raw_' +name + '.tsv'
                table = table.rename(defaultrenamer,axis=1)
                self.tables[name] = [table,filename]
                self.tables[name][0].to_csv(filename, sep='\t',index=False)
                self.default_data = self.tables[name][0]
        
        
    def __str__(self):
        'print current object'
        return 0
    
    
    def savestatus(self,pathlike=None):
        '''
        save the whole content
        '''
        if not pathlike:
            pathlike = self.datapath
            
            
        
    def loadstatus(self,pathlike):
        '''
        load the whole content
        
        '''
        
        
    def average_technical(self,title):
        '''
        average technical replicates, store result
        title -  table key from self.tables
        use self.conditions['technical_replicates'] dictionary to define technical replicates
        stores new table of name 'technical_averaged_'+title and file name technical_averaged_'+title+'.tsv'
        '''
        frame = copy.copy(self.tables[title][0])
        technical_averaged = Processing.average_technical(frame, self.conditions['technical_replicates'])
        technical_averaged.to_csv(self.datapath+'/technical_averaged_'+title+'.tsv', sep='\t',index=False)
        self.tables['technical_averaged_'+title] = [technical_averaged,self.datapath+'/technical_averaged_'+title+'.tsv']
        self.quantcolumns = list(self.conditions['technical_replicates'].keys())
        return technical_averaged
        
    def associate_precursor_peptide(self,peptideobject,title = 'peptide_level'):
        '''associate matching peptidelevel data
        peptideobject - object of type Peptidelevel
        id_precursor  - column containing unique! precursor identifiers
        id_peptide  - column containing unique! precursor identifiers
        title - quantification level title , default 'peptide_level'
        
        sets the quantification level data in self.associatedlevels dictionary with key 'title', data contains a dictionary with keys:
        'object' - Peptidelevel object
        'a_key'  - precursor identifier
        'b_key'  - peptide identifier
        'conection_rule' - ['many','one'] - defines a rule necessary for concatenation of precursor and peptide tables
        
        '''
        pepleveldata = {}
        pepleveldata['object'] = peptideobject
        pepleveldata['a_key'] = "Precursor.Id"
        pepleveldata['b_key'] = "Modified.Sequence"
        pepleveldata['connection_rule'] = ['many','one']
        self.associatedlevels[title] = pepleveldata


    def associate_precursor_protein(self,proteinobject,title = 'protein_level'):
        '''associate matching peptidelevel data
        proteinobject - object of type Peptidelevel

        title - quantification level title , default 'peptide_level'
        
        sets the quantification level data in self.associatedlevels dictionary with key 'title', data contains a dictionary with keys:
        'object' - Peptidelevel object
        'a_key'  - precursor identifier
        'b_key'  - peptide identifier
        'conection_rule' - ['many','one'] - defines a rule necessary for concatenation of precursor and peptide tables
        
        '''
        protleveldata = {}
        protleveldata['object'] = proteinobject
        protleveldata['a_key'] = "Precursor.Id"
        protleveldata['b_key'] = "Protein.Group"
        protleveldata['connection_rule'] = ['many','one']
        self.associatedlevels[title] = protleveldata        
        
        
        
    
    def save_tables(self,which = 'all'):
        '''
        saves  tables to tsv
        save existing DataFrame object to file
        which - list of titles (self.tables keys) to save, or 'all' - save all tables
        
        saving and empty dataframe is prevented
        '''
        if which == 'all':
            tables_to_save = list(self.tables.keys())
        else:
            tables_to_save = which
            
        for tabletitle in tables_to_save:
            if self.tables[tabletitle][0].shape != (0,0):
                self.tables[tabletitle][0].to_csv(self.tables[tabletitle][1],sep='\t',index=False)
            else:
                print('table empty, not saved')
        
    
    def laod_tables(self,which = 'all'):
        '''
        loads  tables to dataframes
        loads .tsv files to respective DataFrame objects
        which - list of titles (self.tables keys) to load, or 'all' - load all tables
        '''
        if which == 'all':
            tables_to_load = list(self.tables.keys())
        else:
            tables_to_load = which
        for tabletitle in tables_to_load:
            self.tables[tabletitle][0] = pd.read_csv(self.tables[tabletitle][1],sep='\t')
            
     
    
    def purge_dataframe_objects(self,which='all'):
        '''saves all tables to files, removes from memory
            will save all tables, then set all dataframe objects to empty dataframe.
            
            which - list of titles (self.tables keys) to purge, or 'all' - purge all tables
            
            
            may be handy to free memory if some frame is exceptionally huge
        
        '''
        if which == 'all':
            tables_to_purge = list(self.tables.keys())
        else:
            tables_to_purge = which
        self.save_tables(tables_to_purge)
        for tabletitle in tables_to_purge:
            self.tables[tabletitle][0] = pd.DataFrame()
        
        
        return 0
    
    def compound_levels(self,levels,level_tables,keypairs,quantcolumns = 'processed',multiple_delim=multiple_item_delimiter):
        
        '''
        match quantification levels for single experiment
        level_tables - level tables 
        keypairs - connection keys for each pair, starting from the left
        
        '''
        #TODO genreate proper column names for keys
        newkeypairs = []
        for num,item in enumerate(keypairs):
            new_a = column_default_dictionary[levels[num]][item[0]]
            new_b = column_default_dictionary[levels[num]][item[1]] 
            newkeypairs.append((new_a,new_b))
        keypairs = newkeypairs
        print(newkeypairs)
        traverse = {}
        
        for i in range(len(keypairs)):
            if levels[i] != 'precursor_level':
                table_a = self.associatedlevels[levels[i]]['object'].tables[level_tables[i]][0]
                print(self.associatedlevels[levels[i]]['object'].tables[level_tables[i]][1])
                
            else:
                table_a = self.tables[level_tables[i]][0]
            if levels[i+1] != 'precursor_level':
                table_b = self.associatedlevels[levels[i+1]]['object'].tables[level_tables[i+1]][0]
                print(self.associatedlevels[levels[i+1]]['object'].tables[level_tables[i+1]][1])
            else:
                table_b = self.tables[level_tables[i+1]][0]
                
            
            key_a = keypairs[i][0]
            key_b = keypairs[i][1]
            links_a_b,links_b_a = Processing.create_links(table_a, table_b, key_a, key_b,multiple_delim)
            traverse[(levels[i],levels[i+1])] = links_a_b
            traverse[(levels[i+1],levels[i])] = links_b_a
            
            
        self.traversion_map = traverse
        
        return self.traversion_map
    
    def synthetize_proteins(self):
        '''
        create artificial proteintable without quant for precursor/peptide linking
        useful if normally quantified proteintable does not make sense
        '''
        
        
        return 0
        
        
        
        
        
class Peptidelevel(Experiment):
    
    def __init__(self,experiment):
        ''''using experiment class created beforehand ensures sharing raw data, metadata, conditions and datapath
        experiment - object of Experiment class, to work properly, should have defined
        mzml_dict,metadata,conditions,datapath,quantcolumns (quantitation columns after averaging of technical replicates)
        and rawquantcolumns(raw quantification columns before averaging)
        '''
        
        self.mzml_dict = experiment.mzml_dict
        self.metadata = copy.copy(experiment.metadata)
        self.conditions = experiment.conditions
        self.datapath = experiment.datapath
        self.quancolumns = experiment.quantcolumns
        self.rawquantcolumns = experiment.rawquantcolumns
        self.randomized_peptidome = None
        self.n_term_cut_motifs = None
        self.c_term_cut_motifs = None
        self.tables = {}
        self.default_data = None
        
    def update_metadata(self,jsonlike):
        jsonlike = open(jsonlike,'r')
        data = json.load(jsonlike)
        metadata = data['metadata']
        self.metadata['name'] = metadata['name']
        self.metadata['columns']=metadata['columns']
        
        
        
    def create_random_peptidome(self,peptidetable,fasta,number,len_range,make_weights = 'from_self',proteometable=None):
        """
        For proper calculation, the number of peptides should be high enough to generate a smooth histogram. 
        (Histograms are plotted when function is called)
        
        
        peptidetable - right now accepts a DIANN quantification result. 
        
        peptidecolumn - column with unmodified peptide sequences
        peptidecolumn - column with protein ids, separated by ; if multiple ids per line
        
        fasta - FASTA file containing protein sequences. For fragpipe, bet to use protein.fas file stored in results
        
        number - number of random peptides to generate
        
        len_range - length range of specified peptides. A tuple.
        
        make_weights - 'from_self' will adjust the frequency of protein selection relying of numper of peptides quantified for protein
                       'from proteome' - will do the same using external tsv file containing at two columns - 'Razor Spectral Count' and 'Protein ID'
                                         (this fit FragPipe search results withot modification, for any other, just rename columns in file, and remove NaNs from 'Razor Spectral Count')
                       False -  do not make weigths, use all proteins from fasta
        create a table with two identical! columns 'Stripped.Sequence' and 'Modified.Sequence' for compatiblity fith other functions
        additionally saves two files:
            .list file with one peptide per line, and .args file with arguments used for creation of random peptidome
        """
        
        table = self.tables[peptidetable][0]
        if table.shape == (0,0):
            self.load_tables(which = [peptidetable])
            table = self.tables[peptidetable][0]
        random_peptidome  = Annotation.generate_random_peptidome(table, 'Stripped.Sequence', 'Protein.Group', fasta, number, len_range,make_weights=make_weights,proteometable = proteometable)
        arguments = (peptidetable, 'Stripped.Sequence', 'Protein.Group', fasta, number, len_range,make_weights,proteometable)
        random_peptidome_frame = pd.DataFrame()
        random_peptidome_frame['Stripped.Sequence'] = random_peptidome
        random_peptidome_frame['Modified.Sequence'] = random_peptidome
        self.randomized_peptidome = random_peptidome_frame
        
        with open(self.datapath+'/random_peptidome_'+peptidetable+'.list','w') as outfile:
            for peptide in random_peptidome:
                outfile.write(peptide +'\n')
                
        with open(self.datapath+'/random_peptidome_'+peptidetable+'.args','w') as outfile:
            outfile.write(str(arguments))
                
                
        random_peptidome_frame.to_csv(self.datapath+'/random_peptidome_'+peptidetable+'.tsv',sep='\t',index=False)
                
    def annotate_peptide_position(self,peptidetable,fasta,title = 'position_annotated_peptidetable'):
        
        '''
        Annotates peptide position in protein sequence(s)
        
        peptidetable - table with peptide sequences 
        peptidecolumn - column with unmodified peptide sequences
        fasta - fasta file
        
        function does not assume I and L as identical
        
        returns table with two additional coolumns - 'Starts','Ends' containing possibe start and end positions of peptide 
        protein sequences that contain this peptide
        
        '''
        
        table = copy.copy(self.tables[peptidetable][0])
        if table.shape == (0,0):
            self.load_tables(which = [peptidetable])
            table = self.tables[peptidetable][0]
            
        position_annotated_table = Annotation.annotate_protein_position(table, 'Stripped.Sequence', fasta)
        arguments = (peptidetable,'Stripped.Sequence',fasta,title)
        tablepath = self.datapath + '/peptide_position_annot_'+peptidetable +'.tsv'
        position_annotated_table.to_csv(tablepath,sep='\t',index=False)
        
        with open(self.datapath + '/peptide_position_annot_'+peptidetable +'.args','w') as outfile:
            outfile.write(str(arguments))
        
        
        self.tables[title] = [position_annotated_table,tablepath]
        
    def annotate_cut_motifs(self,peptidetable,motif_n_extension_in,motif_n_extension_out,motif_c_extension_in,motif_c_extension_out,fasta,i_l_identical = False,title='motif_annotated_peptidetable'):
        
        '''        
        extract possible sequence motifs around n and c termini of peptides
        
        frame -  table with peptides (modified and unmodified sequences)
        peptide_column - columns with unmodified peptide sequences
        modified_peptide_column - column with modified peptide sequences
        motif_n_extension_in,motif_n_extension_out,motif_c_extension_in,motif_c_extension_out -  controls the position
        of extracted sequences around N and C terms, with use of negative numbers a motif totally outsie or totally inside peptide may be defined
        fasta -  fasta file containing protein sequences from which peptides derive
        i_l_identical -  if treat I and L as identical for peptide:protein assignment
        '''
        
        
        
        
        
        table = copy.copy(self.tables[peptidetable][0])
        if table.shape == (0,0):
            self.load_tables(which = [peptidetable])
            table = self.tables[peptidetable][0]
            
        motif_annotated_table = Annotation.annotate_cut_motif(table, 'Stripped.Sequence', 'Modified.Sequence', motif_n_extension_in, motif_n_extension_out, motif_c_extension_in, motif_c_extension_out, fasta,i_l_identical = i_l_identical)
        arguments = (peptidetable,'Stripped.Sequence', 'Modified.Sequence',motif_n_extension_in,motif_n_extension_out,motif_c_extension_in,motif_c_extension_out,fasta,i_l_identical,title)
        tablepath = self.datapath + '/peptide_motif_annot_'+peptidetable +'.tsv'
        motif_annotated_table.to_csv(tablepath,sep='\t',index=False)
        
        with open(self.datapath + '/peptide_motif_annot_'+peptidetable +'.args','w') as outfile:
            outfile.write(str(arguments))
        
        self.tables[title] = [motif_annotated_table,tablepath]
        
        
    def extract_motifs(self, motif_frame,max_motifs_c,max_motifs_n,minlength,mods_to_remove):
        '''
        Extracts motif set from table with previously annotated motifs
        motif_frame -table with motifs annotated by get_motif() function
        modified_peptide_column - column with modified peptide sequences
        max_motifs_c, max_motifs_n - maximum allowed number of motifs per peptide. To prevent
                                    overrepresentation motifs from e.g. protein with multiple isoforms
        minlength - minimum length of motif. some motifs may be shortened due to protein N- or C- term proximity
        mods_to_remove - remove motifs from peptides containing specified modifications
        
        
        returns two lists of n-term and c-term motifs
        '''
        table = self.tables[motif_frame][0]
        if table.shape == (0,0):
            self.load_tables(which = [motif_frame])
            table = self.tables[motif_frame][0]
        n_terms,c_terms = Annotation.extract_motif_set(table, 'Modified.Sequence', max_motifs_c, max_motifs_n, minlength,mods_to_remove = mods_to_remove)
        arguments = (motif_frame,'Modified.Sequence',max_motifs_c,max_motifs_n,minlength,mods_to_remove)
        self.n_term_cut_motifs = n_terms
        self.c_term_cut_motifs = c_terms
        with open(self.datapath + '/' + motif_frame+"_nterms.txt","w") as outfile:
            for item in n_terms:
                outfile.write(item+ "\n")
        with open(self.datapath + '/' + motif_frame+"_cterms.txt","w") as outfile:
            for item in c_terms:
                outfile.write(item+ "\n")
                
        with open(self.datapath + '/' + motif_frame+"_terms.args",'w') as outfile:
            outfile.write(str(arguments))
                
        

    def associate_precursor_peptide(*args):
        return NotImplementedError
    
    
    
    
class Proteinlevel(Experiment):
    
    def __init__(self,experiment):
        ''''using experiment class created beforehand ensures sharing raw data, metadata, conditions and datapath
        experiment - object of Experiment class, to work properly, should have defined
        mzml_dict,metadata,conditions,datapath,quantcolumns (quantitation columns after averaging of technical replicates)
        and rawquantcolumns(raw quantification columns before averaging)
        '''
        
        self.mzml_dict = experiment.mzml_dict
        self.metadata = copy.copy(experiment.metadata)
        self.conditions = experiment.conditions
        self.datapath = experiment.datapath
        self.quancolumns = experiment.quantcolumns
        self.rawquantcolumns = experiment.rawquantcolumns
        self.tables = {}
        self.default_data = None
        
    def update_metadata(self,jsonlike):
        jsonlike = open(jsonlike,'r')
        data = json.load(jsonlike)
        metadata = data['metadata']
        self.metadata['name'] = metadata['name']
        self.metadata['columns']=metadata['columns']
        
        
        
        
        
        
        
        
        
        
        
    
        
        
    




        
    
    