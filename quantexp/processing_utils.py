# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 16:51:55 2023

@author: APirog
"""

from Bio import SeqIO
import copy


class Processing:
    
    @staticmethod
    def average_technical(rawframe,columndict):
        rawcolumns = []
        for key in columndict.keys():
            sampleruns = columndict[key]
            rawcolumns = rawcolumns + sampleruns
            rawframe[key] = rawframe[sampleruns].mean(skipna=True,axis=1)
        dataframe_av = rawframe.drop(rawcolumns,axis=1)
        
        
        return dataframe_av
    
    @staticmethod
    def filter_valid():
        
        return 0
    
    @staticmethod
    def knn_impute():
        
        return 0
    
    @staticmethod
    def normalize_median():
        
        return 0
    
    @staticmethod
    def create_razorfasta(fasta,table,outname='default'):
        '''
        Create a FASTA file containing only proteins detected in experimental table
        protein gropup ids must match protein ids from provided full proteome fasta file
        will fail if identifier from experiment is not present in fasta file
        
        '''
        
        
    def create_links(table_a,table_b,key_a,key_b,multiple_delimiter):
        '''
        I have no idea if this work
        but is probably useful to similarize tracking of precursor-peptide-protein-annotationgroup data 
        in all possible cases imaginable and faster than searching tables if there anre more than 2 in one path to traverse
        sim
    

        '''
        print(key_a)
        print(key_b)
        print(table_a.columns)
        print(table_b.columns)
        table_a_l = copy.copy(table_a)
        table_b_l = copy.copy(table_b)
        table_a[key_a] = table_a_l[key_a].str.split(multiple_delimiter)
        table_a_l = table_a_l.explode(key_a)
        table_b_l[key_b] = table_b_l[key_b].str.split(multiple_delimiter)
        table_b_l = table_b_l.explode(key_b)
        links_a_b = {}
        links_b_a = {}        
        keys_a = list(table_a_l[key_a].unique())
        keys_b = list(table_b_l[key_b].unique())

        for key in keys_a:
            print(key)
            links_a_b[key] = table_b_l[table_b_l[key_b] == key]
        for key in keys_b:
            print(key)

            links_b_a[key] = table_a_l[table_a_l[key_a] == key]
            
            
        
        
        return links_a_b,links_b_a
    
    
    def recreate_links_for_proteinlist(proteinlist,peptidetable,fasta):
        
        
        
        return 0