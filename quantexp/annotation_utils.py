# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 20:49:08 2023

@author: APirog
"""
import matplotlib.pyplot as pyplot
import numpy as np
import scipy.interpolate
import random
import pandas as pd
from Bio import SeqIO
import re
import itertools as it

class Annotation:
    
    
    @staticmethod
    def generate_random_peptidome(peptidetable,peptidecolumn,proteincolumn,fasta,number,len_range,make_weights = 'from_self',proteometable = None):
        '''
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
        
        '''
        
        def make_lengths(peptides,len_range,number_randoms):
            '''
            will generate number random of peptile lengths of siilar distribution to peptides
            specified in peptides file (simple one seq per row). specify lenghth range
            '''
            #peptides = open(example_peptidome).read().splitlines()
            lengths = [len(x) for x in peptides]
            counts, bins = np.histogram(lengths, bins=list(range(len_range[0]-1,len_range[1]+1,1)), density=True)
            pyplot.hist(lengths, list(range(len_range[0]-1,len_range[1]+1,1))) 
            pyplot.title('Source distribution')
            pyplot.show()
            cum_counts = np.cumsum(counts)
            bin_widths = (bins[1:] - bins[:-1])
            x = cum_counts*bin_widths
            y = bins[1:]
            inverse_density_function = scipy.interpolate.interp1d(x, y)
            b = np.zeros(number_randoms)
            for i in range(len(b)):
                u = random.uniform( x[0], x[-1] )
                b[i] = int(inverse_density_function(u))
            return b
        
        
        
        fastadata = SeqIO.parse(fasta,'fasta')
        peptides = list(peptidetable[peptidecolumn])
        peptides = [x for x in peptides if type(x) == str]
        dictionary = {}
        allproteins = []

        for sequence in fastadata:
            dictionary[str(sequence.id).split('|')[1]] = str(sequence.seq)  #dictionary of sequences ready
            allproteins.append(str(sequence.id).split('|')[1])
        if make_weights == 'from_self':
            weigths = peptidetable[proteincolumn].value_counts()
            weigths = pd.DataFrame(weigths)
            weigths = weigths.reset_index()
            weigths[proteincolumn] = weigths[proteincolumn].str.split(';')
            weigths = weigths.explode(proteincolumn)
            weigths.to_csv('test.tsv',sep='\t')
            proteins = list(weigths[proteincolumn])
            counts = list(weigths['count'])
            probablities = [x/sum(counts) for x in counts]
        elif make_weights == 'from_proteome':
            proteome = pd.read_csv(proteometable)
            proteins = list(proteome['Protein ID'])
            weigths = list(proteome['Razor Spectral Count'])
            probablities = [x/sum(weigths) for x in weigths]
            
        else:
            proteins = allproteins
        print(make_weights)
        print(probablities)
            

            
        lengths = make_lengths(peptides,len_range,number)  # 10000 properly distributed lengths ready
        sequences = []
        for i in lengths:
            i = int(i)
            if make_weights:
                protein = np.random.choice(proteins, p = probablities)
            else:
                protein = np.random.choice(proteins)
            
            
            
            
            try:
                sequence = dictionary[protein]
            except KeyError:
                print('Provided fasta does not contain indetifier '+protein)
                sequence = ''
            
            
            if len(sequence) < i+1:
                while len(sequence) < i:
                    if make_weights:
                        protein = np.random.choice(proteins, p = probablities)
                    else:
                        protein = np.random.choice(proteins)
                    try:
                        sequence = dictionary[protein]
                    except:
                        sequence = ''
                    
                    
                    
                    
            if sequence:
                startpos = random.randint(0, len(sequence) - i)
                peptide = sequence[startpos:startpos+i]
                sequences.append(peptide)
                
                
                
        lengths_out = [len(x) for x in sequences]
        pyplot.hist(lengths_out, list(range(len_range[0]-1,len_range[1]+1,1))) 
        pyplot.title('Obtained final length distribution')
        pyplot.show()
        print(sequences[0])

        return sequences
        
        

    
    @staticmethod
    def annotate_cut_motif(frame,peptide_column,modifiedpeptide_column,motif_n_extension_in,motif_n_extension_out,motif_c_extension_in,motif_c_extension_out,fasta,i_l_identical = False):
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
        
        
        
        def get_motif(pep,prot,n_in_ext,n_out_ext,c_in_ext,c_out_ext,makestring = True):
            finds = [[x.start(),x.end()] for x in re.finditer(pep,prot)]
            n_term_motifs = []
            c_term_motifs = []
            n_term_cuts = []
            c_term_cuts = []
            for item in finds:
                if item[0] < n_out_ext:
                    n_out_ext = item[0]
                if n_in_ext >= len(pep):
                    n_in_ext = len(pep)
                n_motif = prot[item[0] - n_out_ext:item[0]+n_in_ext]
                n_term_motifs.append(n_motif)
                n_term_cuts.append(n_out_ext)          # define N-term motif
                if  c_in_ext >= len(pep):
                    c_in_ext = len(pep)
                c_motif = prot[item[1] - c_in_ext:item[1]+c_out_ext]
                c_term_motifs.append(c_motif)
                c_term_cuts.append(c_in_ext)              #define C-term motif
            if makestring:
                n_term_cuts = [str(x) for x in n_term_cuts]
                c_term_cuts = [str(x) for x in c_term_cuts]
                n_term_motifs = ";".join(n_term_motifs)
                c_term_motifs = ";".join(c_term_motifs)
                n_term_cuts = ";".join(n_term_cuts)
                c_term_cuts = ";".join(c_term_cuts)
                
            return [n_term_motifs,c_term_motifs,n_term_cuts,c_term_cuts]
        
        def pep_lookup(serieslike,peptide_column,modifiedpeptide_column,motif_n_extension_in,motif_n_extension_out,motif_c_extension_in,motif_c_extension_out,fastaframe,i_l_identical = False):
            peptide = serieslike[peptide_column]
            peptide = str(peptide)
            modified_peptide = serieslike[modifiedpeptide_column]
            if i_l_identical:
                peptide_forsearch = peptide.replace("I","L")
            else:
                peptide_forsearch = peptide

            proteinframe = fastaframe[fastaframe['Sequence'].str.contains(peptide_forsearch)]
            proteins = list(proteinframe['Sequence'])
            
            allmotifs = []
            for protein in proteins:
                motifs = get_motif(peptide_forsearch,protein,motif_n_extension_in,motif_n_extension_out,motif_c_extension_in,motif_c_extension_out)
                allmotifs.append(motifs)
                
            uniquemotifs = [[],[],[],[]]        
            for motifset in allmotifs:
                for pos,n_motif in enumerate(motifset[0].split(";")):
                    if n_motif not in uniquemotifs[0]:
                        uniquemotifs[0].append(n_motif)
                        uniquemotifs[2].append(motifset[2].split(";")[pos])
                for pos,c_motif in enumerate(motifset[1].split(";")):
                    if c_motif not in uniquemotifs[1]:
                        uniquemotifs[1].append(c_motif)
                        uniquemotifs[3].append(motifset[3].split(";")[pos])
            #if "UniMod:1" in modified_peptide or "UniMod:314" in modified_peptide :
            #      print("no_modified")
            #      uniquemotifs = [[],[],[],[]]      
            uniquemotifs = [";".join(x) for x in uniquemotifs]
            print(uniquemotifs)
            motifs = pd.Series(uniquemotifs)
            
            return motifs
        
        fastadata = list(SeqIO.parse(fasta, 'fasta'))
        fastaframe = pd.DataFrame({'SeqObj':fastadata})
        if i_l_identical:
            fastaframe['Sequence'] = fastaframe['SeqObj'].apply(lambda x: str(x.seq).replace('I','L'))
        else:
            fastaframe['Sequence'] = fastaframe['SeqObj'].apply(lambda x: str(x.seq))
        
        frame[["N_term_motifs","C_term_motifs","N_term_cuts","C_term_cuts"]] = frame.apply(pep_lookup, args = [peptide_column,modifiedpeptide_column,motif_n_extension_in,motif_n_extension_out,motif_c_extension_in,motif_c_extension_out,fastaframe],axis=1)
        
        return frame
    
    
    @staticmethod
    def extract_motif_set(motif_frame,modified_peptide_column,max_motifs_c,max_motifs_n,minlength,mods_to_remove = []):
        
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
        
        
        def makefilter(serieslike,modified_peptide_column,minlength,mods_to_remove):
            modified_peptide = serieslike[modified_peptide_column]
            candidate = True
            if type(serieslike["N_term_motifs"]) == str:
                num_n = len(serieslike["N_term_motifs"].split(";"))
            else:
                num_n=0
            if type(serieslike["C_term_motifs"]) == str:
                num_c = len(serieslike["C_term_motifs"].split(";"))
            else:
                num_c=0
            if num_n > max_motifs_n:
                candidate = False
            if num_c > max_motifs_c:
                candidate = False
            if any([x in modified_peptide for x in mods_to_remove]):
                candidate = False
            return candidate
        
        def extract_fulllength(data,filtercolumn,minlength):
            n_terms = list(data["N_term_motifs"].unique())
            c_terms = list(data["C_term_motifs"].unique())
            n_terms = [x for x in n_terms if type(x) == str]
            c_terms = [x for x in c_terms if type(x) == str]
            n_terms = [x for x in n_terms if "X" not in x]
            c_terms = [x for x in c_terms if "X" not in x]
            n_terms = [x.split(";") for x in n_terms]
            c_terms = [x.split(";") for x in c_terms]
            n_terms = list(it.chain.from_iterable(n_terms))
            c_terms = list(it.chain.from_iterable(c_terms))
            n_terms = [x for x in n_terms if len(x) >= minlength]
            c_terms = [x for x in c_terms if len(x) >= minlength]
            
            return n_terms,c_terms

        
        
        
        motif_frame['for_motif_extraction'] = motif_frame.apply(makefilter,args = [modified_peptide_column,minlength,mods_to_remove],axis=1)
        motif_frame = motif_frame[motif_frame['for_motif_extraction'] == True]
        
        n_terms,c_terms = extract_fulllength(motif_frame,"for_motif_extraction", 6)
        
        return n_terms,c_terms
        
        
    
    
    @staticmethod
    def annotate_protein_position(peptidetable,peptidecolumn,fasta):
        
        '''
        Annotates peptide position in protein sequence(s)
        
        peptidetable - table with peptide sequences 
        peptidecolumn - column with unmodified peptide sequences
        fasta - fasta file
        
        function does not assume I and L as identical
        
        returns table with two additional coolumns - 'Starts','Ends' containing possibe start and end positions of peptide 
        protein sequences that contain this peptide
        
        
        '''

        def checkloc_total_fasta(peptide,fastaf):
            print(peptide)

            
            try:
                peptide = peptide.replace('I','L')
            except:
                peptide = 'ZZZ'
            if peptide == '':
                peptide = 'ZZZ'
            starts = []
            ends = []
            
            proteinframe = fastaf[fastaf['Sequence'].str.contains(peptide)]
            proteins = list(proteinframe['Sequence'])
            for sequence in proteins:
                if peptide in sequence:
                    startpos = sequence.find(peptide)
                    endpos = (startpos+len(peptide)) - 1
                    perc_start = startpos/len(sequence)
                    perc_end = endpos/len(sequence)
                    starts.append(perc_start)
                    ends.append(perc_end)
                    if perc_end < 0:
                        print(startpos)
                        print(endpos)
                        print(perc_start)
                        print(perc_end)
                        print('tupeptyd'+peptide+'tupeptyd')
                        print(sequence)
            starts = [str(x) for x in starts]
            ends = [str(x) for x in ends]
            starts = ';'.join(starts)
            ends = ';'.join(ends)
            returnval = pd.Series([starts,ends])
            
            return returnval
        
        fastadata = list(SeqIO.parse(fasta, 'fasta'))
        fastaframe = pd.DataFrame({'SeqObj':fastadata})
        fastaframe['Sequence'] = fastaframe['SeqObj'].apply(lambda x: str(x.seq).replace('I','L'))
        peptidetable[['Starts','Ends']] = peptidetable[peptidecolumn].apply(checkloc_total_fasta, args = [fastaframe])
        
        
        return peptidetable
    
    @staticmethod
    def annotate_go_term():
        
        
        return 0