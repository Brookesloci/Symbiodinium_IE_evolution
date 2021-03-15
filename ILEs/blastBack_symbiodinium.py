
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 09:40:21 2019

@author: lgozasht
"""


""" 

This script generates a multiple sequence alignment of all sequences for each ILE family using MAFFT, 
then uses a positional frequency matrix to construct a consensus sequence for each family
and blasts the consensus back to its corresponding reference.

Required input:

A file with a list of prefixes for each considered genome (one line per prefix). 

For each genome:
1. A genomic fasta file with the name '{Prefix}.fna' 
2. A genomic annotation file with the name '{Prefix}.gff' 
3. A file containing ILEs the name '{Prefix}.Pass.withcoords' and with the header format provided in the .Pass.withcoords example
all contained within a directory named with the given {Prefix}.

"""



import os
from sequenceAnalyzer import FastAreader
import glob
from pathlib import Path



class csvReader():


    def __init__(self, dataFile = ''):
        """
        tialize blasthit as the blast output file.

        """
        self.dataFile = dataFile

    def csv(self):
            """
            Return results, a dictionary containing data extracted from blasthit for each query.


            Edit this to remove alignments that don't align for 80% of sequence length ie. if queryLength/alignmentLength > .8:......
            """
            assemblies = []
            with open(self.dataFile,'r') as f:
                for line in f:
                    line = line.rstrip()
                    sp = line.split(',')
                    #print(len(sp))
                    
                   


                    if len(sp) == 16:
                        #print('what')
                        refLink = sp[15]
                        refLink = refLink.strip('\"')
                        refLink = refLink.strip('\"')
                        genLink = sp[14]
                        genLink = genLink.strip('\"')
                        genLink = genLink.strip('\"')

                        try:
                            #print(link)
                            refPreAssembly = refLink.split("/")[-1]

                            refAssembly = refPreAssembly.split("_")[0] + "_" + refPreAssembly.split("_")[1]
                           # print(Assembly)
                            refFasta = refLink + '/*_genomic.fna.gz'
                            refAnnotation = refLink + '/*_genomic.gff.gz'
                            genPreAssembly = genLink.split("/")[-1]
                           # print(preAssembly)
                           # preList = preAssembly.split("_")
                           # print(preList[0])
                           # print(preList[1])

                            genAssembly = genPreAssembly.split("_")[0] + "_" + genPreAssembly.split("_")[1]
                           # print(Assembly)
                            genFasta = genLink + '/*_genomic.fna.gz'
                            genAnnotation = genLink + '/*_genomic.gff.gz'

                            assemblyDic = {'Species': sp[0],'RefSeq' : refAssembly, 'refFasta' :refFasta , 'refAnnotation' : refAnnotation, 'GenBank' : genAssembly, 'genFasta' :genFasta , 'genAnnotation' : genAnnotation}
                            assemblies.append(assemblyDic)
                        except IndexError:
                        
                            print("Fuck you")
                            link = sp[14]
                            link = link.strip('\"')
                            link = link.strip('\"')
                            try:
                                #print(link)
                                preAssembly = link.split("/")[-1]
    
                                Assembly = preAssembly.split("_")[0] + "_" + preAssembly.split("_")[1]
                               # print(Assembly)
                               
                                fasta = link + '/*_genomic.fna.gz'
                                annotation = link + '/*_genomic.gff.gz'
    
                                assemblyDic = {'Species': sp[0],'GenBank' : Assembly, 'genFasta' : fasta, 'genAnnotation' : annotation}
                                assemblies.append(assemblyDic)
                            except IndexError:
                                pass

            return assemblies

class Consensus():
    def __init__(self, NAME):
        self.NAME = NAME
        self.consensusDic = {}
        self.famDic = {}       
    def msa(self):
                 
        myReaderIEs= FastAreader('{0}/{0}.Pass.withcoords'.format(self.NAME))
        for header, sequence in myReaderIEs.readFasta():
            try:
                print(header.split('Group+')[1])
                fam = header.split('Group+')[1].split(' ')[0]
            except IndexError:
            
                fam = header.split('Group-')[1].split(' ')[0]

            print(fam)
            self.famDic[fam] = ''       
            with open('{0}/{1}_IEs.fa'.format(self.NAME, str(fam)),'a') as f:
                f.write('>{0}\n{1}\n'.format(header,sequence[20:-20]))
        
        for fam in self.famDic:
            os.system('mafft {0}/{1}_IEs.fa > {0}/{1}_IEs_msa.fa'.format(self.NAME, str(fam)))
            self.consensusDic[fam] = self.consensus(fam)
        #print(self.consensusDic)
        return self.consensusDic


    def consensus(self,x):
        myReaderIEs= FastAreader('{0}/{1}_IEs_msa.fa'.format(self.NAME, x))
        indexDic = {}
        total = 0
        consensus = ''
        for header, sequence in myReaderIEs.readFasta():
            #print(len(sequence))
            total += 1
            for i in range(0,len(sequence)):
                if i not in indexDic:
                    indexDic[i] = {'a':0,'t':0,'g':0,'c':0,'-':0,'n':0}
                if sequence[i] not in indexDic[i]:
                    indexDic[i][sequence[i]] = 0   
                indexDic[i][sequence[i]] += 1
        for index in indexDic:
            for nuc in indexDic[index]:
                if float(indexDic[index][nuc])/float(sum(indexDic[index].values())) > .5 and nuc != '-':
                    consensus += nuc.upper()
        #print(len(consensus))
        return consensus




def main():
    with open('test','r') as f:
        for line in f:
            try: 
                dir = line.strip()
            
                fasta = '{0}.fna'.format(dir)
                os.system('mkdir {0}'.format(dir))
                print(fasta)

                myAlignments = Consensus(dir)
                consensusDic = myAlignments.msa()
                with open('{0}/IEconsensi.fa'.format(dir),'w') as f:
                    for header in consensusDic:
       	       	        f.write('>{0}\n{1}\n'.format(header, consensusDic[header]))
       	        os.system('makeblastdb -dbtype nucl -in {0}/{1} -out {0}/genomeDB'.format(dir,fasta))                    
                os.system('blastn  -db {0}/genomeDB -query {0}/IEconsensi.fa -outfmt 6 -perc_identity 80 -out {0}/blast_back.tsv'.format(dir))
            except FileNotFoundError:
                print(dir)                  
  
            
            
if __name__ == "__main__":
    main() 




