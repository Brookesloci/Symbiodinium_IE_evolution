from sequenceAnalyzer import FastAreader
from pathlib import Path
import os
import glob

"""
Required input:
A file with a list of prefixes for each considered genome (one line per prefix). 

For each genome:
1. A genomic annotation file with the name '{Prefix}.gff' 
2. Output from "blastBack_symbiodinium.py"
all contained within a directory named with the given {Prefix}.
"""




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
                           # print(preAssembly)
                           # preList = preAssembly.split("_")
                           # print(preList[0])
                           # print(preList[1])

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
                               # print(preAssembly)
                               # preList = preAssembly.split("_")
                               # print(preList[0])
                               # print(preList[1])
    
                                Assembly = preAssembly.split("_")[0] + "_" + preAssembly.split("_")[1]
                               # print(Assembly)
                               
                                fasta = link + '/*_genomic.fna.gz'
                                annotation = link + '/*_genomic.gff.gz'
    
                                assemblyDic = {'Species': sp[0],'GenBank' : Assembly, 'genFasta' : fasta, 'genAnnotation' : annotation}
                                assemblies.append(assemblyDic)
                            except IndexError:
                                pass

            return assemblies
      




def main():
    with open('test','r') as f:
        for line in f:
            dir = line.strip()
            
            fasta = '{0}/{0}.fna'.format(dir)
            annotation = '{0}/{0}.gff'.format(dir)

            blastBackDic = {}
            ieDic = {}
            geneDic = {}
            with open('{0}/blast_back.tsv'.format(dir),'r') as f:
                for line in f:
                    sp = line.split('\t')
                    fam = sp[0]
                    scaff = sp[1]
                    start = sp[8]
                    stop = sp[9]
                    if scaff not in blastBackDic:
                        blastBackDic[scaff] = [{'fam':fam,'start':int(start),'stop':int(stop)}]

                    else:
       	       	        blastBackDic[scaff].append({'fam':fam,'start':int(start),'stop':int(stop)})

            for i in range(1,100):
                if i < 10:
                    fam = '0{0}'.format(str(i))
                else:
                    fam = str(i)
                ieFile = Path("{0}/{1}.Pass.withcoords".format(dir,fam))
                if ieFile.is_file():
                    myReaderIEs = FastAreader('{0}/{1}.Pass.withcoords'.format(dir, fam))
                    for header, sequence in myReaderIEs.readFasta():
                        coords = header.split(' ')[-1]
                        scaff = coords.split(':')[0]
                        if '+' in coords:
                            start = coords.split('+')[1].split('-')[0]
                            stop = coords.split('+')[1].split('-')[1]
                        else:
                            
                            start = coords.split('-')[1]
                            stop = coords.split('-')[2]              


                        if scaff not in ieDic:
                            ieDic[scaff]=[{'start':int(start),'stop':int(stop)}]
                        else:
                            ieDic[scaff].append({'start':int(start),'stop':int(stop)})
            with open('{0}'.format(annotation),'r') as f:
                for line in f:
                    if '#' in line:
                        pass
                    else:
                        sp = line.split('\t')
                        if sp[2] == 'mRNA':
                            if sp[0] not in geneDic:
                                if int(sp[4])>int(sp[3]):
                                    geneDic[sp[0]] = [sorted([int(sp[3]),int(sp[4])])]
                                else:
                                    geneDic[sp[0]] = [sorted([int(sp[4]),int(sp[3])])]

                            else:
       	       	       	       	if int(sp[4])>int(sp[3]):

       	       	       	            geneDic[sp[0]].append(sorted([int(sp[3]),int(sp[4])]))
                                else:
                                    geneDic[sp[0]].append(sorted([int(sp[4]),int(sp[3])]))


                            
            insideDic = {}
            outsideDic = {}
            #print(geneDic)
            for scaff in blastBackDic:
                outsideDic[scaff] = []
                insideDic[scaff] = []
                #found = False
                for coords in blastBackDic[scaff]:
                    try:
                        if scaff in ieDic:

                            found = False
                            for ie in ieDic[scaff]:
                           # print(ie,coords)
                                if coords['start'] >= ie['start'] and coords['start'] <= ie['stop']:
                                    found = True

                                    break
       	       	                elif coords['stop'] >= ie['start'] and coords['stop'] <= ie['stop']:
                                    found = True 
                                    break
                            if found == False:
                                g = False
                                for gene in geneDic[scaff]:
                                    if coords['stop'] > coords['start']:
                                        if coords['start'] >= gene[0] and coords['stop'] <= gene[1]:
                                        #print(coords['start'], gene, coords['stop'])
                                            insideDic[scaff].append(coords)
                                            g = True
                                            break
                                    else:
                                        if coords['stop'] >= gene[0] and coords['start'] <= gene[1]:
                                            insideDic[scaff].append(coords)
                                            g = True
                                            break
                                if g == False:
                                    outsideDic[scaff].append(coords)
                        else:
                            g = False
                            for gene in geneDic[scaff]:
                                if coords['stop'] > coords['start']:
                                    if coords['start'] >= gene[0] and coords['stop'] <= gene[1]:
                                        insideDic[scaff].append(coords)
                                        g = True
                                        break
                                else:
                                    if coords['stop'] >= gene[0] and coords['start'] <= gene[1]:
                                        insideDic[scaff].append(coords)
                                        g = True
                                        break
                            if g == False:
                                outsideDic[scaff].append(coords)
                    except KeyError:
                        outsideDic[scaff].append(coords)
            print(outsideDic)
           # print(insideDic)

                    
            

            fastaOutside = {}
            fastaInside = {}           
            myReaderGenome = FastAreader('{0}'.format(fasta))

            for header, sequence in myReaderGenome.readFasta():
                scaff = header.split(' ')[0]
                if scaff in outsideDic:
                    #print(scaff)
                    for coords in outsideDic[scaff]:
                        if coords['stop'] > coords['start']:
                            fastaOutside['{0}_{1}_{2}_{3}'.format(coords['fam'],scaff,str(coords['start']),str(coords['stop']))] = sequence[coords['start']:coords['stop']]
                        else:
       	       	       	    print(sequence[coords['stop']:coords['start']])

                            fastaOutside['{0}_{1}_{2}_{3}'.format(coords['fam'],scaff,str(coords['stop']),str(coords['start']))] = sequence[coords['stop']:coords['start']]
                if scaff in insideDic:

       	       	    #print(scaff)
    
                    for coords in insideDic[scaff]:
                        if coords['stop'] > coords['start']:

                            fastaInside['{0}_{1}_{2}_{3}'.format(coords['fam'],scaff,str(coords['start']),str(coords['stop']))] = sequence[coords['start']:coords['stop']]
                        else:
                            fastaInside['{0}_{1}_{3}_{2}'.format(coords['fam'],scaff,str(coords['start']),str(coords['stop']))] = sequence[coords['stop']:coords['start']]

#            print(fastaInside)

            print(fastaOutside)
            with open('{0}/outsidegenes.fa'.format(dir),'w') as f:
                for ie in fastaOutside:
                    f.write('>{0}\n{1}\n'.format(ie,fastaOutside[ie]))
            with open('{0}/insidegenes.fa'.format(dir),'w') as f:
       	       	for ie in fastaInside:
       	       	    f.write('>{0}\n{1}\n'.format(ie,fastaInside[ie]))

                                
            
            
if __name__ == "__main__":
    main() 


