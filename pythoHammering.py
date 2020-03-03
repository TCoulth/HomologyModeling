#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:57:44 2020

@author: TimResearch
"""

import os
import pandas as pd
import pypdb

DBfolder = '/share/siegellab/tcoulther/DBs/'
RunningFiles = '/share/siegellab/fell/homology-modeling/'
HomeFolder = os.getcwd()


def ConvertCatCST(CatResFile):
    with open(CatResFile) as file:
        newfilecontents = []
        for line in file:
            newnames = []
            pdb , residues = line.split()
            catsite = residues.split(',')
            for item in catsite:
                oldnum = item[3:]
                newnum = MakePDBDict(HomeFolder+'/Templates/Maps/'+pdb+'.map')[oldnum]
                newnames.append(item[:3]+newnum)
            residues = ','.join(newnames)
            newfilecontents.append(pdb+' '+residues)
    with open(CatResFile[:-12]+'.catres' , 'a') as newfile:
        for line in newfilecontents:
            newfile.write(line+'\n')
    os.system('rm '+CatResFile)
        
                
def MakePDBDict(mapfile):
    numbering = {}
    with open(mapfile) as file:
        for line in file:
            newnum , oldnum = line.split()
            numbering[oldnum] = newnum
    return numbering
            
def MakeUniq(columns):
    seenitem = set()
    for item in columns:
        newitem = item
        while newitem in seenitem:
            newitem = "{}_{}".format(item, 'domain')
        yield newitem
        seenitem.add(newitem)

def MakeTemplateList(fasta):
    counter = 0
    DataList = []
    for file in os.listdir('.'):
        if file.endswith('.table'):
            source = file[-9:-6]
            for line in open(file):
                counter += 1
                if counter == 2:
                    headerstring = line
                    headerstring = headerstring.replace('target name' , 'target-name')
                    headerstring = headerstring.replace('query name' , 'query-name')
                    headerstring = headerstring.replace('description of target' , 'description-of-target')
                    headerstring = headerstring[1:].split()
                    headerstring.append('UniProtID')
                    headerstring.append('Source')
                    headerstring = list(MakeUniq(headerstring))
                    
                elif counter >=4:
                    if line[0] == '#':
                        pass
                    else:
                        line2 = line.split()
                        description = ' '.join(line2[18:])
                        line3 = line2[:18]
                        line3.append(description)
                        if source == 'CSA':
                            line3[0]=line3[0][:6].lower().replace(':' , '_')
                        PDBinfo = pypdb.get_all_info(line3[0][:4])['polymer']
                        if isinstance(PDBinfo , list):
                            UniProtID = PDBinfo[0]['macroMolecule']['accession']['@id']
                        else:
                            UniProtID = PDBinfo['macroMolecule']['accession']['@id']
                        line3.append(UniProtID)
                        line3.append(source)
                        DataList.append(line3)
    df = pd.DataFrame(DataList , columns = headerstring )
    df['E-value'] = df['E-value'].astype(float)
    df.sort_values(by = ['Source','E-value'] , inplace = True)
    df.to_csv('HmmerInfo/CSV/'+fasta+'_fullhmmer.csv')
    df.drop_duplicates(subset = 'UniProtID' , inplace = True)
    df.to_csv('HmmerInfo/CSV/'+fasta+'_noduphmmer.csv')
    targetlist = df['target-name'].tolist()
    os.system('mv '+fasta+'.MCSA.table HmmerInfo/MCSAhm')
    os.system('mv '+fasta+'.PDB.table HmmerInfo/PDBhm')
    return targetlist

def HammerTime(fastafile , DB):
    DBname = DB[:-3]
    fasta = fastafile[:-6]
    output = fasta+'.'+DBname+'.table'
    os.system('phmmer --tblout '+output+' Fastas/'+fastafile+' '+DBfolder+DB)

    
        
os.mkdir('HmmerInfo')
os.mkdir('HmmerInfo/CSV')
os.mkdir('HmmerInfo/MCSAhm')
os.mkdir('HmmerInfo/PDBhm') 
os.mkdir('Templates')
os.mkdir('Templates/Lists')
os.mkdir('Templates/Fastas')
os.mkdir('Templates/PDBs')  
os.mkdir('Templates/Maps')
os.mkdir('Run')

PDBList= {}
for fastafile in os.listdir('Fastas/'):
    fasta = fastafile[:-6]
    HammerTime(fastafile , 'MCSA.fa')
    HammerTime(fastafile , 'PDB.fa')
    PDBList[fasta] = MakeTemplateList(fasta)

print(PDBList)
for entry in PDBList:
    for template in PDBList[entry]:
        pdb , ch = template.split('_')
        ch = ch.upper()
        PDB = pdb.upper()
        if pdb+'.pdb' in os.listdir('Templates/PDBs/'):
            pass
        else:
            os.system('grabchainmap '+PDB+' '+ch)
            os.system('mv '+PDB+'.pdb Templates/PDBs/'+pdb+'.pdb')
            with open(pdb+'.fasta' , 'a') as newfasta:
                newfasta.write('>'+pdb+'\n')
                for line in open(PDB+'.fasta'):                    
                    if line.startswith('>'):
                        pass
                    else:
                        newfasta.write(line)
            os.system('mv '+pdb+'.fasta Templates/Fastas/'+pdb+'.fasta')
            os.system('rm '+PDB+'.fasta')
            os.system('mv '+PDB+'.map Templates/Maps/'+pdb+'.map')
                
for targetfasta in os.listdir('Fastas/'):
    target = targetfasta[:-6]
    os.mkdir('Run/'+target)
    os.system('cp Fastas/'+targetfasta+' Run/'+target+'/')
    for entry in PDBList[target]:        
        templatepdb , ch = entry.split('_')
        os.system('cp Templates/Fastas/'+templatepdb+'.fasta Run/'+target+'/')
        os.system('cp Templates/PDBs/'+templatepdb+'.pdb Run/'+target+'/')
        os.system('grep "'+templatepdb+'" '+DBfolder+'/MCSAcst/MCSACatRes.list >>  Run/'+target+'/'+target+'.orig.catres')
    os.chdir('Run/'+target)
    ConvertCatCST(target+'.orig.catres')
#    os.system('Firststep '+target)
    os.chdir('../..')
    
#for targetfasta in os.listdir('Fastas/'):
#    target = targetfasta[:-6]
#    os.system('python '+RunningFiles+'HM_ii_StructuralConst.py')
#    
#for targetfasta in os.listdir('Fastas/'):
#    target = targetfasta[:-6]
#    os.system('python '+RunningFiles+'RosettaEnzCM.py -c '+target+'catres -a '+target+'.family_fasta -n '+target)
#    os.system('cat '+fasta+'.dist_csts >> '+fasta+'/'+fasta+'alignment.grishin.dist_csts')
#
#for targetfasta in os.listdir('Fastas/'):
#    target = targetfasta[:-6]
#    os.system('python '+RunningFiles+'HM_iii_Modeling.py')
  

    
    
    