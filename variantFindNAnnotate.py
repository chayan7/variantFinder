#Author: Chayan Kumar Saha
#Description: Variant calling and annotation from GISAID fasta files


import argparse
import glob
import os
import sys
import os.path
import queue as Queue
from queue import Queue, Empty
import threading

usage = '''  Description: Variant Calling and Annotation'''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument(
    "-d", "--datadir", help=" Path for GISAID fasta files, Default directory is './Data' which is the same directory where the script is located or running from. ")
parser.add_argument("-r", "--refGenomeDir",
                    help=" Path for reference genome directory, use genomic fna sequence, Default directory is './referenceGenome' ")
parser.add_argument(
    "-c", "--cpu", help="Maximum number of parallel CPU workers to use for multithreads. ")
parser.add_argument("-v", "--version", action="version",
                    version='%(prog)s 1.0.0')
args = parser.parse_args()
parser.parse_args()

if args.cpu:
    if int(args.cpu) > 0:
        core = int(args.cpu)
    else:
        print('Please use number eg, 1,2...')
        sys.exit()

def worker_func():
    while not stopped.is_set():
        try:
            # use the get_nowait() method for retrieving a queued item to
            # prevent the thread from blocking when the queue is empty
            com = q.get_nowait()
        except Empty:
            continue
        try:
            os.system(com)
        except Exception as e:
            print('Error running command', str(e))
        finally:
            q.task_done()


current_directory = os.getcwd()+'/'

data_directory = ''

if not args.datadir:
    if os.path.exists(current_directory+'Data/'):
        data_directory = current_directory+'Data/'
        print('Data path : ', data_directory, '\n')
else:
    if os.path.isdir(args.datadir):
        if args.datadir[-1] == '/':
            data_directory = args.datadir
            print('Data path : ', data_directory, '\n')
        else:
            data_directory = args.datadir+'/'
            print('Data path : ', data_directory, '\n')
    else:
        print('Data path is missing, give the directory for GISAID datafiles', '\n')
        sys.exit()


reference_directory = ''

if not args.refGenomeDir:
    if os.path.exists(current_directory+'referenceGenome/'):
        reference_directory = current_directory+'referenceGenome/'
        print('Reference Genome Path : ', reference_directory, '\n')
else:
    if os.path.isdir(args.refGenomeDir):
        if args.refGenomeDir[-1] == '/':
            reference_directory = args.refGenomeDir
            print('Reference Genome Path : ', reference_directory, '\n')
        else:
            reference_directory = args.refGenomeDir+'/'
            print('Reference Genome Path : ', reference_directory, '\n')
    else:
        print('Reference path is missing, try again please', '\n')
        sys.exit()


dataList = []
for files in (glob.glob(data_directory+'*.f*a')):
    dataList.append(files)

refGenome = []
for file in (glob.glob(reference_directory+'*.f*a')):
    refGenome.append(file)

if len(refGenome) > 1:
    print('Multiple reference genome files detected. Please, use one reference genome. ')
    sys.exit()

#refFastaIndex command

refIndexcommand = "/home/chayan/bin/samtools-1.9/samtools faidx %s"%(refGenome[0])
os.system(refIndexcommand)
print('\n', '-- Reference Index created.')
#execute

#samfiles creation with minimap command

def renameFile(fileName, newExt):
    prevExt=fileName.split('.')[-1]
    renamedFile=fileName.replace(prevExt, newExt)
    return renamedFile

def fileCounter(dir, Ext):
    fList=[]
    for files in (glob.glob(dir+'*.'+Ext)):
        fList.append(files)
    if fList:
        return len(fList)
    else:
        return 0


samfiles=[]
if not os.path.exists(current_directory+'SAM_files/'):
    os.makedirs(current_directory+'SAM_files/')
    minimapCom=[]
    for fastafiles in dataList: #trialData/GCF_009858895.2.fasta
        samFile=current_directory+'SAM_files/'+renameFile(fastafiles.split('/')[-1], 'sam')
        samfiles.append(samFile)
        mCom="/home/chayan/bin/minimap2/minimap2 -ax asm5 -t %s %s %s > %s"%(round(core/len(dataList)), refGenome[0], fastafiles, samFile)
        minimapCom.append(mCom)
    stopped = threading.Event()
    q = Queue()
    print('\n-- Processing : SAM_files creation with ' + (str(len(minimapCom)) +
          ' tasks in thread queue with ' + str(core)) + ' thread limit')
    for item in minimapCom:
        q.put(item)
    for x in range(core):
        t = threading.Thread(target=worker_func)
        # t.daemon = True #Enable to run threads as daemons
        t.start()
    q.join()       # block until all tasks are done
    stopped.set()
    print('## Process SAM_files creation Done', '\n')
else:
    if fileCounter(current_directory+'SAM_files/', 'sam')==len(dataList):
        pass


#Bam file create

bamfiles=[]
if not os.path.exists(current_directory+'BAM_files/'):
    os.makedirs(current_directory+'BAM_files/')
    #samtools view -bS ${fasta}.sam > ${fasta}.bam
    #samtools sort -o ${fasta}.sorted.bam ${fasta}.bam
    samtoolViewCom=[]
    samtoolSortCom=[]
    for samitems in samfiles: #trialData/GCF_009858895.2.fasta
        bamFile=current_directory+'BAM_files/'+renameFile(samitems.split('/')[-1], 'bam')
        sortedBam=current_directory+'BAM_files/'+renameFile(samitems.split('/')[-1], 'sorted.bam')
        bamfiles.append(sortedBam)
        svCom="/home/chayan/bin/samtools-1.9/samtools view -bS %s > %s"%(samitems, bamFile)
        ssCom="/home/chayan/bin/samtools-1.9/samtools sort -o %s %s"%(sortedBam, bamFile)
        samtoolViewCom.append(svCom)
        samtoolSortCom.append(ssCom)

    stopped = threading.Event()
    q = Queue()
    print('\n-- Processing : BAM_files creation with ' + (str(len(samtoolViewCom)) +
          ' tasks in thread queue with ' + str(core)) + ' thread limit')
    for item in samtoolViewCom:
        q.put(item)
    for x in range(core):
        t = threading.Thread(target=worker_func)
        # t.daemon = True #Enable to run threads as daemons
        t.start()
    q.join()       # block until all tasks are done
    stopped.set()
    print('## Process BAM_files creation Done', '\n')

    stopped = threading.Event()
    q = Queue()
    print('\n-- Processing : BAM_files sorting with ' + (str(len(samtoolSortCom)) +
          ' tasks in thread queue with ' + str(core)) + ' thread limit')
    for item in samtoolSortCom:
        q.put(item)
    for x in range(core):
        t = threading.Thread(target=worker_func)
        # t.daemon = True #Enable to run threads as daemons
        t.start()
    q.join()       # block until all tasks are done
    stopped.set()
    print('## Process BAM_files sorting Done', '\n')

else:
    if fileCounter(current_directory+'SAM_files/', 'sorted.bam')==len(dataList):
        pass



vcf_files=[]
if not os.path.exists(current_directory+'VCF_files/'):
    os.makedirs(current_directory+'VCF_files/')
    vcfCom=[]
    for bamItems in (glob.glob(current_directory+'BAM_files/'+'*sorted.bam')):
        vcfFile=current_directory+'VCF_files/'+renameFile(bamItems.split('/')[-1].replace('.sorted',''), 'vcf')
        filtervcfFile=current_directory+'VCF_files/'+renameFile(bamItems.split('/')[-1].replace('.sorted',''), 'filtered.vcf')
        vcf_files.append(vcfFile)
        vCom="/usr/bin/freebayes -f %s -C 1 %s > %s"%(refGenome[0], bamItems, vcfFile)
        vcfCom.append(vCom)

    stopped = threading.Event()
    q = Queue()
    print('\n-- Processing : VCF_files creation with ' + (str(len(vcfCom)) +
          ' tasks in thread queue with ' + str(core)) + ' thread limit')
    for item in vcfCom:
        q.put(item)
    for x in range(core):
        t = threading.Thread(target=worker_func)
        t.start()
    q.join()       # block until all tasks are done
    stopped.set()
    print('## Process VCF_files creation Done', '\n')
else:
    if fileCounter(current_directory+'VCF_files/', '.vcf')==len(dataList):
        pass

annotated_vcf_files=[]
if not os.path.exists(current_directory+'Annotated_VCF_files/'):
    os.makedirs(current_directory+'Annotated_VCF_files/')
    annotateCom=[]
    for vcfs in vcf_files:
        ann_vcfFile=current_directory+'Annotated_VCF_files/'+renameFile(vcfs.split('/')[-1], 'txt')
        annotated_vcf_files.append(ann_vcfFile)
        annCom="/home/chayan/bin/snpEff/scripts/snpEff NC_045512.2 %s | grep -v '#' | cut -f 1,2,8 | sed 's/|/\t/g' | cut -f 1,2,4,5,6,12,13 > %s"%(vcfs, ann_vcfFile)
        annotateCom.append(annCom)

    stopped = threading.Event()
    q = Queue()
    print('\n-- Processing : VCF_files Annotation with ' + (str(len(annotateCom)) +
          ' tasks in thread queue with ' + str(core)) + ' thread limit')
    for item in annotateCom:
        q.put(item)
    for x in range(core):
        t = threading.Thread(target=worker_func)
        # t.daemon = True #Enable to run threads as daemons
        t.start()
    q.join()       # block until all tasks are done
    stopped.set()
    print('## Process VCF_files Annotation Done', '\n')
else:
    if fileCounter(current_directory+'Annotated_VCF_files/', '.txt')==len(dataList):
        pass


mutationORF=set()
mutation=set()
sampleInfo_dict={}
for items in (glob.glob(current_directory+'Annotated_VCF_files/'+'*.txt')):
    sampleID=items.split('/')[-1].replace('.txt','')
    infoList=[]
    with open(items, 'r') as fileIn:
        for line in fileIn:
            Line=line.rstrip().split('\t')
            mutationORF.add(Line[4])
            mutation.add(Line[5])
            infoList.append(sampleID+'\t'+line.rstrip())
    sampleInfo_dict[sampleID]=infoList

with open('Merged_Variants.txt','w') as tOut:
    print('#GISAID_ids', 'REF', 'POS', 'Variant_Type', 'Level', 'Gene', 'Ntd_Alteration', 'Prot_Alteration', sep='\t', file=tOut)
    for IDs in sorted(sampleInfo_dict):
        for items in sampleInfo_dict[IDs]:
            print(items, file=tOut)
