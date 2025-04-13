from multiprocessing import Pool
from datetime import datetime
import time
import resource
import pysam
import re
import uuid
import os
import gzip

from optparse import OptionParser
import sys, gzip
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-p","--prefix",action = 'store',type = 'string',dest = 'prefix',help = "")
parser.add_option("-t","--threadN",action = 'store',type = 'int',dest = 'threadN',help = "")
(opt, args) = parser.parse_args()
if opt.prefix == None or opt.threadN == None:
    print('Basic usage')
    print('')
    print('     python extract_read2.py -t 24 -p test')
    print('')
    sys.exit()

prefix = opt.prefix
batchN = opt.threadN

#unique_id = str(uuid.uuid4())
unique_id = '3f580251-90ad-4172-ac87-38edefb8643d'
tmpDir = f'tmp/{unique_id}'

###########################################################################################
def log(msg):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {msg}", flush=True)

def get_softclip_baseN(read):
    if read.cigartuples is None:
        return 0
    return sum(length for op, length in read.cigartuples if op == 4)

def filter_read(batchIDX):
    readIDX = -1

    fin = pysam.AlignmentFile(prefix + '.bam', "r")
    fout_pass = pysam.AlignmentFile(f'{tmpDir}/{prefix}_{batchIDX:05d}.pass.bam', "wb", template=fin)
    fout_fail = pysam.AlignmentFile(f'{tmpDir}/{prefix}_{batchIDX:05d}.fail.bam', "wb", template=fin)


    totalN = 0
    passN = 0
    for read in fin:
        if read.is_secondary or read.is_supplementary: continue

        readIDX += 1
        if readIDX%batchN != batchIDX: continue
        totalN += 1

        #unmapped
        if read.is_unmapped:
            fout_fail.write(read)
            continue
            
        #softclip
        softclip_baseN = get_softclip_baseN(read)
        softclip_rate = softclip_baseN / read.query_length
        if softclip_rate > 0.01:
            fout_fail.write(read)
            continue

        #mismatch
        mismatchN = read.get_tag("NM")
        mismatch_rate = mismatchN / read.query_length
        if mismatch_rate > 0.01:
            fout_fail.write(read)
            continue
        
        fout_pass.write(read)
        passN += 1
    
    fout_pass.close()
    fout_fail.close()
    fin.close()

    return passN, totalN

start_time = time.time()

#make tmpDir
if not os.path.exists(tmpDir):
    os.makedirs(tmpDir)

log(f'read file start: {prefix}.bam')
with Pool(processes=batchN) as pool:
    result_LIST = pool.map(filter_read, range(batchN))
log(f'read file done:  {prefix}.bam')

#show statistics
passN = 0
totalN = 0
for _passN, _totalN in result_LIST:
    passN += _passN
    totalN += _totalN
log(f'filter result - totalN:{totalN}, passN:{passN}, rate:{passN / totalN * 100 : 0.6f}%')

#pass
log(f'write pass file start: {prefix}.pass.bam')
fin = pysam.AlignmentFile(prefix + '.bam', "r")
fout = pysam.AlignmentFile(f'/{prefix}.pass.bam', "wb", template=fin)
fin.close()
for batchIDX in range(batchN):
    fin = gzip.open(f'{tmpDir}/{prefix}_{batchIDX:05d}.pass.bam', 'rb')
    for read in fin:
        fout.write(read)
    fin.close()
fout.close()
log(f'write pass file done:  {prefix}.pass.bam')

#fail
log(f'write pass file start: {prefix}.fail.bam')
fin = pysam.AlignmentFile(prefix + '.bam', "r")
fout = pysam.AlignmentFile(f'/{prefix}.fail.bam', "wb", template=fin)
fin.close()
for batchIDX in range(batchN):
    fin = gzip.open(f'{tmpDir}/{prefix}_{batchIDX:05d}.fail.bam', 'rb')
    for read in fin:
        fout.write(read)
    fin.close()
fout.close()
log(f'write pass file done:  {prefix}.fail.bam')

end_time = time.time()
elapsed = end_time - start_time
peak_memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
peak_memory_gb = peak_memory_kb / 1024 / 1024
log(f"Program finished - Elapsed time: {elapsed:.2f} sec, Peak memory usage: {peak_memory_gb:.2f} GB")