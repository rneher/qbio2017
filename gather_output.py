import glob, argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')
ntop = 10
barcode_dir = ['barcode%02d'%b for b in range(1,13)] + ['unclassified']

def concatenate_reads_barcoded(mask='', out='.'):
    fig = plt.figure()
    ax1 = plt.subplot(111)
    fig = plt.figure()
    ax2 = plt.subplot(111)
    fnames = glob.glob(mask.rstrip('/')+'/workspace/*/*fastq')
    top_reads = [[0,0, '']]

    stats_file = open(out+'/run_statistic.txt', 'w')
    total_size = 0
    total_reads = 0
    for dname in barcode_dir:
        read_lengths = []
        print('####processing barcode %s'%dname)
        with open(out+'/%s.fastq'%dname, 'w') as ofile:
            for fastq_file in fnames:
                if dname in fastq_file:
                    print('adding file %s'%fastq_file)
                    for read in SeqIO.parse(fastq_file, 'fastq'):
                        SeqIO.write(read, ofile, 'fastq')
                        l_read = len(read.seq)
                        read_lengths.append(l_read)
                        if l_read>top_reads[0][0]:
                            top_reads.append([l_read, read, dname])
                            top_reads.sort(key=lambda x:x[0])
                            top_reads = top_reads[-ntop:]
                            

        read_lengths = np.array(read_lengths)
        total_size+=read_lengths.sum()
        total_reads+=read_lengths.shape[0]
        ax1.plot(sorted(read_lengths), len(read_lengths) - np.arange(len(read_lengths)), label=dname)
        ax2.plot(sorted(read_lengths), sum(read_lengths) - np.cumsum(sorted(read_lengths)), label='Total %s=%4.2fMb'%(dname, read_lengths.sum()/1000000))
        
        stats_file.write("\n####  %s\ntotal output: %1.1e\n\n"%(dname, np.sum(read_lengths)))
        for cutoff in [1000, 2000, 3000, 10000, 20000, 30000, 50000, 100000]:
            stats_file.write("sequence in reads above %dbp: %1.1e\n"%(cutoff, np.sum(read_lengths*(read_lengths>cutoff))))

    plt.figure(1)
    plt.ylabel('#reads longer than x')
    plt.xlabel('read length')
    plt.xscale('log')
    plt.legend()
    plt.savefig(out+'/read_length.png')

    plt.figure(2)
    plt.ylabel('#sequence in reads longer than x')
    plt.xlabel('read length')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig(out+'/total_sequence.png')

    stats_file.write("\n\n #### TOTAL OUTPUT: %1.2e"%total_size)
    stats_file.write("\n\n #### TOTAL #READS: %1.2e"%total_reads)
    stats_file.close()
    return top_reads


def concatenate_reads(mask='', out='.'):

    fnames = glob.glob(mask.rstrip('/')+'/workspace/*fastq')
    stats_file = open(out+'/run_statistic.txt', 'w')
    total_size = 0
    total_reads = 0
    read_lengths = []
    top_reads = [[0,0, '']]
    with open(out+'/all_reads.fastq', 'w') as ofile:
        for fastq_file in fnames:
            print('adding file %s'%fastq_file)
            for read in SeqIO.parse(fastq_file, 'fastq'):
                SeqIO.write(read, ofile, 'fastq')
                read_lengths.append(len(read.seq))
                l_read = len(read.seq)
                read_lengths.append(l_read)
                if l_read>top_reads[0][0]:
                    top_reads.append([l_read, read, ''])
                    top_reads.sort(key=lambda x:x[0])
                    top_reads = top_reads[-ntop:]
                    
        read_lengths = np.array(read_lengths)
        total_size+=read_lengths.sum()
        total_reads+=read_lengths.shape[0]

    stats_file.write("\n\n #### TOTAL OUTPUT: %1.2e"%total_size)
    stats_file.write("\n #### TOTAL #READS: %1.2e\n\n"%total_reads)
    for cutoff in [1000, 2000, 3000, 10000, 20000, 30000, 50000, 100000]:
        stats_file.write("sequence in reads above %dbp: %1.1e\n"%(cutoff, np.sum(read_lengths*(read_lengths>cutoff))))

    plt.figure()
    plt.plot(sorted(read_lengths), len(read_lengths) - np.arange(len(read_lengths)))
    plt.ylabel('#reads longer than x')
    plt.xlabel('read length')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(out+'/read_length.png')

    plt.figure()
    plt.plot(sorted(read_lengths), sum(read_lengths) - np.cumsum(sorted(read_lengths)), label='Total=%4.2fMb'%(total_size/1000000))
    plt.ylabel('#sequence in reads longer than x')
    plt.xlabel('read length')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig(out+'/total_sequence.png')
    stats_file.close()

    return top_reads

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "gather fastq from different albacore runs")
    parser.add_argument("--mask", type=str, help="mask to select for input directory")
    parser.add_argument("--out", default='.', type=str, help="output directory")
    parser.add_argument("--barcoded", action="store_true", default=False, help="barcoded library")
    params = parser.parse_args()

    mask=params.mask
    barcoded=params.barcoded
    if barcoded:
        top_reads = concatenate_reads_barcoded(mask=mask, out=params.out)
    else:
        top_reads = concatenate_reads(mask=mask, out=params.out)
        

    with open("longest_reads.fasta", 'w') as seq_file, open('longest_reads.txt', 'w') as stats_file:
        for (l, read, dname) in top_reads:
            SeqIO.write(read, seq_file, 'fasta')
            seq_array = np.array(read.seq)
            nc = [np.mean(seq_array==a) for a in 'ACGT']
            stats_file.write('%dbp\t\t%1.3f A\t%1.3f C\t%1.3f G\t%1.3f T\t%s\n'%(len(read.seq), nc[0], nc[1], nc[2], nc[3], dname))
