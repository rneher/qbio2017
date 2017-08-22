#python stage_and_call.py --raw basel_pseudo/fast5 --out basel_pseudo/calls --flowcell FLO-MIN106 --kit SQK-RAD002 --interval 60
from __future__ import print_function
import glob, os, sys, argparse, shutil, time
import numpy as np
from datetime import datetime

def call_albacore(raw_dir, out_dir, flowcell, kit, min_reads, submit_script):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    uncalled = glob.glob("%s/*fast5"%raw_dir)
    uncalled.extend(glob.glob("%s/*/*fast5"%raw_dir))
    uncalled.extend(glob.glob("%s/*/*/*fast5"%raw_dir))
    suffix=""

    if len(uncalled)>params.min_reads:
        suffix = "".join(map(chr, np.random.randint(65,89, size=20)))
        stage_dir = "%s/stage_%s"%(out_dir, suffix) 
        os.makedirs(stage_dir) 
        
        # move reads
        for f in uncalled:
            dest = stage_dir+'/'+f.replace(raw_dir,"")
            dirname = os.path.dirname(dest)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            shutil.move(f, dest)                
                
        # start albacore run on the cluster
        call = ['qsub', submit_script,  '-F', 
                '"%s %s %s %s"'%(params.flowcell, params.kit, stage_dir, out_dir+'/calls_'+suffix)]
        print(" ".join(call))
        os.system(" ".join(call))
        
    return len(uncalled), suffix


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "stage reads and call bases")
    parser.add_argument("--raw", type=str, help="directory with raw reads")
    parser.add_argument("--out", type=str, help="directory for called reads")
    parser.add_argument("--flowcell", type=str, help="flowcell")
    parser.add_argument("--kit", type=str, help="library kit")
    parser.add_argument("--script", default='rbk', type=str, help="ligation, rbk")
    parser.add_argument("--interval", type=int, default=600, help="number of second delay in the cycle")
    parser.add_argument("--min_reads", type=int, default=1000, help="minimal reads to run albacore on")
    params = parser.parse_args()
    
    raw_dir = params.raw.rstrip('/')
    out_dir = params.out.rstrip('/')

    if params.script=='rbk':
        submit_script= '/home/qbiodata/nanopore/submit_albacore_rbk.sh'
    elif params.script=='ligation':
        submit_script= '/home/qbiodata/nanopore/submit_albacore_ligation.sh'
    elif params.script=='ligation_barcoding':
        submit_script= '/home/qbiodata/nanopore/submit_albacore_ligation_barcoding.sh'

    
    staged_reads = []
    while True:
        n, suffix = call_albacore(raw_dir, out_dir, params.flowcell, params.kit, params.min_reads, submit_script)
        if n>params.min_reads:
            print(datetime.now().isoformat(), ": submitted albacore job with %d reads, directory id %s"%(n,suffix))
        else:
            print(datetime.now().isoformat(), ": insufficient number of reads, %d out of %d"%(n, params.min_reads))

        time.sleep(params.interval)
