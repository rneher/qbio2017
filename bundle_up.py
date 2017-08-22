import glob, gzip, os

for a in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    print('tar czf base_called_reads_%s.tar.gz calls_%s*'%(a,a))
    os.system('tar czf base_called_reads_%s.tar.gz calls_%s*'%(a,a))


