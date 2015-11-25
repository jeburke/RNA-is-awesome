__author__ = 'jordanburke'

def count_read_stack_end((sortedbamfile,feature_di)):
    bamfile=HTSeq.BAM_Reader(sortedbamfile)
    print "Counting Reads in ",bamfile.filename,"..."
    cntr=0
    read_counter=collections.Counter()
    for name,end in feature_di.iteritems():
        cntr=cntr+1
        if cntr%100==0: print cntr, "features analyzed"
        for almt in bamfile[end]:
            read_counter[name]+=1