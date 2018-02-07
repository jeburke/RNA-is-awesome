import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome')
import GeneTools as GT

bam = sys.argv[1]

GT.PE_fragment_size(bam)