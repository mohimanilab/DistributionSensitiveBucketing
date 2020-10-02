import sys
import numpy as np
from collections import defaultdict
import argparse

class GroundTruth:
    def __init__(self, query, q_start, q_end, t_start, t_end, direction):
        self.query = query
        self.q_start = q_start
        self.q_end = q_end
        self.t_start = t_start
        self.t_end = t_end
        self.length = abs(t_end - t_start)
        self.dir = direction
    
    def __str__(self):
        return "{} {} {} {} {} {}".format(self.query, self.dir,
            self.q_start, self.q_end, self.t_start, self.t_end)

# analyze mapped information for a query
# returns all significant alignments found (>=50bp)
def analyze_mapped(mapped):
    alignments = []
    cur_q, cur_t = [-1,-1], [-1,-1] # cur_alignments
    
    # mapped is sorted by second value (bp position on genome)
    for m in mapped:
        q_pos, t_pos = m

        # initialization
        if cur_q[0] == -1:
            cur_q = [q_pos, q_pos]
        if cur_t[0] == -1:
            cur_t = [t_pos, t_pos]
        
        # since mapped is sorted by second position, we check:
        # 1) if t_pos is within a tolerable range from last t_pos;
        # 2) if (1) is true, if q_pos is within a tolerable range from last
        # q_pos
        #
        # if (1) is false, see if cur_alignment is significant; if so, add it
        # to alignments; clear cur_alignment and start over from q_pos, t_pos
        # if (1) is true, (2) is false, continue checking without updating
        # cur_alignment
        # if both true, update cur_alignments
        
        # case 1
        # original: 700 for search range, 100 for length threshold
        if t_pos - cur_t[1] <= 700:
            # case 2
            if q_pos - cur_q[1] <= 700:
                # update end of cur_q, cur_t
                cur_q[1] = q_pos
                cur_t[1] = t_pos
            else:
                continue
        else:
            if cur_t[1] - cur_t[0] >= 100 and cur_q[1] - cur_q[0] >= 100:
                alignments.append([cur_q, cur_t])
            cur_q, cur_t = [-1,-1], [-1,-1] # cur_alignments
    if cur_q[0] != -1 and cur_t[0] != -1:
        if cur_t[1] - cur_t[0] >= 100 and cur_q[1] - cur_q[0] >= 100:
            alignments.append([cur_q, cur_t])
    return alignments

def main():
    parser = argparse.ArgumentParser(description="Process DSB raw output.")
    parser.add_argument('-i', '--input', required=True, type=str, dest="infile",
            metavar="[FILE]", help="The raw output file from DSB.")

    infile = parser.parse_args().infile
    outfile = ('.'.join(infile.split('.')[:-1]) + '_mapped.' + 
                infile.split('.')[-1])
    
    # map result from dsb to ground truth
    try:
        d = open("{}".format(infile), 'r')
        output = open("{}_mapped.out".format(outfile), 'w')
        lines = d.read().split('\n')[:-1]
        runtime = lines[-1]
        lines = lines[:-1]
    
        for line in lines:
            if (len(line) == 0):
                continue
            query = int(line.split(':')[0].split(',')[0]) + 1
            target = int(line.split(':')[0].split(',')[1]) + 1
            print("on query {}\ttarget {}".format(query, target), end='\r')

            mapped = [eval(x) for x in line.split(':')[1].split(';')[:-1]]
            mapped = sorted(mapped, key=lambda s: s[1])
            
            # compile the results to get actual alignment on genome
            alignments = analyze_mapped(mapped)
            
            # write to output
            for item in alignments:
                # [query id] [direction] [q start] [q end] [t start] [t end]
                output.write("{} {} {} {} {} {}\n".format(query, target,
                    item[0][0], item[0][1], item[1][0], item[1][1]))
        output.write(runtime + '\n')
        output.close()
    except FileNotFoundError:
        print("File < {} > does not exist.".format(infile))
        exit(1)
        
if __name__ == "__main__":
    main()
