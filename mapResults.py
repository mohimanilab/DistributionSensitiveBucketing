import sys
import numpy as np
from collections import defaultdict

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
    return alignments

def main():
    assert(len(sys.argv) > 1)
    n = sys.argv[1]

    # a map to store all ground truth, for each query sequence (can have
    # multiple alignments)
    #gts = defaultdict(list)
    #o = open("ground_truth_alignments_{}.out".format(n), 'r')
    #for line in o.read().split('\n'):
    #    if (len(line) == 0):
    #        continue
    #    # split line to critical parts
    #    parts = line.split(' ')

    #    direction = 0
    #    query = int(parts[0].split('/')[0])
    #    t_start, t_end = int(parts[6]), int(parts[7])
    #    q_start, q_end = int(parts[9]), int(parts[10])
    #    if parts[2] != parts[3]:
    #        direction = 1
    #    
    #    # append new groundtruth class to list
    #    gts[query].append(GroundTruth(query, q_start, q_end, t_start, t_end, direction))
    
    # map result from dsb to ground truth
    d = open("dsb_{}.out".format(n), 'r')
    output = open("dsb_{}_mapped.out".format(n), 'w')
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
        
if __name__ == "__main__":
    main()
