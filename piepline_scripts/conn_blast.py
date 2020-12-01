#!/usr/bin/env python

import fileinput
import os, sys
import collections
import subprocess

prev = []
result = []

if len(sys.argv) < 2: 
    print "Usage: conn_blast.py blast_fmt6 blast_fmt6.conn"
    print "Required format1: Blastn outfmt6"
    print 'Required format2: Run blasn with -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids" '
    sys.exit()

input_file = open(sys.argv[1])
output_file = sys.argv[2]

oot = open(output_file, 'w')

## Required functions
def convert(str_a):
    try:
        str_a = float(str_a)
    except ValueError:
        str_a = str_a
    return(str_a)

def reverse_se(x):
    a = x[6]
    b = x[7]
    if  int(a) > int(b):
        a,b = b,a
        x[6] = a
        x[7] = b
    return(x)

def out_fmt(alist, fmt):
    prev_qname, prev_sname, identity, qcov_all, q_range, s_range, score, prev_tax_id = alist
    identity = round(identity, 2)
    qcov_all = int(qcov_all)
    score = int(score)
    q_final_start = min(q_range)
    q_final_end = max(q_range)
    s_final_start = min(s_range)
    s_final_end = max(s_range)

    if prev_tax_id != "NA":
        try:
            prev_tax_id = int(prev_tax_id)
        except ValueError:
            prev_tax_id = str(prev_tax_id)
        out = [prev_qname, prev_sname, identity, qcov_all, q_final_start, q_final_end, s_final_start, s_final_end, score, prev_tax_id]
    else:
        out = [prev_qname, prev_sname, identity, qcov_all, q_final_start, q_final_end, s_final_start, s_final_end, score]
    return out


## load table
line_list = []
for line in input_file:
    line = line.strip().split("\t")
    line = map(lambda x:convert(x), line)
    line_list.append(line)

# print line_list
## reset the start and sort lines (reverse the oeder when sstart < send)
line_list_reset = map(lambda x:reverse_se(x), line_list)
sort_lines = sorted(line_list_reset, key = lambda x:[x[0], x[1], x[6], x[7], -x[11]])

## Determine output format
if len(line_list_reset[1]) == 12:
    fmt = "fmt6"
elif len(line_list_reset[1]) == 13:
    fmt = "tmt6plustaxid"
else:
    print "Not required blast output format, exit." 
    sys.exit()

## ...
for s in sort_lines:
    blast_id = s[2]
    bitscore = int(s[11])
    qname = s[0] 
    sname = s[1]
    q_start = int(s[6])
    q_end = int(s[7])
    s_start = min([int(i) for i in s[8:10]])
    s_end = max([int(i) for i in s[8:10]])
    qcovlen = (q_end - q_start) + 1
    scovlen = (s_end - s_start) + 1
    raw_cov = s[3]
    ## For different blast output format
    if fmt == "tmt6plustaxid":
        tax_id = s[12]
    else:
        tax_id = "NA"

    if len(prev) != 0:
        if qname == prev_qname and sname == prev_sname: #spanned event
            if q_start == prev_q_end:
                identity = (identity*qcov_all + qcovlen*blast_id)/(qcov_all+qcovlen)
                # score += (score*qcov_all + bitscore*qcovlen)/(qcov_all+qcovlen)
                score += bitscore
                qcov_all += qcovlen - 1
            elif q_start > prev_q_end:
                identity = (identity*qcov_all + qcovlen*blast_id)/(qcov_all+qcovlen)
                score += bitscore
                # score += (score*qcov_all + bitscore*qcovlen)/(qcov_all+qcovlen)
                qcov_all += qcovlen
            elif q_start <= prev_q_end:
                if q_end > prev_q_end:
                    identity = (identity*qcov_all + qcovlen*blast_id)/(qcov_all+qcovlen)
                    extension_pct = (float(q_end) - float(prev_q_end))/((float(q_end) - float(q_start)))
                    # print extension_pct
                    score += bitscore*extension_pct
                    # score += (score*qcov_all + bitscore*qcovlen)/(qcov_all+qcovlen)
                    qcov_all += q_end - prev_q_end
                else:
                    continue
            # identity += blast_id
            # score += bitscore
            counter += 1
            last_line = "apart"
            s_range += [s_start, s_end]
            q_range += [q_start, q_end]
        else:
            if counter == 1:
                # qcov_all = prev_q_end -  prev_q_start + 1 
                qcov_all = p_raw_cov
            # mean_id = float(identity)/float(counter)
            # print s_range
            out_in = [prev_qname, prev_sname, identity, qcov_all, q_range, s_range, score, prev_tax_id]
            out = out_fmt(out_in, fmt)
            out = map(lambda x:str(x), out)
            out_str = "\t".join(out)
            oot.write(out_str + "\n")
            identity = blast_id ## reset value
            score = bitscore
            counter = 1
            qcov_all = qcovlen
            s_range = [s_start, s_end]
            q_range = [q_start, q_end]

    else:
        s_range = [s_start, s_end]
        q_range = [q_start, q_end]
        identity = blast_id
        score = bitscore
        counter = 1
        qcov_all = qcovlen
    # print s, qcov_all, qcovlen
    #assign prevs
    prev = s
    p_raw_cov = prev[3]
    prev_s_start = s_start
    prev_s_end = s_end
    prev_q_start = int(q_start)
    prev_q_end = int(q_end)
    prev_id = s[2]
    prev_qname = s[0]
    prev_sname = s[1]
    prev_score = int(s[11])
    if fmt == "tmt6plustaxid":
        prev_tax_id = prev[12]
    else:
        prev_tax_id = "NA"

# For last line
if last_line != "apart":
    qcov_all = p_raw_cov

out_in = [prev_qname, prev_sname, identity, qcov_all, q_range, s_range, score, prev_tax_id]
out = out_fmt(out_in, fmt)
## write to file
out = map(lambda x:str(x), out)
out_str = "\t".join(out)
oot.write(out_str + "\n")

oot.close()
## sort outfiles
sort_cmd = "sort -k1,1  -k9,9nr " + output_file + " -o " + output_file
# print sort_cmd
subprocess.call(sort_cmd, shell=True)