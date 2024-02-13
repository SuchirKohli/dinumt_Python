#!/usr/bin/env python

import math
from datetime import datetime
import re
import random
import time
import subprocess
import os
import sys
import getopt

version = "0.0.23"

opts = {
    'len_cluster_include': 600,
    'len_cluster_link': 800,
    'min_reads_cluster': 1,
    'min_clipped_seq': 5,
    'clipped_flank': 50,
    'max_num_clipped': 5,
    'include_mask': 0,
    'min_evidence': 3,
    'min_map_qual': 10,
    'filter_qual': 13,
    'filter_depth': 5,
    'max_read_cov': 200,
    'use_priors': 0,
    'mask_filename': "/home2/remills/projects/numts/numtS.bed",
    'info_filename': "/scratch/remills_flux/remills/sampleInfo.txt",
    'reference': "/nfs/remills-scratch/reference/hs37d5/hs37d5.fa",
    'mt_filename': "/nfs/remills-scratch/reference/GRCh37/Sequence/Chromosomes/MT.fa",
    'samtools': "/home2/remills/bin/samtools",
    'exonerate': "/home2/remills/apps/exonerate-2.2.0/bin/exonerate",
    'dir_tmp': "/tmp",
    'prefix': "numt",
    'len_mt': 16596,
    'ploidy': 2
}

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:m:n:", ["input_filename=", "output_filename=", "mask_filename=", "info_filename=",
                                                         "mt_filename=", "dir_tmp=", "chr=", "include_mask", "len_cluster_include=",
                                                         "len_cluster_link=", "min_reads_cluster=", "min_evidence=", "min_clipped_seq=",
                                                         "max_num_clipped=", "min_map_qual=", "max_read_cov=", "read_groups",
                                                         "breakpoint", "reference=", "samtools=", "exonerate=", "use_priors",
                                                         "by_chr_dir", "prefix=", "ucsc", "help", "verbose"])
except getopt.GetoptError:
    print("Invalid option")
    sys.exit(2)

input_filename = output_filename = mask_filename = info_filename = mt_filename = dir_tmp = chr = None
include_mask = read_groups = breakpoint = use_priors = by_chr_dir = ucsc = help = verbose = False

for opt, arg in opts:
    if opt in ("-i", "--input_filename"):
        input_filename = arg
    elif opt in ("-o", "--output_filename"):
        output_filename = arg
    elif opt in ("-m", "--mask_filename"):
        mask_filename = arg
    elif opt in ("-n", "--info_filename"):
        info_filename = arg
    elif opt == "--mt_filename":
        mt_filename = arg
    elif opt == "--dir_tmp":
        dir_tmp = arg
    elif opt == "--chr":
        chr = arg
    elif opt == "--include_mask":
        include_mask = True
    elif opt == "--len_cluster_include":
        opts['len_cluster_include'] = int(arg)
    elif opt == "--len_cluster_link":
        opts['len_cluster_link'] = int(arg)
    elif opt == "--min_reads_cluster":
        opts['min_reads_cluster'] = int(arg)
    elif opt == "--min_evidence":
        opts['min_evidence'] = int(arg)
    elif opt == "--min_clipped_seq":
        opts['min_clipped_seq'] = int(arg)
    elif opt == "--max_num_clipped":
        opts['max_num_clipped'] = int(arg)
    elif opt == "--min_map_qual":
        opts['min_map_qual'] = int(arg)
    elif opt == "--max_read_cov":
        opts['max_read_cov'] = int(arg)
    elif opt == "--read_groups":
        read_groups = True
    elif opt == "--breakpoint":
        breakpoint = True
    elif opt == "--reference":
        opts['reference'] = arg
    elif opt == "--samtools":
        opts['samtools'] = arg
    elif opt == "--exonerate":
        opts['exonerate'] = arg
    elif opt == "--use_priors":
        use_priors = True
    elif opt == "--by_chr_dir":
        by_chr_dir = True
    elif opt == "--prefix":
        opts['prefix'] = arg
    elif opt == "--ucsc":
        ucsc = True
    elif opt == "--help":
        help = True
    elif opt == "--verbose":
        verbose = True

seq_num = 0
seq_hash = {}

sorted_hash = {}

i = 1
sample_hash = {}
infile_hash = {}
group_hash = {}
outfile_hash = {}
mask_hash = {}
data_hash = {}


###############################################

def mappedClippedSeq(clipSeq):
    if 'N' in clipSeq:
        return 0

    random_num = str(random.randint(0, 100000))
    time_data = time.localtime(time.time())
    timestamp = ''.join(map(str, time_data))
    file_clp = f"{opts['dir_tmp']}/clp{timestamp}{random_num}.fa"

    with open(file_clp, 'w') as clp_file:
        clp_file.write(f">clp\n{clipSeq}\n")

    results = subprocess.getoutput(
        f"{opts['exonerate']} --model affine:local --exhaustive yes --percent 70 --showvulgar no --showalignment no --showcigar yes {file_clp} {opts['mt_filename']} 2> /dev/null"
    )
    match = re.search(r'cigar:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|-)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|-)\s+(\d+)\s+(.*?)\n', results)

    os.remove(file_clp)

    if match:
        return 1
    else:
        return 0

def assessBreaks(infile_hash, sample_hash, data_hash):
    print("Entering assessBreaks()")
    for var in data_hash:
        print(f"-variant {var}")
        clipped_pos = {}
        sum_clipped = 0
        win_start = 1e10
        win_end = 0
        max_first = 0
        max_second = 0

        for sample in data_hash[var]:
            for c_pos in data_hash[var][sample]['clipPos']:
                if c_pos == -1:
                    continue
                if data_hash[var][sample]['clipPos'][c_pos][1] is not None:
                    if mappedClippedSeq(data_hash[var][sample]['clipSeq'][c_pos]):
                        clipped_pos[c_pos] = clipped_pos.get(c_pos, 0) + data_hash[var][sample]['clipPos'][c_pos][1]
                        sum_clipped += data_hash[var][sample]['clipPos'][c_pos][1]

        win_len = len(clipped_pos)
        sorted_pos = sorted(clipped_pos, key=lambda x: clipped_pos[x], reverse=True)

        if win_len == 0:
            continue
        elif win_len == 1:
            infile_hash[var]['leftBkpt'] = sorted_pos[0]
            infile_hash[var]['rightBkpt'] = sorted_pos[0] + 1
        elif win_len == 2:
            infile_hash[var]['leftBkpt'] = sorted_pos[0]
            infile_hash[var]['rightBkpt'] = sorted_pos[0] + 1
            if sorted_pos[1] > sorted_pos[0]:
                infile_hash[var]['rightBkpt'] = sorted_pos[1]
            else:
                infile_hash[var]['leftBkpt'] = sorted_pos[1]
                infile_hash[var]['rightBkpt'] -= 1
        else:
            mean_clipped = sum_clipped / win_len
            sum_squares = sum((clipped_pos[p] - mean_clipped) ** 2 for p in clipped_pos)
            print(f"\tmean clipped: {mean_clipped}")
            sd_clipped = (sum_squares / (win_len - 1)) ** 0.5
            print(f"\tsd clipped: {sd_clipped}")
            print(f"\t1st clipped: {clipped_pos.get(sorted_pos[0], 0)}")
            print(f"\t2nd clipped: {clipped_pos.get(sorted_pos[1], 0)}")

            if clipped_pos.get(sorted_pos[0], 0) > mean_clipped + 1 * sd_clipped:
                print("\t*updating to 1st")
                infile_hash[var]['leftBkpt'] = sorted_pos[0]
                infile_hash[var]['rightBkpt'] = sorted_pos[0] + 1

            if clipped_pos.get(sorted_pos[1], 0) > mean_clipped + 1 * sd_clipped:
                print("\t*updating to 1st and 2nd")
                if sorted_pos[1] > sorted_pos[0]:
                    infile_hash[var]['rightBkpt'] = sorted_pos[1]
                else:
                    infile_hash[var]['leftBkpt'] = sorted_pos[1]
                    infile_hash[var]['rightBkpt'] -= 1

        print(f"-variant {var} completed")
    print("Exiting assessBreaks()\n")

def scoreData(data_hash, sample_hash):
    print("Entering scoreData()")
    for var in data_hash:
        print(f"Variant: {var}")
        priors = {}
        geno_freq = {}
        
        for geno in range(opts['ploidy'] + 1):
            priors[geno] = 1 / (opts['ploidy'] + 1)
            geno_freq[geno] = {'old': 0, 'new': 0}

        num_iteration = 0
        sum_geno_freq = 0

        while num_iteration < 10:
            num_iteration += 1

            for sample in sample_hash:
                print(f"\tSample: {sample}")

                # Zero out alternative supporting evidence if below threshold
                if (data_hash[var][sample]['numAltRP'] + data_hash[var][sample]['numAltSR'] <
                        opts['min_evidence']):
                    data_hash[var][sample]['numAltRP'] = 0
                    data_hash[var][sample]['numAltSR'] = 0
                    data_hash[var][sample]['qualAltRP'] = []
                    data_hash[var][sample]['qualAltSR'] = []

                num_ref_rp = data_hash[var][sample]['numRefRP']
                num_alt_rp = data_hash[var][sample]['numAltRP']
                num_ref_sr = data_hash[var][sample]['numRefSR']
                num_alt_sr = data_hash[var][sample]['numAltSR']

                for geno in range(opts['ploidy'] + 1):
                    data_hash[var][sample]['pl'][geno] = 0
                    data_hash[var][sample]['gl'][geno] = 0
                    data_hash[var][sample]['gl0'][geno] = 0

                data_hash[var][sample]['gq'] = 0
                data_hash[var][sample]['gt'] = "./."
                data_hash[var][sample]['ft'] = "LowQual"

                if num_alt_rp + num_ref_rp + num_alt_sr + num_ref_sr == 0:
                    continue

                if data_hash[var][sample]['avgQ'] <= 0 or data_hash[var][sample]['avgQ'] > 1:
                    data_hash[var][sample]['avgQ'] = 0.999999

                data_hash[var][sample]['avgQ'] = 0.00001

                for g in range(opts['ploidy'] + 1):
                    geno = opts['ploidy'] - g
                    if num_alt_rp + num_ref_rp > 0 and 1 / opts['ploidy'] ** (num_alt_rp + num_ref_rp) > 0:
                        data_hash[var][sample]['gl0'][geno] += calcGl(
                            opts['ploidy'], g, num_alt_rp + num_ref_rp, num_ref_rp,
                            data_hash[var][sample]['qualRefRP'], data_hash[var][sample]['qualAltRP']
                        )
                    if num_alt_sr + num_ref_sr > 0 and opts['breakpoint'] and \
                            1 / opts['ploidy'] ** (num_alt_sr + num_ref_sr) > 0:
                        data_hash[var][sample]['gl0'][geno] += calcGl(
                            opts['ploidy'], g, num_alt_sr + num_ref_sr, num_ref_sr,
                            data_hash[var][sample]['qualRefSR'], data_hash[var][sample]['qualAltSR']
                        )
                    print(f"\tgl0 returned from calcGL() for geno {geno}: "
                          f"{data_hash[var][sample]['gl0'][geno]}")

                    if data_hash[var][sample]['gl0'][geno] < -255:
                        data_hash[var][sample]['gl0'][geno] = -255

                    data_hash[var][sample]['gl'][geno] = math.log10(priors[geno]) + \
                                                          data_hash[var][sample]['gl0'][geno]

                sorted_geno = sorted(data_hash[var][sample]['gl'], key=lambda x: data_hash[var][sample]['gl'][x], reverse=True)

                # Calculate PL from GL
                for geno in range(opts['ploidy'] + 1):
                    print(f"\tCalculating pl from geno ({data_hash[var][sample]['gl'][geno]})")
                    data_hash[var][sample]['pl'][geno] = int(-10 * data_hash[var][sample]['gl'][geno])
                    if data_hash[var][sample]['pl'][geno] > 255:
                        data_hash[var][sample]['pl'][geno] = 255
                    print(f"\t...{data_hash[var][sample]['pl'][geno]}")

                # Normalize PL to most likely genotype
                for geno in range(opts['ploidy'] + 1):
                    data_hash[var][sample]['pl'][geno] -= data_hash[var][sample]['pl'][sorted_geno[0]]

                # Determine genotype quality
                data_hash[var][sample]['gq'] = int(10 * (data_hash[var][sample]['gl'][sorted_geno[0]] -
                                                          data_hash[var][sample]['gl'][sorted_geno[1]]))
                print(f"\t...{data_hash[var][sample]['gq']}")

                gt = "0/0"
                if sorted_geno[0] == 1:
                    gt = "0/1"
                elif sorted_geno[0] == 2:
                    gt = "1/1"
                data_hash[var][sample]['gt'] = gt

                if num_alt_rp + num_ref_rp + num_alt_sr + num_ref_sr < opts['filter_depth'] or \
                        data_hash[var][sample]['gq'] < opts['filter_qual']:
                    data_hash[var][sample]['ft'] = "LowQual"
                else:
                    data_hash[var][sample]['ft'] = "PASS"

                geno_freq[sorted_geno[0]]['new'] += 1
                sum_geno_freq += 1

            is_converged = 1

            if sum_geno_freq == 0:
                print("No genotypes passing filters")
                break

            for geno in range(opts['ploidy'] + 1):
                if geno_freq[geno]['new'] == 0:
                    priors[geno] = 1 / sum_geno_freq
                else:
                    priors[geno] = geno_freq[geno]['new'] / sum_geno_freq

                if geno_freq[geno]['old'] != geno_freq[geno]['new']:
                    is_converged = 0

                geno_freq[geno]['old'] = geno_freq[geno]['new']
                geno_freq[geno]['new'] = 0

            if is_converged or not opts['use_priors']:
                break

    print("Exiting scoreData()\n\n")

def get_date():
    now = datetime.now()
    year = now.year
    month = now.month
    day = now.day

    fmonth = f"{month:02d}"
    fday = f"{day:02d}"

    return f"{year}{fmonth}{fday}"

def report(infile_hash, data_hash):
    print("Entering report()")
    
    # Open output file
    if opts.get("output_filename"):
        with open(opts["output_filename"], "w") as foutname1:
            filedate = get_date()
            foutname1.write(f"##fileformat=VCFv4.1\n")
            # ... (other header lines)
            foutname1.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            for sample in sorted(data_hash[list(data_hash.keys())[0]]):
                foutname1.write(f"\t{sample}")
            foutname1.write("\n")

            vars_sorted = sorted(infile_hash.keys(), key=lambda v: (infile_hash[v]["chr"], infile_hash[v]["start"]))
            samples_sorted = sorted(data_hash[vars_sorted[0]])

            index = 1
            for var in vars_sorted:
                chrom = infile_hash[var]["chr"]
                id_ = infile_hash[var]["id"]
                alt = "<INS:MT>"
                qual = infile_hash[var]["qual"]
                filter_ = infile_hash[var]["filter"]
                info = {"IMPRECISE": None, "CIPOS": "0,0", "CIEND": "0,0", "END": infile_hash[var]["end"], "SVTYPE": "INS"}
                refline = subprocess.check_output([opts["samtools"], "faidx", opts["reference"], f"{chrom}:{infile_hash[var]['start']}-{infile_hash[var]['start']}"], text=True).split("\n")[1]
                ref = refline if refline else "N"
                format_ = "GT:FT:GL0:GQ:PL"

                info_str = ";".join([f"{key}={value}" if value is not None else key for key, value in info.items()])

                foutname1.write(f"{chrom}\t{infile_hash[var]['start']}\t{id_}\t{ref}\t{alt}\t{qual}\t{filter_}\t{info_str}\t{format_}")
                for sample in samples_sorted:
                    gls = [f"{float(data_hash[var][sample]['gl0'][geno]):.2f}" for geno in range(opts["ploidy"] + 1)]
                    pls = [str(data_hash[var][sample]['pl'][geno]) for geno in range(opts["ploidy"] + 1)]
                    gl = ",".join(gls)
                    pl = ",".join(pls)
                    foutname1.write(f"\t{data_hash[var][sample]['gt']}:{data_hash[var][sample]['ft']}:{gl}:{data_hash[var][sample]['gq']}:{pl}")
                foutname1.write("\n")
                index += 1
    print("Exiting report()\n")

def calcGl(m, g, k, l, er, ea):
    print("in calcGl():") if opts["verbose"] else None
    print(f"\t{m}\t{g}\t{k}\t{l}\t{er}\t{ea}") if opts["verbose"] else None

    if 1 / m**k <= 0:
        raise ValueError(f"problem in calcGL 1, \t{m}\t{g}\t{k}\t{l}\t{er}\t{ea}")

    gl = math.log10(1 / (m**k))

    for e in er:
        if (m - g) * e + (1 - e) * g <= 0:
            raise ValueError(f"problem in calcGL 2, \t{m}\t{g}\t{k}\t{l}\t{e}")
        gl += math.log10((m - g) * e + (1 - e) * g)

    for e in ea:
        if (m - g) * (1 - e) + g * e <= 0:
            raise ValueError(f"problem in calcGL 3, \t{m}\t{g}\t{k}\t{l}\t{e}")
        gl += math.log10((m - g) * (1 - e) + g * e)

    return gl

def get_input(infile_hash, mask_hash, sample_hash):
    print("Entering getInput()") if opts["verbose"] else None

    # Input sample information
    with open(opts["info_filename"], "r") as info_file:
        header = []
        fields = {}
        for line in info_file:
            line = line.strip()
            if line.startswith("sample"):
                header = line.split("\t")
                continue
            row = line.split("\t")
            for i in range(len(header)):
                fields[header[i]] = row[i]

            sample = fields["sample"]
            sample_hash[sample]["pop"] = fields["pop"]
            sample_hash[sample]["filename"] = fields["filename"]
            if fields["median_insert_size"] != "NA":
                sample_hash[sample]["winlen"] = int(fields["median_insert_size"]) + 3 * int(fields["median_absolute_deviation"])
                sample_hash[sample]["median_insert_size"] = int(fields["median_insert_size"])
                sample_hash[sample]["median_absolute_deviation"] = int(fields["median_absolute_deviation"])
            else:
                sample_hash[sample]["winlen"] = 500
                sample_hash[sample]["median_insert_size"] = 400
                sample_hash[sample]["median_absolute_deviation"] = 40
            sample_hash[sample]["meancoverage"] = float(fields["mean_coverage"])
            if "read_groups" in fields:
                rg_list = fields["read_groups"].split(",")
                sample_hash[sample]["read_groups"] = {rg: 1 for rg in rg_list}

    # Input mask coordinates
    with open(opts["mask_filename"], "r") as mask_file:
        for line in mask_file:
            chr, start, end, _ = line.strip().split("\t")
            chr = chr.replace("chr", "")
            mask_hash[chr][int(start)] = int(end)

    # Input variant coordinates
    with open(opts["input_filename"], "r") as vars_file:
        varnum = 1

        # BED FORMAT
        # for line in vars_file:
        #     chr, start, end = line.strip().split("\t")
        #     chr = chr.replace("chr", "")
        #     infile_hash[varnum]["chr"] = chr
        #     infile_hash[varnum]["start"] = int(start)
        #     infile_hash[varnum]["end"] = int(end)
        #     varnum += 1

        # VCF FORMAT
        for line in vars_file:
            if line.startswith("#"):
                continue
            chr, pos, id, ref, alt, qual, filter, info = line.strip().split("\t")
            chr = chr.replace("chr", "")
            if "chr" in opts and chr != opts["chr"]:
                continue
            infile_hash[varnum]["chr"] = chr
            infile_hash[varnum]["id"] = id
            infile_hash[varnum]["start"] = int(pos) + 1
            end_match = re.search(r"END=(\d+)", info)
            mstart_match = re.search(r"MSTART=(\d+)", info)
            mend_match = re.search(r"MEND=(\d+)", info)
            mlen_match = re.search(r"MLEN=(\d+)", info)
            samples_match = re.search(r"SAMPLES=(\w+?);", info)

            if end_match:
                infile_hash[varnum]["end"] = int(end_match.group(1))
            if mlen_match:
                infile_hash[varnum]["mlen"] = int(mlen_match.group(1))
                infile_hash[varnum]["mstart"] = int(mstart_match.group(1))
                infile_hash[varnum]["mend"] = int(mend_match.group(1))
            if samples_match:
                infile_hash[varnum]["samples"] = samples_match.group(1)
            
            infile_hash[varnum]["filter"] = filter
            infile_hash[varnum]["qual"] = qual
            varnum += 1

    if opts["verbose"]:
        print("Exiting getInput()\n\n") 
    
def get_mate_info(qname, rnext, pnext, readgroup_hash, filename):
    command = ""
    if opts["by_chr_dir"]:
        if opts["ucsc"]:
            command = f"{opts['samtools']} view {filename}chr{rnext}.*bam chr{rnext}:{pnext}-{pnext} |"
        else:
            command = f"{opts['samtools']} view {filename}{rnext}.*bam {rnext}:{pnext}-{pnext} |"
    else:
        if opts["ucsc"]:
            command = f"{opts['samtools']} view {filename} chr{rnext}:{pnext}-{pnext} |"
        else:
            command = f"{opts['samtools']} view {filename} {rnext}:{pnext}-{pnext} |"
    
    with subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
        output, error = proc.communicate()

    cFlag = 0
    cPos = -1
    clipside = "n"
    clipsize = -1
    clipseq = ""
    seqLen = 0
    seq = ""
    matchLen = 0

    for line in output.splitlines():
        fields = line.strip().split("\t")
        m_qname, m_flag, m_rname, m_pos, m_mapq, m_cigar, m_rnext, m_pnext, m_tlen, m_seq, m_qual, opt = fields[:12]
        read_group_match = re.search(r'RG:Z:(\S+)', line)
        read_group = read_group_match.group(1) if read_group_match else None

        if m_qname != qname:
            continue

        if opts["read_groups"] and read_group is None:
            continue
        elif opts["read_groups"] and read_group not in readgroup_hash:
            continue

        cFlag, cPos, clipside, clipsize, clipseq = get_soft_clip_info(int(m_pos), m_cigar, m_qual, m_seq)
        seqLen = sum(map(int, re.findall(r'(\d+)M', m_cigar)))
        matchLen = sum(map(int, re.findall(r'(\d+)M', m_cigar)))
        seq = m_seq

        for length in re.findall(r'(\d+)N', m_cigar):
            seqLen += int(length)

        for length in re.findall(r'(\d+)D', m_cigar):
            seqLen += int(length)

    return cFlag, cPos, clipside, clipsize, seqLen, matchLen, seq

def refine_data(infile_hash, mask_hash, sample_hash, data_hash):
    print("Entering refineData()\n") if opts["verbose"] else None

    num_samp = 0
    for var in infile_hash:
        if infile_hash[var]["leftBkpt"] is None:
            continue
        for sample in sample_hash:
            num_samp += 1

            chr_var = infile_hash[var]["chr"]
            l_start = infile_hash[var]["leftBkpt"] - sample_hash[sample]["winlen"]
            l_end = infile_hash[var]["leftBkpt"]
            r_start = infile_hash[var]["rightBkpt"]
            r_end = infile_hash[var]["rightBkpt"] + sample_hash[sample]["winlen"]
            command = ""
            chr_m = ""

            print(f"REFINED: {var}\t{chr_var}\t{l_start}\t{l_end}\t{r_start}\t{r_end}\n") if opts["verbose"] else None
            if opts["by_chr_dir"]:
                if opts["ucsc"]:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']}chr{chr_var}.*bam chr{chr_var}:{l_start}-{r_end} |"
                    chr_m = "M"
                else:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']}{chr_var}.*bam {chr_var}:{l_start}-{r_end} |"
                    chr_m = "MT"
            else:
                if opts["ucsc"]:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']} chr{chr_var}:{l_start}-{r_end} |"
                    chr_m = "M"
                else:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']} {chr_var}:{l_start}-{r_end} |"
                    chr_m = "MT"

            data_hash[var][sample]["numRefRP"] = 0
            data_hash[var][sample]["numAltRP"] = 0
            data_hash[var][sample]["numRefSR"] = 0
            data_hash[var][sample]["numAltSR"] = 0
            data_hash[var][sample]["qualRefRP"] = []
            data_hash[var][sample]["qualAltRP"] = []
            data_hash[var][sample]["qualRefSR"] = []
            data_hash[var][sample]["qualAltSR"] = []
            data_hash[var][sample]["avgQ"] = 0

            found = {}

            with subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
                output, error = proc.communicate()

            for line in output.splitlines():
                fields = line.strip().split("\t")
                qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fields[:11]
                rname = rname.replace("chr", "")
                rnext = rnext.replace("chr", "")

                mapE = 10**(-1 * float(mapq) / 10)
                read_group_match = re.search(r'RG:Z:(\S+)', line)
                read_group = read_group_match.group(1) if read_group_match else None

                if opts["read_groups"] and read_group is None:
                    continue
                elif opts["read_groups"] and read_group not in sample_hash[sample]["read_groups"]:
                    continue

                if float(mapq) < opts["min_map_qual"]:
                    continue

                dir = 0 if int(flag) & 16 == 0 else 1
                dnext = 0 if int(flag) & 32 == 0 else 1

                seqLen = 0
                matchLen = 0
                for m in re.finditer(r'(\d+)M', cigar):
                    matchLen += int(m.group(1))
                    seqLen += int(m.group(1))
                for n in re.finditer(r'(\d+)N', cigar):
                    seqLen += int(n.group(1))
                for d in re.finditer(r'(\d+)D', cigar):
                    seqLen += int(d.group(1))

                cFlag, cPos, clipside, clipsize, clipseq = get_soft_clip_info(int(pos), cigar, qual, seq)
                m_cFlag, m_cPos, m_clipside, m_clipsize, m_seqLen, m_matchLen, m_seq = get_mate_info(qname, rname, pnext, sample_hash[sample]["read_groups"], sample_hash[sample]["filename"])

                if cFlag == 1 and (abs(cPos - infile_hash[var]["leftBkpt"]) <= opts["min_clipped_seq"] or abs(cPos - infile_hash[var]["rightBkpt"]) <= opts["min_clipped_seq"]):
                    data_hash[var][sample]["numAltSR"] += 1
                    data_hash[var][sample]["qualAltSR"].append(mapE)
                elif matchLen >= 0.95 * len(seq) and ((pos < infile_hash[var]["leftBkpt"] and pos + seqLen - 1 > infile_hash[var]["leftBkpt"] and abs(infile_hash[var]["leftBkpt"] - pos) > opts["min_clipped_seq"] and abs(infile_hash[var]["leftBkpt"] - (pos + seqLen - 1)) > opts["min_clipped_seq"]) or (pos < infile_hash[var]["rightBkpt"] and pos + seqLen - 1 > infile_hash[var]["rightBkpt"] and abs(infile_hash[var]["rightBkpt"] - pos) > opts["min_clipped_seq"] and abs(infile_hash[var]["rightBkpt"] - (pos + seqLen - 1)) > opts["min_clipped_seq"])):
                    data_hash[var][sample]["numRefSR"] += 1
                    data_hash[var][sample]["qualRefSR"].append(mapE)

                if rnext == chr_m or check_mask_overlap(chr_var, int(pos), mask_hash):
                    if (dir == 0 and pos >= l_start and pos <= l_end) or (dir == 1 and pos >= r_start and pos <= r_end):
                        data_hash[var][sample]["numAltRP"] += 1
                        data_hash[var][sample]["qualAltRP"].append(mapE)
                elif dir == 0 and dnext == 1 and rnext == "=" and qname not in found:
                    found[qname] = 1
                    if abs(int(tlen)) > sample_hash[sample]["median_insert_size"] + 3 * sample_hash[sample]["median_absolute_deviation"]:
                        continue
                    if abs(int(tlen)) < sample_hash[sample]["median_insert_size"] - 3 * sample_hash[sample]["median_absolute_deviation"]:
                        continue

                    if cPos > -1 and (abs(cPos - l_end) <= opts["min_clipped_seq"] and cFlag == 1):
                        data_hash[var][sample]["numAltRP"] += 1
                        data_hash[var][sample]["qualAltRP"].append(mapE)
                    elif m_cPos > -1 and (abs(m_cPos - l_end) <= opts["min_clipped_seq"] and m_cFlag == 1):
                        data_hash[var][sample]["numAltRP"] += 1
                        data_hash[var][sample]["qualAltRP"].append(mapE)
                    elif cPos > -1 and (abs(cPos - r_start) <= opts["min_clipped_seq"] and cFlag == 1):
                        data_hash[var][sample]["numAltRP"] += 1
                        data_hash[var][sample]["qualAltRP"].append(mapE)
                    elif m_cPos > -1 and (abs(m_cPos - r_start) <= opts["min_clipped_seq"] and m_cFlag == 1):
                        data_hash[var][sample]["numAltRP"] += 1
                        data_hash[var][sample]["qualAltRP"].append(mapE)
                    elif (pos <= l_end and pnext + m_seqLen >= l_end) or (pos <= r_start and pnext + m_seqLen >= r_start):
                        if (abs(pos - l_end) > opts["min_clipped_seq"] and abs(pos + seqLen - 1 - l_end) > opts["min_clipped_seq"] and abs(pnext - l_end) > opts["min_clipped_seq"] and abs(pnext + m_seqLen - 1 - l_end) > opts["min_clipped_seq"] and abs(pos - r_start) > opts["min_clipped_seq"] and abs(pos + seqLen - 1 - r_start) > opts["min_clipped_seq"] and abs(pnext - r_start) > opts["min_clipped_seq"] and abs(pnext + m_seqLen - 1 - r_start) > opts["min_clipped_seq"] and matchLen >= 0.95 * len(seq) and m_matchLen >= 0.95 * len(m_seq)):
                            data_hash[var][sample]["numRefRP"] += 1
                            data_hash[var][sample]["qualRefRP"].append(mapE)

                data_hash[var][sample]["avgQ"] += mapE

            if data_hash[var][sample]["numRefRP"] + data_hash[var][sample]["numAltRP"] > 0:
                data_hash[var][sample]["avgQ"] /= (data_hash[var][sample]["numRefRP"] + data_hash[var][sample]["numAltRP"])
            else:
                data_hash[var][sample]["avgQ"] = 0.00001

            print(f"{sample}\t{var}\t{data_hash[var][sample]['numRefRP']}\t{data_hash[var][sample]['numRefSR']}\t{data_hash[var][sample]['numAltRP']}\t{data_hash[var][sample]['numAltSR']}\t{data_hash[var][sample]['avgQ']}\n") if opts["verbose"] else None

    print("Exiting refineData()\n\n") if opts["verbose"] else None

def get_data(infile_hash, mask_hash, sample_hash, data_hash):
    print("Entering getData()\n") if opts["verbose"] else None

    num_samp = 0
    for var in infile_hash:
        for sample in sample_hash:
            num_samp += 1
            chr_var = infile_hash[var]["chr"]
            l_start = infile_hash[var]["start"] - sample_hash[sample]["winlen"]
            l_end = infile_hash[var]["start"]
            r_start = infile_hash[var]["end"]
            r_end = infile_hash[var]["end"] + sample_hash[sample]["winlen"]
            command = ""
            chr_m = ""

            if opts["by_chr_dir"]:
                if opts["ucsc"]:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']}chr{chr_var}.*bam chr{chr_var}:{l_start}-{r_end} |"
                    chr_m = "M"
                else:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']}{chr_var}.*bam {chr_var}:{l_start}-{r_end} |"
                    chr_m = "MT"
            else:
                if opts["ucsc"]:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']} chr{chr_var}:{l_start}-{r_end} |"
                    chr_m = "M"
                else:
                    command = f"{opts['samtools']} view {sample_hash[sample]['filename']} {chr_var}:{l_start}-{r_end} |"
                    chr_m = "MT"

            data_hash[var][sample]["numRefRP"] = 0
            data_hash[var][sample]["numAltRP"] = 0
            data_hash[var][sample]["numRefSR"] = 0
            data_hash[var][sample]["numAltSR"] = 0
            data_hash[var][sample]["qualRefRP"] = []
            data_hash[var][sample]["qualAltRP"] = []
            data_hash[var][sample]["qualRefSR"] = []
            data_hash[var][sample]["qualAltSR"] = []
            data_hash[var][sample]["avgQ"] = 0.000001

            with open(command, "r") as sam:
                for line in sam:
                    fields = line.strip().split("\t")
                    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fields[:11]
                    rname = rname.replace("chr", "")
                    rnext = rnext.replace("chr", "")

                    mapE = 10**(-1 * float(mapq) / 10)
                    read_group_match = re.search(r'RG:Z:(\S+)', line)
                    read_group = read_group_match.group(1) if read_group_match else None

                    if opts["read_groups"] and read_group is None:
                        continue
                    elif opts["read_groups"] and read_group not in sample_hash[sample]["read_groups"]:
                        continue

                    if float(mapq) < opts["min_map_qual"]:
                        continue

                    dir = 0 if int(flag) & 16 == 0 else 1
                    dnext = 0 if int(flag) & 32 == 0 else 1

                    seqLen = 0
                    matchLen = 0
                    for m in re.finditer(r'(\d+)M', cigar):
                        matchLen += int(m.group(1))
                        seqLen += int(m.group(1))
                    for n in re.finditer(r'(\d+)N', cigar):
                        seqLen += int(n.group(1))
                    for d in re.finditer(r'(\d+)D', cigar):
                        seqLen += int(d.group(1))

                    cFlag, cPos, clipside, clipsize, clipseq = get_soft_clip_info(int(pos), cigar, qual, seq)
                    m_cFlag, m_cPos, m_clipside, m_clipsize, m_seqLen, m_matchLen, m_seq = get_mate_info(qname, rname, pnext, sample_hash[sample]["read_groups"], sample_hash[sample]["filename"])

                    if cFlag == 1 and (abs(cPos - infile_hash[var]["start"]) <= opts["min_clipped_seq"] or abs(cPos - infile_hash[var]["end"]) <= opts["min_clipped_seq"]):
                        data_hash[var][sample]["numAltSR"] += 1
                        data_hash[var][sample]["qualAltSR"].append(mapE)
                    elif matchLen >= 0.95 * len(seq) and ((pos < infile_hash[var]["start"] and pos + seqLen - 1 > infile_hash[var]["start"] and abs(infile_hash[var]["start"] - pos) > opts["min_clipped_seq"] and abs(infile_hash[var]["start"] - (pos + seqLen - 1)) > opts["min_clipped_seq"]) or (pos < infile_hash[var]["end"] and pos + seqLen - 1 > infile_hash[var]["end"] and abs(infile_hash[var]["end"] - pos) > opts["min_clipped_seq"] and abs(infile_hash[var]["end"] - (pos + seqLen - 1)) > opts["min_clipped_seq"])):
                        data_hash[var][sample]["numRefSR"] += 1
                        data_hash[var][sample]["qualRefSR"].append(mapE)

                    if rnext == chr_m or check_mask_overlap(chr_var, int(pos), mask_hash):
                        if (dir == 0 and pos >= l_start and pos <= l_end) or (dir == 1 and pos >= r_start and pos <= r_end):
                            data_hash[var][sample]["numAltRP"] += 1
                            data_hash[var][sample]["qualAltRP"].append(mapE)
                    elif dir == 0 and dnext == 1 and rnext == "=" and qname not in found:
                        found[qname] = 1
                        if abs(int(tlen)) > sample_hash[sample]["median_insert_size"] + 3 * sample_hash[sample]["median_absolute_deviation"]:
                            continue
                        if abs(int(tlen)) < sample_hash[sample]["median_insert_size"] - 3 * sample_hash[sample]["median_absolute_deviation"]:
                            continue

                        if cPos > -1 and (abs(cPos - l_end) <= opts["min_clipped_seq"] and cFlag == 1):
                            data_hash[var][sample]["numAltRP"] += 1
                            data_hash[var][sample]["qualAltRP"].append(mapE)
                        elif m_cPos > -1 and (abs(m_cPos - l_end) <= opts["min_clipped_seq"] and m_cFlag == 1):
                            data_hash[var][sample]["numAltRP"] += 1
                            data_hash[var][sample]["qualAltRP"].append(mapE)
                        elif cPos > -1 and (abs(cPos - r_start) <= opts["min_clipped_seq"] and cFlag == 1):
                            data_hash[var][sample]["numAltRP"] += 1
                            data_hash[var][sample]["qualAltRP"].append(mapE)
                        elif m_cPos > -1 and (abs(m_cPos - r_start) <= opts["min_clipped_seq"] and m_cFlag == 1):
                            data_hash[var][sample]["numAltRP"] += 1
                            data_hash[var][sample]["qualAltRP"].append(mapE)
                        elif (pos <= l_end and pnext + m_seqLen >= l_end) or (pos <= r_start and pnext + m_seqLen >= r_start):
                            if (abs(pos - l_end) > opts["min_clipped_seq"] and abs(pos + seqLen - 1 - l_end) > opts["min_clipped_seq"] and abs(pnext - l_end) > opts["min_clipped_seq"] and abs(pnext + m_seqLen - 1 - l_end) > opts["min_clipped_seq"] and abs(pos - r_start) > opts["min_clipped_seq"] and abs(pos + seqLen - 1 - r_start) > opts["min_clipped_seq"] and abs(pnext - r_start) > opts["min_clipped_seq"] and abs(pnext + m_seqLen - 1 - r_start) > opts["min_clipped_seq"] and matchLen >= 0.95 * len(seq) and m_matchLen >= 0.95 * len(m_seq)):
                                data_hash[var][sample]["numRefRP"] += 1
                                data_hash[var][sample]["qualRefRP"].append(mapE)

                    data_hash[var][sample]["avgQ"] += mapE

                if data_hash[var][sample]["numRefRP"] + data_hash[var][sample]["numAltRP"] > 0:
                    data_hash[var][sample]["avgQ"] /= (data_hash[var][sample]["numRefRP"] + data_hash[var][sample]["numAltRP"])
                else:
                    data_hash[var][sample]["avgQ"] = 0.00001

                print(f"{sample}\t{var}\t{data_hash[var][sample]['numRefRP']}\t{data_hash[var][sample]['numRefSR']}\t{data_hash[var][sample]['numAltRP']}\t{data_hash[var][sample]['numAltSR']}\t{data_hash[var][sample]['avgQ']}\n") if opts["verbose"] else None

    print("Exiting getData()\n\n") if opts["verbose"] else None

def check_mask_overlap(chr_var, pos, mask_hash):
    is_mask_overlap = 0
    if chr_var in mask_hash:
        for mask_start, mask_end in mask_hash[chr_var].items():
            if pos >= mask_start and pos <= mask_end:
                is_mask_overlap = 1
                break
    return is_mask_overlap

def get_soft_clip_info(pos, cigar, qual, seq, opts):
    clip_side = ""
    clip_size = 0
    clip_seq = ""
    clip_pos = -1
    avg_qual = -1
    clip_flag = 1

    if match := re.match(r'^(\d+)S.*M.*?(\d+)S$', cigar):
        if int(match.group(1)) > int(match.group(2)):
            clip_pos = pos
            clip_side = "l"
            clip_size = int(match.group(1))
        else:
            clip_pos = pos - 1
            clip_side = "r"
            clip_size = int(match.group(2))
            for m in re.finditer(r'(\d+)M', cigar):
                clip_pos += int(m.group(1))
            for m in re.finditer(r'(\d+)I', cigar):
                clip_pos -= int(m.group(1))
            for m in re.finditer(r'(\d+)D', cigar):
                clip_pos += int(m.group(1))

    # upstream soft clip only
    elif match := re.match(r'^(\d+)S.*M', cigar):
        clip_pos = pos
        clip_side = "l"
        clip_size = int(match.group(1))

    # downstream soft clip only
    elif match := re.search(r'M.*?(\d+)S', cigar):
        clip_pos = pos - 1
        clip_side = "r"
        clip_size = int(match.group(1))
        for m in re.finditer(r'(\d+)M', cigar):
            clip_pos += int(m.group(1))
        for m in re.finditer(r'(\d+)I', cigar):
            clip_pos -= int(m.group(1))
        for m in re.finditer(r'(\d+)D', cigar):
            clip_pos += int(m.group(1))

    # Check quality of clipped sequence and alignment to reference
    if clip_pos > -1:
        clipped_quals = ""

        if clip_side == "r":
            clip_seq = seq[len(qual) - clip_size - 1:]
            clipped_quals = qual[len(qual) - clip_size - 1:]
        else:
            clip_seq = seq[:clip_size]
            clipped_quals = qual[:clip_size]

        avg_qual_sum = sum(ord(q) - 33 for q in clipped_quals)
        avg_qual_num = len(clipped_quals)
        avg_qual = avg_qual_sum / avg_qual_num if avg_qual_num > 0 else 0

    if avg_qual < 10:
        clip_flag = 0
    if clip_size < opts['min_clipped_seq']:
        clip_flag = 0

    return clip_flag, clip_pos, clip_side, clip_size, clip_seq

def usage(version):
    print("\n")
    print(f"{'%-9s' % 'Program:':<9} gnomit.py")
    print(f"{'%-9s' % 'Version:':<9} {version}")
    print("\n")
    print(f"{'%-9s' % 'Usage:':<9} gnomit.py [options]")
    print("\n")
    print(f"{'%-9s' % 'Options:':<9} {'--input_filename=[filename]':<35} Input alignment file in BAM format")
    print(f"{'%-9s' % '':<9} {'--info_filename=[filename]':<35} Input file with per-sample information (required)")
    print(f"{'%-9s' % '':<9} {'--output_filename=[filename]':<35} Output file (default stdout)")
    print(f"{'%-9s' % '':<9} {'--mask_filename=[filename]':<35} Mask file for reference numts in BED format (optional)")
    print(f"{'%-9s' % '':<9} {'--reference=[filename]':<35} Reference file")
    print(f"{'%-9s' % '':<9} {'--include_mask':<35} Include aberrant reads mapped to mask regions in clustering")
    print(f"{'%-9s' % '':<9} {'--breakpoint':<35} Include soft-clipped reads in likelihood calculation")
    print(f"{'%-9s' % '':<9} {'--len_cluster_include=[integer]':<35} Maximum distance to be included in the cluster (default 600)")
    print(f"{'%-9s' % '':<9} {'--len_cluster_link=[integer]':<35} Maximum distance to link clusters (default 800)")
    print(f"{'%-9s' % '':<9} {'--min_reads_cluster=[integer]':<35} Minimum number of reads to link a cluster (default 1)")
    print(f"{'%-9s' % '':<9} {'--min_evidence=[integer]':<35} Minimum evidence to consider an insertion event for genotyping (default 3)")
    print(f"{'%-9s' % '':<9} {'--min_map_qual=[integer]':<35} Minimum mapping quality for read consideration (default 10)")
    print(f"{'%-9s' % '':<9} {'--max_read_cov=[integer]':<35} Maximum read coverage allowed for breakpoint searching (default 200)")
    print(f"{'%-9s' % '':<9} {'--min_clipped_seq=[integer]':<35} Minimum clipped sequence required to consider as a putative breakpoint (default 5)")
    print(f"{'%-9s' % '':<9} {'--max_num_clipped=[integer]':<35} Maximum number of clipped sequences observed before removing from evidence consideration (default 5)")
    print(f"{'%-9s' % '':<9} {'--read_groups':<35} Stratify analysis to specified read group(s) indicated in info_filename (optional)")
    print(f"{'%-9s' % '':<9} {'--use_priors':<35} Estimate priors using EM framework")
    print(f"{'%-9s' % '':<9} {'--by_chr_dir':<35} If set, expects to find chr specific BAM files info_filename indicated directory")
    print(f"{'%-9s' % '':<9} {'--prefix=[string]':<35} Prepend label in report output")
    print(f"{'%-9s' % '':<9} {'--ucsc':<35} Use UCSC genome formatting (e.g. chrM)")
    print("\n")

def check_options(opt_result, opts, version):
    if not opt_result or opts.get("help"):
        usage(version)
        exit()

    if "input_filename" not in opts and "by_chr_dir" not in opts:
        print("\n***ERROR***\t--input_filename or --by_chr_dir is required\n")
        usage(version)
        exit()
    elif "by_chr_dir" not in opts and not os.path.exists(opts.get("input_filename", "")):
        print("\n***ERROR***\t--input_filename does not exist\n")
        usage(version)
        exit()
    elif "by_chr_dir" in opts and not os.path.isdir(opts["by_chr_dir"]):
        print("\n***ERROR***\t--by_chr_dir does not exist\n")
        usage(version)
        exit()

    if "mask_filename" in opts and not os.path.exists(opts["mask_filename"]):
        print("\n***ERROR***\t--mask_filename does not exist\n")
        usage(version)
        exit()

    if not opts.get("include_mask") and "mask_filename" not in opts:
        print("\n***ERROR***\t--mask_filename is necessary with --include_mask option\n")
        usage(version)
        exit()
