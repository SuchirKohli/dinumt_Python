import sys
import os
from datetime import datetime
import math
import argparse
import subprocess

version = "0.0.23"

# Create the parser
parser = argparse.ArgumentParser()

# Define default values
defaults = {
    'len_cluster_include': 600,
    'len_cluster_link': 800,
    'filter_quality': 50,
    'filter_evidence': 4,
    'filter_depth': 5,
    'min_reads_cluster': 1,
    'min_clipped_seq': 5,
    'max_num_clipped': 5,
    'include_mask': False,
    'min_evidence': 4,
    'min_map_qual': 10,
    'max_read_cov': 200,
    'mask_filename': 'numtS.bed',
    'samtools': 'samtools',
    'prefix': 'numt',
    'len_mt': 16569,
    'ploidy': 2,
    'output_support': False,
    'output_gl': False,
}

# len_cluster_include: The width of the window to consider anchor reads as part of the same cluster.
# len_cluster_link: The width of the window to link two clusters of anchor reads in proper orientation.
# filter_quality: The quality threshold for filtering.
# filter_evidence: The evidence threshold for filtering.
# filter_depth: The depth threshold for filtering.
# min_reads_cluster: The minimum number of reads to form a cluster.
# min_clipped_seq: The minimum length of a clipped sequence.
# max_num_clipped: The maximum number of clipped sequences.
# include_mask: Whether to include masked regions in the analysis.
# min_evidence: The minimum evidence required.
# min_map_qual: The minimum mapping quality required.
# max_read_cov: The maximum read depth at potential breakpoint location.
# mask_filename: The filename of the mask file.
# samtools: The command to run samtools.
# prefix: The prefix for the output files.
# len_mt: The length of the mitochondrial genome. The comment suggests that this should eventually be read from the BAM header.
# ploidy: The ploidy of the organism.
# output_support: Whether to output supporting evidence.
# output_gl: Whether to output genotype likelihoods.

# Add arguments to the parser
for opt, default in defaults.items():
    if isinstance(default, bool):
        parser.add_argument(f'--{opt}', action='store_true', default=default)
    else:
        parser.add_argument(f'--{opt}', type=type(default), default=default)

# Parse the arguments
args = parser.parse_args()

# Convert the Namespace to a dictionary
opts = vars(args)

# def check_options(opt_result, opts, version):
    # Assuming checkOptions is a function that takes three parameters

seq_num = 0
seq_hash = {}

sorted_hash = {}
readgroup_hash = {}

i = 1
infile_hash = {}
group_hash = {}
outfile_hash = {}
mask_hash = {}
mt_hash = {}

if 'read_groups' in opts:
    rgs = opts['read_groups'].split(',')
    readgroup_hash = {rg: 1 for rg in rgs}

if 'mt_names' in opts:
    mts = opts['mt_names'].split(',')
    mt_hash = {mt: 1 for mt in mts}

def score_data(outfile_hash):
    print("entering scoreData()")
    
    for group in outfile_hash:
        # This loop iterates over each group in the outfile_hash.
        print("Group:", group)
        sumGP = 0
        numRef = outfile_hash[group]["numRefRP"]
        numAlt = outfile_hash[group]["numAltRP"]
        if outfile_hash[group]["numAltSR"] > 0:
            numRef += outfile_hash[group]["numRefSR"]
            numAlt += outfile_hash[group]["numAltSR"]
        print("\t", numRef,"\t", numAlt)
        
        # This Loop iterates over genotypes from 0 to the specified ploidy.
        for g in range(opts["ploidy"] + 1):
            geno = opts["ploidy"] - g      #need to reverse as calculation is reference allele based
            if numAlt + numRef > 0 and 1 / opts["ploidy"] ** (numAlt + numRef) > 0:
                outfile_hash[group]["gl"][geno] = calc_gl(opts["ploidy"], g, numAlt + numRef, numRef, outfile_hash[group]["avgQ"])
                outfile_hash[group]["gp"][geno] = 10 ** outfile_hash[group]["gl"][geno]
                sumGP += outfile_hash[group]["gp"][geno]
                print("\t", geno, outfile_hash[group]["gl"][geno], outfile_hash[group]["gp"][geno])
        
        print("\tsumGP:", sumGP)
        if sumGP == 0:
            # If sumGP is zero, it sets all genotype likelihoods and probabilities to zero, and sets the genotype (gt) to "./." and the filter (ft) to "NC".
            for geno in range(opts["ploidy"] + 1):
                outfile_hash[group]["pl"][geno] = 0
                outfile_hash[group]["gl"][geno] = 0
            outfile_hash[group]["gq"] = 0
            outfile_hash[group]["gt"] = "./."
            outfile_hash[group]["ft"] = "NC"
        else:
            # If sumGP is not zero, it normalizes the genotype probabilities by sumGP, calculates Phred-scaled likelihoods (pl), and identifies the genotype with the maximum probability.
            maxGP = 0
            maxGeno = 0
            for geno in range(opts["ploidy"] + 1):
                if outfile_hash[group]["gp"][geno] == 0:
                    outfile_hash[group]["gp"][geno] = 1e-200
                outfile_hash[group]["gp"][geno] /= sumGP
                outfile_hash[group]["pl"][geno] = int(-10 * math.log10(outfile_hash[group]["gp"][geno]))
                if outfile_hash[group]["gp"][geno] > maxGP:
                    maxGP = outfile_hash[group]["gp"][geno]
                    maxGeno = geno
            
            maxGP = 1 - outfile_hash[group]["gp"][0]

            # It then calculates the genotype quality (gq) as the Phred-scaled probability of the most likely genotype not being the true genotype.
            if 1 - maxGP == 0:
                outfile_hash[group]["gq"] = 199
            else:
                outfile_hash[group]["gq"] = int(-10 * math.log10(1 - maxGP))
            
            # It sets the genotype (gt) based on the most likely genotype.
            gt = "0/0"
            if maxGeno == 1:
                gt = "0/1"
            elif maxGeno == 2:
                gt = "1/1"
            
            outfile_hash[group]["gt"] = gt
            
            # It applies several filters based on genotype quality, evidence, and depth, and sets the filter (ft) to either the failed filters or "PASS" if all filters are passed.
            filters = []
            if outfile_hash[group]["gq"] < opts["filter_quality"]:
                filters.append("q" + str(opts["filter_quality"]))
            if numAlt < opts["filter_evidence"]:
                filters.append("e" + str(opts["filter_evidence"]))
            if numAlt + numRef < opts["filter_depth"]:
                filters.append("d" + str(opts["filter_depth"]))
            
            outfile_hash[group]["ft"] = ";".join(filters) if filters else "PASS"
    
    print("exiting scoreData()")

def get_date():
    now = datetime.now()
    year = now.year
    month = now.month
    day = now.day

    fmonth = f"{month:02d}"
    fday = f"{day:02d}"

    return f"{year}{fmonth}{fday}"


def report(outfile_hash):
    print("entering report()") if opts.get('verbose', False)

    # Open output file
    if 'output_filename' in opts and opts['output_filename']:
        with open(opts['output_filename'], 'w') as foutname1:
            pass  
    else:
        foutname1 = sys.stdout

    if opts.get('output_support', False):
        # Open support file
        with open(opts['support_filename'], 'w') as support1:
            pass  

    filedate = get_date()
    header = f'''##fileformat=VCFv4.1
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS:MT,Description="Nuclear Mitochondrial Insertion">
##FILTER=<ID=q{opts['filter_quality']},Description="Phred-scaled quality filter">
##FILTER=<ID=e{opts['filter_evidence']},Description="Support reads filter">
##FILTER=<ID=d{opts['filter_depth']},Description="Sequence depth filter">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description="Copy number likelihoods">
##FORMAT=<ID=CNL0,Number=.,Type=Float,Description="Copy number likelihoods with no frequency prior">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=MSTART,Number=1,Type=Flag,Description="Mitochondrial start coordinate of inserted sequence">
##INFO=<ID=MEND,Number=1,Type=Flag,Description="Mitochondrial end coordinate of inserted sequence">
##INFO=<ID=MLEN,Number=1,Type=Flag,Description="Estimated length of mitochondrial insert">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##fileDate={filedate}
##reference={opts['reference']}
##source=dinumt-{version}
'''

    foutname1.write(header)

    vars_sorted = sorted(outfile_hash.keys(), key=lambda var: (outfile_hash[var]['chr'], outfile_hash[var]['leftBkpt']))
    if opts.get('output_gl', False):
        foutname1.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + opts['prefix'] + "\n")
    else:
        foutname1.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    index = 1
    for group in vars_sorted:
        print(f"{group} in report()") # if opts.get('verbose', False)
        if outfile_hash[group]['gt'] == "0/0" or outfile_hash[group]['gt'] == "./.":
            print("\thomref or nc, skipping") # if opts.get('verbose', False)
            continue

        chrom = outfile_hash[group]['chr']
        if not opts.get('ucsc', False) and not opts.get('ensembl', False):
            chrom = chrom.replace("chr", "")

        identifier = f"{opts['prefix']}_{index}"
        alt = "<INS:MT>"
        qual = outfile_hash[group]['gq']
        ftr = outfile_hash[group]['ft']
        info = {
            'IMPRECISE': None,
            'CIPOS': f"0,{outfile_hash[group]['rightBkpt'] - outfile_hash[group]['leftBkpt'] + 1}",
            'CIEND': f"-{outfile_hash[group]['rightBkpt'] - outfile_hash[group]['leftBkpt'] + 1},0",
            'END': outfile_hash[group]['rightBkpt'],
            'SVTYPE': 'INS'
        }

        if outfile_hash[group]['m_len'] != "NA":
            info['MSTART'] = outfile_hash[group]['l_m_pos']
            info['MEND'] = outfile_hash[group]['r_m_pos']
            info['MLEN'] = outfile_hash[group]['m_len']

        refline = subprocess.getoutput(f"{opts['samtools']} faidx {opts['reference']} {chrom}:{outfile_hash[group]['leftBkpt']}-{outfile_hash[group]['leftBkpt']}")
        ref = refline.split('\n')[1] if refline else "N"

        format_str = "GT:FT:GL:GQ:PL"
        info_str = ";".join([f"{key}={info[key]}" if info[key] is not None else key for key in sorted(info.keys())])

        if opts.get('output_gl', False):
            gls = [f"{geno:.2f}" for geno in outfile_hash[group]['gl'].values()]
            pls = [str(pl) for pl in outfile_hash[group]['pl'].values()]
            gl = ",".join(gls)
            pl = ",".join(pls)
            foutname1.write(f"{chrom}\t{outfile_hash[group]['leftBkpt']}\t{identifier}\t{ref}\t{alt}\t{qual}\t{ftr}\t{info_str}\t{format_str}\t{outfile_hash[group]['gt']}:{outfile_hash[group]['ft']}:{gl}:{outfile_hash[group]['gq']}:{pl}\n")
        else:
            foutname1.write(f"{chrom}\t{outfile_hash[group]['leftBkpt']}\t{identifier}\t{ref}\t{alt}\t{qual}\t{ftr}\t{info_str}\n")

        if opts.get('output_support', False):
            with open(opts['support_filename'], 'a') as support1:
                support1.write(f"{outfile_hash[group]['support']}")

        index += 1

    if opts.get('output_support', False):
        support1.close()

    foutname1.close()
    if opts.get('verbose', False):
        print("exiting report()") 
    

def calc_gl(m, g, k, l, e):
    if opts.get('verbose', False):
        print("in calcGl():") 
    if opts.get('verbose', False):
        print(f"\t{m}\t{g}\t{k}\t{l}\t{e}") 

    if 1 / m**k <= 0:
        raise ValueError("Problem in calcGL 1, \t{m}\t{g}\t{k}\t{l}\t{e}")

    gl = math.log10(1 / m**k)

    if ((m - g) * e + (1 - e) * g) <= 0:
        raise ValueError("Problem in calcGL 2, \t{m}\t{g}\t{k}\t{l}\t{e}")

    for i in range(1, l + 1):
        gl += math.log10((m - g) * e + (1 - e) * g)

    if ((m - g) * (1 - e) + g * e) <= 0:
        raise ValueError("Problem in calcGL 3, \t{m}\t{g}\t{k}\t{l}\t{e}")

    for i in range(l + 1, k + 1):
        gl += math.log10((m - g) * (1 - e) + g * e)

    return gl

def get_input(infile_hash, readgroup_hash, mask_hash, mt_hash):
    input_lines = []
    print("Reading input files...") # if opts.get('verbose', False)

    # Open input file
    if opts.get('by_chr_dir', False):
        if opts.get('mt_names', None):
            for mt_name in mt_hash.keys():
                input_lines.append(f"samtools view {opts['by_chr_dir']}/{mt_name}.*bam |")
        elif opts.get('ucsc', False):
            input_lines.append(f"samtools view {opts['by_chr_dir']}/chrM.*bam |")
        elif opts.get('ensembl', False):
            input_lines.append(f"samtools view {opts['by_chr_dir']}/chrMT.*bam |")
        else:
            input_lines.append(f"samtools view {opts['by_chr_dir']}/MT*bam |")
    else:
        if opts.get('mt_names', None):
            for mt_name in mt_hash.keys():
                input_lines.append(f"samtools view {opts['input_filename']} {mt_name} |")
        elif opts.get('ucsc', False):
            input_lines.append(f"samtools view {opts['input_filename']} chrM |")
        elif opts.get('ensembl', False):
            input_lines.append(f"samtools view {opts['input_filename']} chrMT |")
        else:
            input_lines.append(f"samtools view {opts['input_filename']} MT |")

    # Input mask coordinates
    if opts.get('include_mask', False):
        with open(opts['mask_filename'], 'r') as mask1:
            for line in mask1:
                parts = line.strip().split('\t')
                chr_name = parts[0].replace("chr", "")
                mask_start, mask_end = int(parts[1]), int(parts[2])
                mask_hash[chr_name][mask_start] = mask_end

                if opts.get('include_mask', False):
                    if opts.get('by_chr_dir', False):
                        if opts.get('ucsc', False) or opts.get('ensembl', False):
                            chr_name = "chr" + chr_name
                        input_lines.append(f"samtools view {opts['by_chr_dir']}/{chr_name}*bam {chr_name}:{mask_start}-{mask_end} |")
                    else:
                        if opts.get('ucsc', False) or opts.get('ensembl', False):
                            chr_name = "chr" + chr_name
                        input_lines.append(f"samtools view {opts['input_filename']} {chr_name}:{mask_start}-{mask_end} |")

    for input_line in input_lines:
        print(f"command: {input_line}") # if opts.get('verbose', False)
        with subprocess.Popen(input_line, shell=True, stdout=subprocess.PIPE) as proc:
            for line1 in proc.stdout:
                seq_num = i
                line1 = line1.decode('utf-8').strip()
                parts = line1.split('\t')
                qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = parts[:11]
                pnextend = int(pnext) + len(seq)

                read_group = None
                if 'RG:Z:' in line1:
                    read_group = line1.split('RG:Z:')[1].split('\t')[0]

                if opts.get('read_groups', False) and read_group is None:
                    continue
                elif opts.get('read_groups', False) and read_group not in readgroup_hash:
                    continue

                if rnext == '=' or rnext == '*':
                    continue

                if opts.get('mt_names', None):
                    if rname not in mt_hash and rnext in mt_hash:
                        continue
                else:
                    if 'M' not in rname and 'M' in rnext:
                        continue

                dnext = 1 if int(flag) & 32 else 0
                direction = 1 if int(flag) & 16 else 0

                is_mask_overlap = 0
                if opts.get('include_mask', False):
                    for mask_start, mask_end in mask_hash[rnext].items():
                        if pnext >= mask_start and pnext <= mask_end:
                            is_mask_overlap = 1
                            break
                        elif pnextend >= mask_start and pnextend <= mask_end:
                            is_mask_overlap = 1
                            break
                        elif pnext <= mask_start and pnextend >= mask_end:
                            is_mask_overlap = 1
                            break

                if is_mask_overlap:
                    continue

                with subprocess.Popen(f"samtools view {' '.join(input_line.split()[2:])} |", shell=True, stdout=subprocess.PIPE) as sam_proc:
                    c_mapq = 0
                    cnext = 0
                    for sam_line in sam_proc.stdout:
                        sam_line = sam_line.decode('utf-8').strip().split('\t')
                        if sam_line[0] != qname:
                            continue
                        c_mapq = int(sam_line[4])
                        cnext = sam_line[5]

                if c_mapq < opts['min_map_qual']:
                    continue

                infile_hash[seq_num] = {
                    'group': 0,
                    'seq': seq,
                    'dir': direction,
                    'qname': qname,
                    'rname': rname,
                    'pos': pos,
                    'cigar': cigar,
                    'cnext': cnext,
                    'rnext': rnext,
                    'pnext': pnext,
                    'dnext': dnext,
                    'tlen': tlen,
                    'qual': qual,
                    'line': line1
                }

                i += 1

    return infile_hash

def find_breakpoint(outfile_hash, readgroup_hash, mask_hash):
    print("entering findBreakpoint()")# if opts.get('verbose', False)
    for group, values in outfile_hash.items():
        print(f"\tAssessing group {group} for breakpoints") # if opts.get('verbose', False)
        l_pos = values['l_pos']
        r_pos = values['r_pos']
        chro = values['chr']
        values['numAltSR'] = 0
        values['numRefSR'] = 0
        values['leftBkpt'] = l_pos
        values['rightBkpt'] = r_pos

        # Compare to masked regions
        is_mask_overlap = False
        if opts.get('include_mask', False):
            for mask_start, mask_end in mask_hash[chro].items():
                if (l_pos >= mask_start and l_pos <= mask_end) or (r_pos >= mask_start and r_pos <= mask_end) or (l_pos <= mask_start and r_pos >= mask_end):
                    is_mask_overlap = True
                    break

        if is_mask_overlap:
            continue

        clipped_pos = {}

        # Open input file
        if opts.get('by_chr_dir', False):
            with subprocess.Popen(f"samtools view {opts['by_chr_dir']}/{chro}.*bam {chro}:{l_pos}-{r_pos} |", shell=True, stdout=subprocess.PIPE) as sam_proc:
                process_sam_output(sam_proc, values, group)
        else:
            with subprocess.Popen(f"samtools view {opts['input_filename']} {chro}:{l_pos}-{r_pos} |", shell=True, stdout=subprocess.PIPE) as sam_proc:
                process_sam_output(sam_proc, values, group)

        if group not in outfile_hash:
            continue

        bkpts = {}
        num_bkpts = 0
        for c_pos, count in clipped_pos.items():
            if count > 1:
                bkpts[c_pos] = count
                num_bkpts += 1

        num_bkpt_support = 0
        num_ref_support = 0
        left_bkpt = l_pos
        right_bkpt = r_pos

        if num_bkpts > 0 and len(clipped_pos) <= opts.get('max_num_clipped', 0):

            # Take two most prevalent breaks for now
            sorted_bkpts = sorted(bkpts.keys(), key=lambda x: bkpts[x], reverse=True)
            if num_bkpts == 1:
                left_bkpt = sorted_bkpts[0]
                right_bkpt = left_bkpt + 1
                num_bkpt_support = bkpts[sorted_bkpts[0]]
                num_ref_support = values['cnt'][sorted_bkpts[0]] - num_bkpt_support
            else:
                if sorted_bkpts[0] < sorted_bkpts[1]:
                    left_bkpt = sorted_bkpts[0]
                    right_bkpt = sorted_bkpts[1]
                else:
                    left_bkpt = sorted_bkpts[1]
                    right_bkpt = sorted_bkpts[0]
                num_bkpt_support = bkpts[sorted_bkpts[0]] + bkpts[sorted_bkpts[1]]
                num_ref_support = values['cnt'][sorted_bkpts[0]] + values['cnt'][sorted_bkpts[1]] - num_bkpt_support

        values['leftBkpt'] = left_bkpt
        values['rightBkpt'] = right_bkpt
        values['numAltSR'] = num_bkpt_support
        values['numRefSR'] = num_ref_support

    print("exiting findBreakpoints()") if opts.get('verbose', False)

def seq_cluster(infile_hash):
    k = 0
    d = {
        0: {'k': 0, 'pnext': 0, 'last': [], 'rnext': 0},
        1: {'k': 0, 'pnext': 0, 'last': [], 'rnext': 0}
    }

    sorted_keys = sorted(infile_hash.keys(), key=lambda x: (infile_hash[x]['rnext'], infile_hash[x]['pnext']))
    print(f"{len(sorted_keys)} total reads to process for clustering")# if opts.get('verbose', False)

    for c_seq_num in sorted_keys:
        c_pnext = infile_hash[c_seq_num]['pnext']
        c_dnext = infile_hash[c_seq_num]['dnext']
        c_rnext = infile_hash[c_seq_num]['rnext']
        c_qname = infile_hash[c_seq_num]['qname']

        print(c_dnext)# if opts.get('verbose', False)

        if c_pnext - d[c_dnext]['pnext'] > opts.get('len_cluster_include', 0) or d[c_dnext]['k'] == 0 or c_rnext != d[c_dnext]['rnext']:
            if d[c_dnext]['k'] > 0 and len(d[c_dnext]['last']) < opts.get('min_reads_cluster', 0):
                for seq_num in d[c_dnext]['last']:
                    infile_hash.pop(seq_num, None)

            k += 1
            infile_hash[c_seq_num]['group'] = k
            d[c_dnext]['k'] = k
            d[c_dnext]['last'] = []
        else:
            infile_hash[c_seq_num]['group'] = d[c_dnext]['k']

        print(f"{d[c_dnext]['k']}\t{c_qname}\t{c_rnext}\t{c_pnext}") if opts.get('verbose', False)

        d[c_dnext]['last'].append(c_seq_num)
        d[c_dnext]['pnext'] = c_pnext
        d[c_dnext]['rnext'] = c_rnext

    return infile_hash

def link_cluster(infile_hash):
    # this can link multiple F's to a single leftmost R
    sorted_keys = sorted(infile_hash.keys(), key=lambda x: (infile_hash[x]['rnext'], infile_hash[x]['pnext']))
    print(f"{len(sorted_keys)} total reads to process for linking clusters") if opts.get('verbose', False)

    for c in range(len(sorted_keys)):
        c_seq_num = sorted_keys[c]
        infile_hash[c_seq_num]['link'] = 0
        c_pnext = infile_hash[c_seq_num]['pnext']
        c_dnext = infile_hash[c_seq_num]['dnext']
        c_rnext = infile_hash[c_seq_num]['rnext']
        c_dir = infile_hash[c_seq_num]['dir']

        if c_dnext == 1:
            continue

        for d in range(c + 1, len(sorted_keys)):
            d_seq_num = sorted_keys[d]
            d_pnext = infile_hash[d_seq_num]['pnext']
            d_dnext = infile_hash[d_seq_num]['dnext']
            d_rnext = infile_hash[d_seq_num]['rnext']
            d_dir = infile_hash[d_seq_num]['dir']

            if d_dnext == 0:
                continue

            delta = d_pnext - c_pnext

            if (
                delta < opts.get('len_cluster_link', 0)
                and c_rnext == d_rnext
                and (
                    (c_dir == 0 and c_dnext == 1 and d_dnext == 0 and d_dir == 1)
                    or (c_dir == 0 and c_dnext == 0 and d_dnext == 1 and d_dir == 1)
                    or (c_dir == 1 and c_dnext == 0 and d_dnext == 1 and d_dir == 0)
                )
            ):
                infile_hash[c_seq_num]['link'] = infile_hash[d_seq_num]['group']
                infile_hash[d_seq_num]['link'] = infile_hash[c_seq_num]['group']
            elif infile_hash[c_seq_num]['link'] is None:
                infile_hash[c_seq_num]['link'] = 0

def map_cluster(infile_hash, outfile_hash, readgroup_hash):
    l_linked_groups = {}
    linked_group_pnext = {}

    for key in sorted(infile_hash.keys(), key=lambda x: infile_hash[x]['group']):
        if (
            infile_hash[key]['group'] > 0
            and infile_hash[key]['dnext'] == 0
            and infile_hash[key]['link'] > 0
        ):
            l_linked_groups[infile_hash[key]['group']] = infile_hash[key]['link']

    i = 1
    for group, link in l_linked_groups.items():
        print(f"group = {group} ;; link = {link}") if opts.get('verbose', False)
        l_rnext = [infile_hash[key]['rnext'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_rname = [infile_hash[key]['rname'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_pos = [infile_hash[key]['pos'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_dir = [infile_hash[key]['dir'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_qname = [infile_hash[key]['qname'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_pnext = [infile_hash[key]['pnext'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_dnext = [infile_hash[key]['dnext'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_cigar = [infile_hash[key]['cigar'] for key in infile_hash if infile_hash[key]['group'] == group]
        l_cnext = [infile_hash[key]['cnext'] for key in infile_hash if infile_hash[key]['group'] == group]
        r_rnext = [infile_hash[key]['rnext'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_pos = [infile_hash[key]['pos'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_dir = [infile_hash[key]['dir'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_rname = [infile_hash[key]['rname'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_pnext = [infile_hash[key]['pnext'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_dnext = [infile_hash[key]['dnext'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_qname = [infile_hash[key]['qname'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_cigar = [infile_hash[key]['cigar'] for key in infile_hash if infile_hash[key]['group'] == link]
        r_cnext = [infile_hash[key]['cnext'] for key in infile_hash if infile_hash[key]['group'] == link]

        rc_pnext = []
        lc_pnext = []
        chr_val = l_rnext[0]

        for c in range(len(r_cnext)):
            cigar = r_cnext[c]
            rc_pnext.append(r_pnext[c])

            while "(\d+)M" in cigar:
                match = re.search("(\d+)M", cigar)
                rc_pnext[c] += int(match.group(1))
                cigar = cigar[match.end():]
            while "(\d+)N" in cigar:
                match = re.search("(\d+)N", cigar)
                rc_pnext[c] += int(match.group(1))
                cigar = cigar[match.end():]
            while "(\d+)D" in cigar:
                match = re.search("(\d+)D", cigar)
                rc_pnext[c] += int(match.group(1))
                cigar = cigar[match.end():]

        for c in range(len(l_cnext)):
            cigar = l_cnext[c]
            lc_pnext.append(l_pnext[c])

            while "(\d+)M" in cigar:
                match = re.search("(\d+)M", cigar)
                lc_pnext[c] += int(match.group(1))
                cigar = cigar[match.end():]
            while "(\d+)N" in cigar:
                match = re.search("(\d+)N", cigar)
                lc_pnext[c] += int(match.group(1))
                cigar = cigar[match.end():]
            while "(\d+)D" in cigar:
                match = re.search("(\d+)D", cigar)
                lc_pnext[c] += int(match.group(1))
                cigar = cigar[match.end():]

        s_l_pnext = sorted(l_pnext)
        s_r_pnext = sorted(r_pnext)
        s_lc_pnext = sorted(lc_pnext)
        s_rc_pnext = sorted(rc_pnext)

        l_brk_point = s_l_pnext[-1]
        r_brk_point = s_rc_pnext[0]

        win_l_s = s_l_pnext[0]
        win_l_e = s_lc_pnext[-1]
        win_r_s = s_r_pnext[0]
        win_r_e = s_rc_pnext[-1]

        linked_group_pnext[i] = {'l_group': group, 'r_group': link}

        if l_brk_point > r_brk_point:
            l_brk_point, r_brk_point = r_brk_point, l_brk_point

        command_left = ""
        command_right = ""
        if opts.get('by_chr_dir'):
            command_left = f"samtools view {opts['by_chr_dir']}/{chr_val}.*bam {chr_val}:{win_l_s}-{win_l_e} |"
            command_right = f"samtools view {opts['by_chr_dir']}/{chr_val}.*bam {chr_val}:{win_r_s}-{win_r_e} |"
        else:
            command_left = f"samtools view {opts['input_filename']} {chr_val}:{win_l_s}-{win_l_e} |"
            command_right = f"samtools view {opts['input_filename']} {chr_val}:{win_r_s}-{win_r_e} |"

        num_ref_rp = 0
        num_alt_rp = len(l_pnext) + len(r_pnext)
        sum_e = 0

        with os.popen(command_left) as sam_file:
            for line in sam_file:
                qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip().split('\t')
                map_e = 10 ** (-1 * int(mapq) / 10)
                read_group = re.search(r'RG:Z:(\S+)', line)
                if opts.get('read_groups') and not read_group:
                    continue
                elif opts.get('read_groups') and read_group.group(1) not in readgroup_hash:
                    continue
                if int(mapq) < opts.get('min_map_qual', 0):
                    continue
                dir_val = 0 if int(flag) & 16 == 0 else 1
                if int(pos) >= win_l_s and int(pos) <= win_l_e and dir_val == 0:
                    num_ref_rp += 1
                    sum_e += map_e

        with os.popen(command_right) as sam_file:
            for line in sam_file:
                qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip().split('\t')
                map_e = 10 ** (-1 * int(mapq) / 10)
                read_group = re.search(r'RG:Z:(\S+)', line)
                if opts.get('read_groups') and not read_group:
                    continue
                elif opts.get('read_groups') and read_group.group(1) not in readgroup_hash:
                    continue
                if int(mapq) < opts.get('min_map_qual', 0):
                    continue
                dir_val = 0 if int(flag) & 16 == 0 else 1
                if int(pos) >= win_r_s and int(pos) <= win_r_e and dir_val == 1:
                    num_ref_rp += 1
                    sum_e += map_e

        avg_q = sum_e / num_ref_rp
        num_ref_rp -= num_alt_rp

        l_m_min = 1e10
        l_m_max = 0
        l_m_min_i = -1
        l_m_max_i = -1
        l_n_dir = -1
        l_m_dir = -1

        for i in range(len(l_qname)):
            if 'M' not in l_rname[i]:
                continue
            if l_pos[i] < l_m_min:
                l_m_min = l_pos[i]
                l_m_min_i = i
            if l_pos[i] > l_m_max:
                l_m_max = l_pos[i]
                l_m_max_i = i
            l_n_dir = l_dnext[i]
            l_m_dir = l_dir[i]

        r_m_min = 1e10
        r_m_max = 0
        r_m_min_i = -1
        r_m_max_i = -1
        r_n_dir = -1
        r_m_dir = -1

        for i in range(len(r_qname)):
            if 'M' not in r_rname[i]:
                continue
            if r_pos[i] < r_m_min:
                r_m_min = r_pos[i]
                r_m_min_i = i
            if r_pos[i] > r_m_max:
                r_m_max = r_pos[i]
                r_m_max_i = i
            r_n_dir = r_dnext[i]
            r_m_dir = r_dir[i]

        outfile_hash[group] = {
            'l_m_pos': "NA",
            'r_m_pos': "NA",
            'm_len': "NA",
            'avgQ': avg_q,
            'l_pos': l_brk_point,
            'r_pos': r_brk_point,
            'chr': chr_val,
            'numRefRP': num_ref_rp,
            'numAltRP': num_alt_rp,
            'support': "",
        }

        if l_m_dir > -1 and r_m_dir > -1:
            if l_n_dir == 0 and l_m_dir == 1 and r_m_dir == 0 and r_n_dir == 1:
                cigar = r_cigar[r_m_max_i]
                while "(\d+)M" in cigar:
                    match = re.search("(\d+)M", cigar)
                    r_m_max += int(match.group(1))
                    cigar = cigar[match.end():]
                while "(\d+)N" in cigar:
                    match = re.search("(\d+)N", cigar)
                    r_m_max += int(match.group(1))
                    cigar = cigar[match.end():]
                while "(\d+)D" in cigar:
                    match = re.search("(\d+)D", cigar)
                    r_m_max += int(match.group(1))
                    cigar = cigar[match.end():]
            elif l_n_dir == 0 and l_m_dir == 0 and r_m_dir == 1 and r_n_dir == 1:
                cigar = l_cigar[l_m_max_i]
                while "(\d+)M" in cigar:
                    match = re.search("(\d+)M", cigar)
                    l_m_max += int(match.group(1))
                    cigar = cigar[match.end():]
                while "(\d+)N" in cigar:
                    match = re.search("(\d+)N", cigar)
                    l_m_max += int(match.group(1))
                    cigar = cigar[match.end():]
                while "(\d+)D" in cigar:
                    match = re.search("(\d+)D", cigar)
                    l_m_max += int(match.group(1))
                    cigar = cigar[match.end():]
            
            outfile_hash[group]['l_m_pos'] = l_m_min
            outfile_hash[group]['r_m_pos'] = r_m_max
            if r_m_max > l_m_min:
                outfile_hash[group]['l_m_pos'] = l_m_min
                outfile_hash[group]['r_m_pos'] = r_m_max
            else:
                outfile_hash[group]['l_m_pos'] = r_m_min
                outfile_hash[group]['r_m_pos'] = l_m_max
            
            #currently assumes smallest sequence possible due to circular nature of mt dna
            outfile_hash[group]['m_len'] = outfile_hash[group]['r_m_pos'] - outfile_hash[group]['l_m_pos'] + 1
            lenAlt = opts['len_mt'] - outfile_hash[group]['r_m_pos'] + outfile_hash[group]['l_m_pos'] + 1
            if lenAlt < outfile_hash[group]['m_len']:
                outfile_hash[group]['m_len'] = lenAlt

        outfile_hash[group]['avgQ'] = avgQ
        outfile_hash[group]['l_pos'] = l_brk_point
        outfile_hash[group]['r_pos'] = r_brk_point
        outfile_hash[group]['chr'] = chr
        outfile_hash[group]['numRefRP'] = numRefRP
        outfile_hash[group]['numAltRP'] = numAltRP

        if opts['output_support']:
            for key in [k for k, v in infile_hash.items() if v['group'] == group]:
                outfile_hash[group]['support'] += infile_hash[key]['line'] + "\n"
            for key in [k for k, v in infile_hash.items() if v['group'] == link]:
                outfile_hash[group]['support'] += infile_hash[key]['line'] + "\n"

        i += 1

        if opts['verbose']:
            print(f"\t{chr}\t{outfile_hash[group]['l_pos']}\t{outfile_hash[group]['r_pos']}\t{outfile_hash[group]['numRefRP']}\t{outfile_hash[group]['numAltRP']}\t{outfile_hash[group]['avgQ']}\t{outfile_hash[group]['l_m_pos']}\t{outfile_hash[group]['r_m_pos']}\t{outfile_hash[group]['m_len']}")            

import re

def get_soft_clip_info(pos, cigar, qual, opts):
    clip_side = ""
    clip_size = 0
    c_pos = -1
    avg_qual = -1

    # Check for soft clipping at both ends
    if cigar.find('S') > 0 and cigar.find('M', cigar.find('S')) > 0:
        split_cigar = cigar.split('M')
        left_clip, right_clip = split_cigar[0], split_cigar[1].split('S')[1]
        
        if int(left_clip) > int(right_clip):
            c_pos = pos
            clip_side = "l"
            clip_size = int(left_clip[:-1])  # exclude 'S' from the size
        else:
            c_pos = pos - 1
            clip_side = "r"
            clip_size = int(right_clip)

            for m in re.finditer(r'(\d+)M', cigar):
                c_pos += int(m.group(1))
            for i in re.finditer(r'(\d+)I', cigar):
                c_pos -= int(i.group(1))
            for d in re.finditer(r'(\d+)D', cigar):
                c_pos += int(d.group(1))

    # Check for only upstream soft clip
    elif cigar.find('S') > 0 and cigar.find('M') == -1:
        c_pos = pos
        clip_side = "l"
        clip_size = int(cigar.split('S')[0])

    # Check for only downstream soft clip
    elif cigar.find('M') > 0 and cigar.find('S', cigar.find('M')) > 0:
        split_cigar = cigar.split('S')
        clip_size = int(split_cigar[1].split('M')[0])
        c_pos = pos - 1
        clip_side = "r"

        for m in re.finditer(r'(\d+)M', cigar):
            c_pos += int(m.group(1))
        for i in re.finditer(r'(\d+)I', cigar):
            c_pos -= int(i.group(1))
        for d in re.finditer(r'(\d+)D', cigar):
            c_pos += int(d.group(1))

    if c_pos > -1:
        clipped_quals = ""

        if clip_side == "r":
            clipped_quals = qual[-clip_size - 1:]
        else:
            clipped_quals = qual[:clip_size]

        avg_qual_sum = sum(ord(q) - 33 for q in clipped_quals)
        avg_qual_num = len(clipped_quals)
        avg_qual = avg_qual_sum / avg_qual_num if avg_qual_num > 0 else 0

    if avg_qual < 10:
        c_pos = -1
    if clip_size < opts['min_clipped_seq']:
        c_pos = -1

    return c_pos, clip_side, clip_size

def usage(version):
    print("\n")
    print(f"{'Program:':<9} {'dinumt.py':<35}")
    print(f"{'Version:':<9} {version}")
    print("\n")
    print(f"{'Usage:':<9} {'dinumt.py [options]'}")
    print("\n")
    print(f"{'Options:':<9} {'--input_filename=[filename]':<35} {'Input alignment file in BAM format'}")
    print(f"{'':<9} {'--output_filename=[filename]':<35} {'Output file (default stdout)'}")
    print(f"{'':<9} {'--mask_filename=[filename]':<35} {'Mask file for reference numts in BED format (optional)'}")
    print(f"{'':<9} {'--include_mask':<35} {'Include aberrant reads mapped to mask regions in clustering'}")
    print(f"{'':<9} {'--len_cluster_include=[integer]':<35} {'Maximum distance to be included in cluster (default 600)'}")
    print(f"{'':<9} {'--len_cluster_link=[integer]':<35} {'Maximum distance to link clusters (default 800)'}")
    print(f"{'':<9} {'--min_reads_cluster=[integer]':<35} {'Minimum number of reads to link a cluster (default 1)'}")
    print(f"{'':<9} {'--min_evidence=[integer]':<35} {'Minimum evidence to consider an insertion event (default 4)'}")
    print(f"{'':<9} {'--min_map_qual=[integer]':<35} {'Minimum mapping quality for read consideration (default 10)'}")
    print(f"{'':<9} {'--max_read_cov=[integer]':<35} {'Maximum read coverage allowed for breakpoint searching (default 200)'}")
    print(f"{'':<9} {'--min_clipped_seq=[integer]':<35} {'Minimum clipped sequence required to consider as putative breakpoint (default 5)'}")
    print(f"{'':<9} {'--max_num_clipped=[integer]':<35} {'Maximum number of clipped sequences observed before removing from evidence consideration (default 5)'}")
    print(f"{'':<9} {'--read_groups=[read_group1],...':<35} {'Limit analysis to specified read group(s)'}")
    print(f"{'':<9} {'--mt_names=[mt_name1],...':<35} {'Limit analysis to specified mitochondrial sequence names'}")
    print(f"{'':<9} {'--by_chr_dir=[directory]':<35} {'If set, expects to find chr specific BAM files in indicated directory'}")
    print(f"{'':<9} {'--prefix=[string]':<35} {'Prepend label in report output'}")
    print(f"{'':<9} {'--ucsc':<35} {'Use UCSC genome formatting (e.g. chrM)'}")
    print(f"{'':<9} {'--ensembl':<35} {'Use Ensembl genome formatting (e.g. chrMT)'}")
    print(f"{'':<9} {'--output_gl':<35} {'Output genotype likelihood information'}")
    print("\n")


def check_options(opt_result, opts, version):
    if not opt_result or opts.get('help'):
        usage(version)
        sys.exit()

    if not opts.get('input_filename') and not opts.get('by_chr_dir'):
        print("\n***ERROR***\t--input_filename or --by_chr_dir is required\n")
        usage(version)
        sys.exit()
    elif not opts.get('by_chr_dir') and not os.path.exists(opts.get('input_filename')):
        print("\n***ERROR***\t--input_filename does not exist\n")
        usage(version)
        sys.exit()
    elif opts.get('by_chr_dir') and not os.path.isdir(opts.get('by_chr_dir')):
        print("\n***ERROR***\t--by_chr_dir does not exist\n")
        usage(version)
        sys.exit()
    
    if opts.get('mask_filename') and not os.path.exists(opts.get('mask_filename')):
        print("\n***ERROR***\t--mask_filename does not exist\n")
        usage(version)
        sys.exit()
    
    if not opts.get('reference'):
        print("\n***ERROR***\t--reference is required\n")
        usage(version)
        sys.exit()
    
    if opts.get('include_mask') and not opts.get('mask_filename'):
        print("\n***ERROR***\t--mask_filename is necessary with --include_mask option\n")
        usage(version)
        sys.exit()
    
    if opts.get('output_support') and not opts.get('output_filename'):
        print("\n***ERROR***\t--output_filename is necessary with --output_support option\n")
        usage(version)
        sys.exit()
    
    if opts.get('ucsc') and opts.get('ensembl'):
        print("\n***ERROR***\t--ucsc and --ensembl are mutually exclusive options\n")
        usage(version)
        sys.exit()



get_input(infile_hash, readgroup_hash, mask_hash, mt_hash)
seq_cluster(infile_hash)
link_cluster(infile_hash)
map_cluster(infile_hash, outfile_hash, readgroup_hash)
find_breakpoint(outfile_hash, readgroup_hash, mask_hash)
score_data(outfile_hash)
report(outfile_hash)
