"""
New York Genome Center
SOFTWARE COPYRIGHT NOTICE AGREEMENT
This software and its documentation are copyright (2016) by the New York
Genome Center. All rights are reserved. This software is supplied without
any warranty or guaranteed support whatsoever. The New York Genome Center
cannot be responsible for its use, misuse, or functionality.
Version: 1.0
Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)
"""
import sys
import os
import argparse
import imp
import pkg_resources
from collections import defaultdict
import numpy as np
from math import pow
from conpair import contamination_model

HOMOZYGOUS_P_VALUE_THRESHOLD = 0.999


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


def get_parser():
    desc = 'Program to estimate tumor-normal sample contaminationz'
    parser = argparse.ArgumentParser(description=desc)
    required = parser.add_argument_group('required arguments')
    required.add_argument('-T', '--tumor_pileup', help='TUMOR PILEUP FILE', required=True)
    required.add_argument('-N', '--normal_pileup', help='NORMAL PILEUP FILE', required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-M', '--markers', help='MARKER FILE [default: markers for GRCh37]')
    optional.add_argument('-O', '--outfile', help='TXT OUTPUT FILE [%(default)s]', default="-")
    optional.add_argument('-G', '--grid', help='GRID INTERVAL [%(default)s]', default=0.01, type=float)
    optional.add_argument('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY [%(default)s]', default=10, type=int)
    return parser


def main(args=None):
    parser = get_parser()
    args = parser.parse_args(args)

    if not os.path.exists(args.tumor_pileup):
        err = 'ERROR: Input tumor file {0} cannot be find.'.format(args.tumor_pileup)
        parser.error(err)
    
    if not os.path.exists(args.normal_pileup):
        err = 'ERROR: Input normal file {0} cannot be find.'.format(args.normal_pileup)
        parser.error(err)
        
    if args.markers:
        MARKER_FILE = args.markers
    else:
        path = os.path.join('data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')
        MARKER_FILE = pkg_resources.resource_filename('conpair', path)

    if not os.path.exists(MARKER_FILE):
        err = 'ERROR: Marker file {0} cannot be find.'.format(MARKER_FILE)
        parser.error(err)


    grid_precision = args.grid
    MMQ = args.min_mapping_quality


    Markers = contamination_marker.get_markers(MARKER_FILE)

    Normal_homozygous_genotype = defaultdict(lambda: defaultdict())

    checkpoints = [i for i in drange(0.0, 1.0, grid_precision)]
    checkpoints.append(0.5)
    Scores = contamination_model.create_conditional_likelihood_of_base_dict(checkpoints)
    checkpoints = [i for i in drange(0.0, 0.5, grid_precision)]
    checkpoints.append(0.5)

    if args.outfile != "-":
        outfile = open(args.outfile, 'w')
        
    ### PARSING THE NORMAL PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

    file = open(args.normal_pileup)
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = contamination_marker.parse_mpileup_line(line, min_map_quality=MMQ)
        try:
            marker = Markers[pileup.chrom + ":" + pileup.pos]
        except:
            continue
        
        if pileup.Quals[marker.ref] == [] and pileup.Quals[marker.alt] == []:
            continue

        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]
        
        AA_likelihood, AB_likelihood, BB_likelihood = genotypes.compute_genotype_likelihood(pileup.Quals[marker.ref], pileup.Quals[marker.alt], normalize=True)

        if AA_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {'genotype': marker.ref, 'AA_likelihood': AA_likelihood, 'AB_likelihood': AB_likelihood, 'BB_likelihood': BB_likelihood}
        elif BB_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {'genotype': marker.alt, 'AA_likelihood': AA_likelihood, 'AB_likelihood': AB_likelihood, 'BB_likelihood': BB_likelihood}
            
        
        p_AA, p_AB, p_BB = genotypes.RAF2genotypeProb(RAF)
        lPAA = math_operations.log10p(p_AA)
        lPAB = math_operations.log10p(p_AB)
        lPBB = math_operations.log10p(p_BB)

        
        priors = [lPAA*2, lPAA+lPBB, lPAA+lPAB, lPAB*2, lPAB+lPAA, lPAB+lPBB, lPBB*2, lPBB+lPAA, lPBB+lPAB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)
        
    file.close()
    D = contamination_model.calculate_contamination_likelihood(checkpoints, Data, Scores)
    ARGMAX = np.argmax(D)
    cont = checkpoints[ARGMAX]

    x1 = max(cont-grid_precision, 0.0)
    x2 = cont
    x3 = min(cont+grid_precision, 1.0)

    if x2 == 0.0:
        x2 += grid_precision/100
    elif x2 == 1.0:
        x2 -= grid_precision/100

    ### SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

    optimal_val = contamination_model.apply_brents_algorithm(Data, Scores, x1, x2, x3)

    ### PRINTING THE NORMAL RESULTS
        
    if args.outfile == "-":
        print("Normal sample contamination level: " + str(round(100.0*optimal_val, 3)) + "%")
    else:
        outfile.write("Normal sample contamination level: " + str(round(100.0*optimal_val, 3)) + "%\n")


    ### PARSING THE TUMOR PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

    file = open(args.tumor_pileup)

    checkpoints = [i for i in drange(0.0, 1.0, grid_precision)]
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"): 
            continue
        pileup = contamination_marker.parse_mpileup_line(line, min_map_quality=MMQ)
        
        try:
            normal_hom_genotype = Normal_homozygous_genotype[pileup.chrom][pileup.pos]['genotype']
        except:
            continue

        try:
            marker = Markers[pileup.chrom + ":" + pileup.pos]
        except:
            continue
        
        if pileup.Quals[marker.ref] == [] and pileup.Quals[marker.alt] == []:
            continue
        
        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]
        
        Normal_info = Normal_homozygous_genotype[pileup.chrom][pileup.pos]
        AA_likelihood = Normal_info['AA_likelihood']
        AB_likelihood = Normal_info['AB_likelihood']
        BB_likelihood = Normal_info['BB_likelihood']
        nlPAA = math_operations.log10p(AA_likelihood)
        nlPAB = math_operations.log10p(AB_likelihood)
        nlPBB = math_operations.log10p(BB_likelihood)
        
        
        p_AA, p_AB, p_BB = genotypes.RAF2genotypeProb(RAF)
        lPAA = math_operations.log10p(p_AA)
        lPAB = math_operations.log10p(p_AB)
        lPBB = math_operations.log10p(p_BB)
        priors = [lPAA+nlPAA, lPBB+nlPAA,lPAB+nlPAA, lPAB+nlPAB, lPAA+nlPAB, lPBB+nlPAB, lPBB+nlPBB, lPAA+nlPBB, lPAB+nlPBB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)
        
    file.close()

    D = contamination_model.calculate_contamination_likelihood(checkpoints, Data, Scores)
    ARGMAX = np.argmax(D)
    cont = checkpoints[ARGMAX]

    x1 = max(cont-grid_precision, 0.0)
    x2 = cont
    x3 = min(cont+grid_precision, 1.0)

    if x2 == 0.0:
        x2 += grid_precision/100
    elif x2 == 1.0:
        x2 -= grid_precision/100

    ### SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

    optimal_val = contamination_model.apply_brents_algorithm(Data, Scores, x1, x2, x3)

    ### PRINTING THE TUMOR RESULTS
        
    if args.outfile == "-":
        print("Tumor sample contamination level: " + str(round(100.0*optimal_val,3)) + "%")
    else:
        outfile.write("Tumor sample contamination level: " + str(round(100.0*optimal_val,3)) + "%\n")
        outfile.close()


if __name__ == '__main__':
    main()
