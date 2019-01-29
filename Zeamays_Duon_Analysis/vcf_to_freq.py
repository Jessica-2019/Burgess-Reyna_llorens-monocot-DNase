#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 16:03:02 2018

@author: srs62
"""
from __future__ import print_function
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

vcf = sys.argv[1]

out = open("{}.frq".format(vcf), "w")

with open(vcf) as f:
    for line in f:
        line_s = line.split()
        chrom = line_s[0]
        position = line_s[1]
        ref = line_s[2]
        alt = line_s[3].split(",")
        flags = line_s[4].split(";")
        alleles = flags[2].split("=")[1].split(",")
        alleles = [float(i) for i in alleles]
        n_alleles = len(alleles)
        n_chr = int(flags[0].split("=")[1])
        if n_chr != 0:
            allele1 = alleles[0]
            allele1_freq = "{}:{}".format(ref, round(float(allele1/n_chr), 5))
            allele2 = alleles[1]
            allele2_freq = "{}:{}".format(alt[0], round(float(allele2/n_chr), 5))
            freqs = allele1_freq + "\t" + allele2_freq
            try:
                allele3 = alleles[2]
                allele3_freq = "{}:{}".format(alt[1], round(float(allele3/n_chr), 5))
                freqs = freqs + "\t" + allele3_freq
            except:
                pass
            print(chrom, position, n_alleles, n_chr, freqs, sep = "\t", file = out)
