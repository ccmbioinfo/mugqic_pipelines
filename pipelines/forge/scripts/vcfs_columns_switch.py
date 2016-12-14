#!/bin/env python2

import os, sys, re, string
import datetime, subprocess


# Purpose: To change the INFO column of a VCF file with another while keeping the rest of the fields intact
# Input: Two input file (in this case "variant/sample.0002.vcf" and "alignment/sample.hc.vcf"
# Output: The first VCF file but with the QUAL and INFO fields of the second file


MultiToIndiv_VCF = sys.argv [ 1 ]
SingleSampleVCF = sys.argv [ 2 ]
output_file = sys.argv [ 3 ]

output=open(output_file,"w")

VCFinput=open(SingleSampleVCF,"r")
VCFinput_list=VCFinput.readlines()

for line in VCFinput_list:
    if line.startswith("#"):
	output.write(line)
    else:
	break

last_pos=0


with open(MultiToIndiv_VCF,'r') as MultiToIndiv_input:
    for MTI_line in MultiToIndiv_input:
        if MTI_line.startswith('#'):
            #output.write(MTI_line)
	    continue
	MTI_spline = MTI_line.split("\t")

	for k, SampleVCFline in enumerate(VCFinput_list[last_pos:]):
	    if SampleVCFline.startswith('#'):
	        continue
	    Samplespline = SampleVCFline.split("\t")
	
	    # If the chromosome and position in the two files match; obtain the QUAL and INFO columns from the Sample.hc.vcf file and everything else from the Sample.PASSvqsr.vcf file	   
	    if MTI_spline[0]==Samplespline[0] and MTI_spline[1]==Samplespline[1]:


		SampleInfoFields=Samplespline[7].split(";")
		MTIInfoFields=MTI_spline[7].split(";")


		# We want to keep the INFO field of "SingleSampleVCF" file but AC, AN, AF, MLEAC, and MLEAF that correspond to alles should
		# come from the "MultiToIndiv_VCF" file. AC, AN, and AF appear in the beginnin of the INFO files so IN THE "HybridInfo" below
		# we get the first three fields.
		# For MLEAC and MLEAF, we first have to obtain the correct ("corr") below and then replace them with the "wrong" ones in the 
		# HybridInfoFields list. The list is then joined together with ";" character to make up the regular INFO field format
		MLEAC_corr_value=''
		MLEAF_corr_value=''
		for el in MTIInfoFields:
			if "MLEAC" in el: MLEAC_corr_value=el
			elif "MLEAF" in el: MLEAF_corr_value=el

		HybridInfo=MTIInfoFields[0]+";"+MTIInfoFields[1]+";"+MTIInfoFields[2]+";"+";".join(SampleInfoFields[3:])
		HybridInfoFields = HybridInfo.split(";")

		MLEAC_wrong_value=''
		MLEAF_wrong_value=''
		for element in SampleInfoFields:
			if "MLEAC" in element: MLEAC_wrong_value=element
			elif "MLEAF" in element: MLEAF_wrong_value=element
	
		# replace the old/wrong MLEAC and MLEAF with new/correct ones
		HybridInfoFields = [w.replace(MLEAC_wrong_value, MLEAC_corr_value) for w in HybridInfoFields]
		HybridInfoFields = [w.replace(MLEAF_wrong_value, MLEAF_corr_value) for w in HybridInfoFields]

		# convert the list to a string	
		HybridInfo=";".join(HybridInfoFields)

	        output.write(MTI_spline[0]+"\t"+MTI_spline[1]+"\t"+MTI_spline[2]+"\t"+MTI_spline[3]+"\t"+MTI_spline[4]+"\t"+Samplespline[5]+"\t"+MTI_spline[6]+"\t"+HybridInfo+"\t"+MTI_spline[8]+"\t"+MTI_spline[9])
	        last_pos+=k
		break

output.close()
VCFinput.close()
