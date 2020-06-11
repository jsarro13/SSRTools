#!/usr/bin/python
#Fromat Matrix  03-2020
#Written by Joseph Sarro
#Format read 1 matrix file to  
from __future__ import division
import numpy as np
import os
import argparse
from operator import itemgetter
contig_list=["001-QR_002", "002-CD_001", "003-CD_003", "004-CD_005", "005-CD_011", "006-CD_013", "007-CD_014", "008-CD_015", "009-CD_017", "010-CD_019", "011-CD_020", "012-CD_024", "013-CD_026", "014-CD_028", "015-CD_033", "016-CD_042", "017-CD_043", "018-CD_049", "019-CD_053", "020-CD_056", "021-CD_057", "022-CD_058", "023-CD_059", "024-CD_061", "025-CD_063", "026-CD_066", "027-CD_067", "029-CD_077", "030-CD_082", "031-CD_084", "032-CD_087", "033-CD_088", "035-CD_095", "036-CD_100", "037-CD_101", "038-CD_104", "039-CD_105", "047-CD_122", "049-QR_003", "050-CM_001", "051-CM_003", "052-CM_007", "053-CM_014", "057-QR_012", "058-QR_041", "059-QR_046", "060-QR_055", "061-QR_056", "062-QR_059", "063-QR_062", "064-QR_069", "065-QR_077", "066-QR_088", "067-QR_092", "068-QR_098", "069-QR_122", "070-QR_150", "071-QR_181", "072-QR_186", "073-QR_192", "074-QA_001", "075-QA_004"]
sample_list=["1628-C1_S3_L001", "1631-C2_S11_L001", "1632-C3_S19_L001", "1633-C4_S27_L001", "1635-C5_S35_L001", "1660-E1_S5_L001", "1664-E2_S13_L001", "1681-E3_S21_L001", "1698-E4_S29_L001", "1706-E5_S37_L001", "2077-C6_S43_L001", "2078-C7_S51_L001", "2080-H5_S40_L001", "2148-C8_S59_L001", "2153-C9_S67_L001", "2162-C10_S75_L001", "2171-A10_S73_L001", "2173-C11_S83_L001", "2179-E6_S45_L001", "2183-A4_S25_L001", "2187-E7_S53_L001", "2199-C12_S91_L001", "2204-D1_S4_L001", "2211-F8_S62_L001", "2213-E9_S69_L001", "2220-D2_S12_L001", "2225-A5_S33_L001", "2226-A11_S81_L001", "2227-D3_S20_L001", "2228-D4_S28_L001", "2229-A6_S41_L001", "2240-A12_S89_L001", "2243-B1_S2_L001", "2245-B2_S10_L001", "2266-E10_S77_L001", "2289-A7_S49_L001", "2290-D5_S36_L001", "2295-E11_S85_L001", "2297-E12_S93_L001", "2305-F1_S6_L001", "2306-F2_S14_L001", "2307-F3_S22_L001", "2328-F4_S30_L001", "2329-B3_S18_L001", "2330-B4_S26_L001", "2331-B5_S34_L001", "2336-B10_S74_L001", "2352-D6_S44_L001", "2355-D7_S52_L001", "2364-B6_S42_L001", "2375-D8_S60_L001", "2376-D9_S68_L001", "2377-D10_S76_L001", "2380-F5_S38_L001", "2395-F6_S46_L001", "2410-D11_S84_L001", "2417-B7_S50_L001", "2420-D12_S92_L001", "2421-B8_S58_L001", "2605-A8_S57_L001", "2609-F7_S54_L001", "2619-A2_S9_L001", "3125-A3_S17_L001", "3134-E8_S61_L001", "3140-F9_S70_L001", "3156-F10_S78_L001", "3207-H1_S8_L001", "3272-A9_S65_L001", "3274-G1_S7_L001", "3275-G2_S15_L001", "3277-H3_S24_L001", "3278-G4_S31_L001", "3279-G6_S47_L001", "3280-G5_S39_L001", "3281-G7_S55_L001", "3282-G8_S63_L001", "3283-B11_S82_L001", "3284-G9_S71_L001", "3285-G11_S87_L001", "3286-G12_S95_L001", "3288-H2_S16_L001", "3289-G3_S23_L001", "3290-H4_S32_L001", "3291-G10_S79_L001", "3292-H10_S80_L001", "3292-H6_S48_L001", "3293-H7_S56_L001", "3294-H8_S64_L001", "3295-H9_S72_L001", "3297-H11_S88_L001", "3298-H12_S96_L001", "3299-F11_S86_L001", "3300-F12_S94_L001", "3301-B12_S90_L001", "3303-B9_S66_L001", "ND-CHC-A1_S1_L001"]
#sample_list=["ND-CHC-A1_S1_L001", "ND-CHC-A1_S1_L001"]
G= open("SNP_Matrix.txt",'w')
print >> G, "Sample	",
for i in contig_list:
	print >> G, i,"		",  
print >> G
for i in sample_list:
	os.system("cp "+i+".sam.Matrix.1.txt.trimmed.txt temp")
	os.system("sed temp -e '1,1d' > tempA")
	os.system("sed -i -e 's/([0-9]\+)//g' tempA")
	os.system("rm temp")
	#ta = open("tempA", 'r')
	print >> G, i, 
	for j in contig_list:
		found= False
		ta = open("tempA", 'r')
		for lines in ta:
			fields = lines.split("	")
			contig=fields[0]
			if j.rstrip(" ")==contig.rstrip(" "):
				found=True
				ssr=fields[1]
                        	ssr = ssr.rstrip("\n")
                       		ssrs=ssr.split(',')
				#ssr_list=ssr.split(',')
			#	print contig
				#print ssr
				#print ssrs
				#print len(ssrs)
				#print ssrs[0]
				if len(ssrs)==1:
					print >> G, "	", ssrs[0],"	",ssrs[0],
				elif len(ssrs)>1:
					print >> G, "	",ssrs[0],"	",ssrs[1],
				#else:
				#	print  >> G, "	",ssrs[0],"	",ssrs[1],"	",ssrs[2],
		if found == False:
			print >> G,"	NULL	NULL",  	
		ta.close()
	print >> G
#print >> G, "Total_Library","  ",total_count_A,"       ",matched_count_A,"     ",'{:.1%}'.format(matched_count_A/total_count_A)
os.system("rm tempA")
G.close()
