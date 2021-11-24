#!/usr/bin/python
#Fromat Matrix  03-2020
#Written by Joseph Sarro
#Format read 1 matrix file to  
from __future__ import division
import numpy as np
import os
import argparse
from operator import itemgetter
contig_list=["001-QR_002", "002-CD_001", "003-CD_003", "004-CD_005", "005-CD_011", "006-CD_013", "014-CD_028", "007-CD_014", "008-CD_015", "009-CD_017", "010-CD_019", "011-CD_020", "012-CD_024", "013-CD_026", "015-CD_033", "016-CD_042", "017-CD_043", "018-CD_049", "019-CD_053", "020-CD_056", "049-QR_003", "021-CD_057", "022-CD_058", "023-CD_059", "024-CD_061", "025-CD_063", "026-CD_066", "027-CD_067", "029-CD_077", "030-CD_082", "031-CD_084", "032-CD_087", "033-CD_088", "050-CM_001", "051-CM_003", "052-CM_007", "053-CM_014", "075-QA_004", "057-QR_012", "058-QR_041", "059-QR_046", "060-QR_055", "061-QR_056", "074-QA_001", "062-QR_059", "063-QR_062", "064-QR_069", "065-QR_077", "066-QR_088", "067-QR_092", "068-QR_098", "069-QR_122", "070-QR_150", "071-QR_181", "072-QR_186", "073-QR_192", "035-CD_095", "036-CD_100", "037-CD_101", "038-CD_104", "039-CD_105", "047-CD_122", "054-CM_015"]
#contig_list=["001-QR_002", "002-CD_001", "003-CD_003", "004-CD_005", "005-CD_011", "006-CD_013", "007-CD_014", "008-CD_015", "009-CD_017", "010-CD_019", "011-CD_020", "012-CD_024", "013-CD_026", "014-CD_028", "015-CD_033", "016-CD_042", "017-CD_043", "018-CD_049", "019-CD_053", "020-CD_056", "021-CD_057", "022-CD_058", "023-CD_059", "024-CD_061", "025-CD_063", "026-CD_066", "027-CD_067", "029-CD_077", "030-CD_082", "031-CD_084", "032-CD_087", "033-CD_088", "035-CD_095", "036-CD_100", "037-CD_101", "038-CD_104", "039-CD_105", "047-CD_122", "049-QR_003", "050-CM_001", "051-CM_003", "052-CM_007", "053-CM_014", "057-QR_012", "058-QR_041", "059-QR_046", "060-QR_055", "061-QR_056", "062-QR_059", "063-QR_062", "064-QR_069", "065-QR_077", "066-QR_088", "067-QR_092", "068-QR_098", "069-QR_122", "070-QR_150", "071-QR_181", "072-QR_186", "073-QR_192", "074-QA_001", "075-QA_004"]
#sample_list=["2081-D3_S20", "2091-H12_S96", "2094-A6_S41", "2095-G6_S47", "2099-F6_S46", "2100-A4_S25", "2102-F5_S38", "2104-E6_S45", "2110-D6_S44", "2114-E5_S37", "2115-D5_S36", "2118-A7_S49", "2119-G7_S55", "2120-F7_S54", "2122-E7_S53", "2123-D7_S52", "2147-H8_S64", "2150-C1_S3", "2157-D1_S4", "2160-C6_S43", "2161-E2_S13", "2163-E1_S5", "2172-G2_S15", "2174-H2_S16", "2185-H4_S32", "2189-C5_S35", "2191-H6_S48", "2194-B5_S34", "2201-A2_S9", "2209-B2_S10", "2219-F1_S6", "2221-A3_S17", "2233-C2_S11", "2234-D2_S12", "2239-B6_S42", "2244-B3_S18", "2246-C12_S91", "2259-F2_S14", "2260-G1_S7", "2264-H1_S8", "2288-E12_S93", "2293-F12_S94", "2298-D12_S92", "2326-G8_S63", "2332-G9_S71", "2333-F8_S62", "2336-E8_S61", "2361-H7_S56", "2369-B4_S26", "2381-D8_S60", "2385-C8_S59", "2389-B8_S58", "2393-G12_S95", "2394-G5_S39", "2402-C4_S27", "2404-F4_S30", "2408-D4_S28", "2409-E4_S29", "2412-A5_S33", "2416-A8_S57", "2425-H9_S72", "2430-A1_S1", "2435-B1_S2", "2566-G4_S31", "2569-C3_S19", "3119-F9_S70", "3120-E9_S69", "3121-B9_S66", "3123-D9_S68", "3125-C9_S67", "3126-E3_S21", "3127-H3_S24", "3128-F3_S22", "3130-H5_S40", "3135-G3_S23", "3255-A9_S65", "3256-H10_S80", "3257-G10_S79", "3258-F10_S78", "3259-E10_S77", "3260-D10_S76", "3261-C10_S75", "3263-B10_S74", "3264-A10_S73", "3265-H11_S88", "3266-G11_S87", "3267-F11_S86", "3268-E11_S85", "3269-D11_S84", "3270-1-B7_S50", "3270-2-C7_S51", "3271-B12_S90", "3273-C11_S83", "3276-B11_S82", "3302-A11_S81", "ND-CHC-A12_S89"]
#sample_list=["ND-CHC-A1_S1_L001", "ND-CHC-A1_S1_L001"]
sample_list=["CHC3309_S1_L001", "CHC3310_S2_L001", "CHC3311_S3_L001", "CHC3312_S4_L001", "CHC3313_S5_L001", "CHC3314_S6_L001", "CHC3315_S7_L001", "CHC3317_S8_L001", "CHC3318_S9_L001", "CHC3319_S10_L001", "CHC3320_S11_L001", "CHC3321_S12_L001", "CHC3322_S13_L001", "CHC3323_S14_L001", "CHC3325_S15_L001", "CHC3326_S16_L001", "CHC3327_S17_L001", "CHC3328_S18_L001", "CHC3329_S19_L001", "CHC3331_S20_L001", "CHC3332_S21_L001", "CHC3333_S22_L001", "CHC3335_S23_L001", "CHC3336_S24_L001", "CHC3337_S25_L001", "CHC3338_S26_L001", "CHC3339_S27_L001", "CHC3340_S28_L001", "CHC3341_S29_L001", "CHC3342_S30_L001", "CHC3343_S31_L001", "CHC3344_S32_L001", "CHC3345_S33_L001", "CHC3346_S34_L001", "CHC3347_S35_L001", "CHC3349_S36_L001", "CHC3351_S37_L001", "CHC3352_S38_L001", "CHC3356_S39_L001", "CHC3357_S40_L001", "CHC3358_S41_L001", "CHC3359_S42_L001", "CHC3360_S43_L001", "CHC3361_S44_L001", "CHC3364_S45_L001", "CHC3366_S46_L001", "CHC3367_S47_L001", "CHC3368_S48_L001", "CHC3369_S49_L001", "CHC3371_S50_L001", "CHC3372_S51_L001", "CHC3373_S52_L001", "CHC3374_S53_L001", "CHC3375_S54_L001", "CHC3376_S55_L001", "CHC3377_S56_L001", "CHC3379_S57_L001", "CHC3380_S58_L001", "CHC3381_S59_L001", "CHC3382_S60_L001", "CHC3383_S61_L001", "CHC3384_S62_L001", "CHC3386_S63_L001", "CHC3387_S64_L001", "CHC3388_S65_L001", "CHC3389_S66_L001", "CHC3390_S67_L001", "CHC3391_S68_L001", "CHC3392_S69_L001", "CHC3394_S70_L001", "CHC3395_S71_L001", "CHC3396_S72_L001", "CHC3397_S73_L001", "CHC3398_S74_L001", "CHC3399_S75_L001", "CHC3401_S76_L001", "CHC3402_S77_L001", "CHC3403_S78_L001", "CHC3404_S79_L001", "CHC3405_S80_L001", "CHC3406_S81_L001", "CHC3408_S82_L001", "CHC3409_S83_L001", "CHC3410_S84_L001", "CHC3411_S85_L001", "CHC3412_S86_L001", "CHC3413_S87_L001", "CHC3414_S88_L001", "CHC3415_S89_L001", "CHC3416_S90_L001", "CHC3417_S91_L001", "CHC3418_S92_L001", "CHC3419_S93_L001", "CHC3420_S94_L001", "CHC3422_S95_L001", "NegCtrl_S96_L001"]
G= open("VCF_Matrix.txt",'w')
print >> G, "Sample	",
for i in contig_list:
	print >> G, i,"		",  
print >> G
for i in sample_list:
	#os.system("cp "+i+".sam.Matrix.1.txt.trimmed.txt temp")
	#os.system("sed temp -e '1,1d' > tempA")
	#os.system("sed -i -e 's/([0-9]\+)//g' tempA")
	#os.system("rm temp")
	##ta = open("tempA", 'r')
	print >> G, i, 
	for j in contig_list:
		found= False
		prev="NA"
		ta = open(i+".NH.vcf", 'r')
		#print ta
		for lines in ta:
			fields = lines.split("	")
			contig=fields[0]
		#	print contig
			pos=fields[1]
                        alt=fields[4]
                        qual=fields[5]
			if prev.rstrip(" ")==contig.rstrip(" ") and j.rstrip(" ")==contig.rstrip(" "):
				#print "here"
				if float(qual)>=30:
                                        print >> G,",",alt,"(",pos,")",
			elif j.rstrip(" ")==contig.rstrip(" "):
				prev=contig
                        	#ssr = ssr.rstrip("\n")
                       		#ssrs=ssr.split(',')
				#print >> G
				if float(qual)>=30:
					found=True
                                        print >> G,"	",alt,"(",pos,")",
				#else:
				#	print  >> G, "	",ssrs[0],"	",ssrs[1],"	",ssrs[2],
		if found == False:
			print >> G,"	NULL",  	
		ta.close()
	print >> G
#print >> G, "Total_Library","  ",total_count_A,"       ",matched_count_A,"     ",'{:.1%}'.format(matched_count_A/total_count_A)
#os.system("rm tempA")
G.close()
os.system("sed -i -e  's/\t/_NULL_/g' VCF_Matrix.txt")
os.system("sed -i -e  's/\s//g' VCF_Matrix.txt")
os.system("sed -i -e  's/_NULL_/\t/g' VCF_Matrix.txt")
