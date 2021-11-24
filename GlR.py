#!/usr/bin/python
#SnipMatrix  03-2020
#Written by Joseph Sarro
#Generate Length Ratios for read1 and read2 of all samples in a library 
from __future__ import division
import numpy as np
import os
import argparse
from operator import itemgetter
my_list=["1628-C1_S3_L001", "1631-C2_S11_L001", "1632-C3_S19_L001", "1633-C4_S27_L001", "1635-C5_S35_L001", "1660-E1_S5_L001", "1664-E2_S13_L001", "1681-E3_S21_L001", "1698-E4_S29_L001", "1706-E5_S37_L001", "2077-C6_S43_L001", "2078-C7_S51_L001", "2080-H5_S40_L001", "2148-C8_S59_L001", "2153-C9_S67_L001", "2162-C10_S75_L001", "2171-A10_S73_L001", "2173-C11_S83_L001", "2179-E6_S45_L001", "2183-A4_S25_L001", "2187-E7_S53_L001", "2199-C12_S91_L001", "2204-D1_S4_L001", "2211-F8_S62_L001", "2213-E9_S69_L001", "2220-D2_S12_L001", "2225-A5_S33_L001", "2226-A11_S81_L001", "2227-D3_S20_L001", "2228-D4_S28_L001", "2229-A6_S41_L001", "2240-A12_S89_L001", "2243-B1_S2_L001", "2245-B2_S10_L001", "2266-E10_S77_L001", "2289-A7_S49_L001", "2290-D5_S36_L001", "2295-E11_S85_L001", "2297-E12_S93_L001", "2305-F1_S6_L001", "2306-F2_S14_L001", "2307-F3_S22_L001", "2328-F4_S30_L001", "2329-B3_S18_L001", "2330-B4_S26_L001", "2331-B5_S34_L001", "2336-B10_S74_L001", "2352-D6_S44_L001", "2355-D7_S52_L001", "2364-B6_S42_L001", "2375-D8_S60_L001", "2376-D9_S68_L001", "2377-D10_S76_L001", "2380-F5_S38_L001", "2395-F6_S46_L001", "2410-D11_S84_L001", "2417-B7_S50_L001", "2420-D12_S92_L001", "2421-B8_S58_L001", "2605-A8_S57_L001", "2609-F7_S54_L001", "2619-A2_S9_L001", "3125-A3_S17_L001", "3134-E8_S61_L001", "3140-F9_S70_L001", "3156-F10_S78_L001", "3207-H1_S8_L001", "3272-A9_S65_L001", "3274-G1_S7_L001", "3275-G2_S15_L001", "3277-H3_S24_L001", "3278-G4_S31_L001", "3279-G6_S47_L001", "3280-G5_S39_L001", "3281-G7_S55_L001", "3282-G8_S63_L001", "3283-B11_S82_L001", "3284-G9_S71_L001", "3285-G11_S87_L001", "3286-G12_S95_L001", "3288-H2_S16_L001", "3289-G3_S23_L001", "3290-H4_S32_L001", "3291-G10_S79_L001", "3292-H10_S80_L001", "3292-H6_S48_L001", "3293-H7_S56_L001", "3294-H8_S64_L001", "3295-H9_S72_L001", "3297-H11_S88_L001", "3298-H12_S96_L001", "3299-F11_S86_L001", "3300-F12_S94_L001", "3301-B12_S90_L001", "3303-B9_S66_L001", "ND-CHC-A1_S1_L001"]
#my_list=["ND-CHC-A1_S1_L001", "ND-CHC-A1_S1_L001"]
G= open("percent_same.R1.txt",'w')
G2= open("percent_same.R2.txt",'w')
print >> G, "Sample","	","SSR_Total"," ","SSR_Matched","	","Matched_Percent"
print >> G2, "Sample","  ","SSR_Total"," ","SSR_Matched","       ","Matched_Percent"
total_count_A=0
matched_count_A=0
total_count_B=0
matched_count_B=0
for i in my_list:
	contig_total_count=0
	contig_matched_count=0
	os.system("cp "+i+".sam.Matrix.1.txt.trimmed.txt temp")
	os.system("sed temp -e '1,1d' > tempA")
	os.system("sed -i -e 's/([0-9]\+)//g' tempA")
	os.system("cp "+i+".sam.Matrix.2.txt.trimmed.txt temp")
        os.system("sed temp -e '1,1d' > tempB")
        os.system("sed -i -e 's/([0-9]\+)//g' tempB")
	os.system("rm temp")
	ta = open("tempA", 'r')
	#print >> G, field[0],"	",field[1].rstrip("\n")
	for lines in ta:
		#total_count_A += 1
		fields = lines.split("	")
		contig=fields[0]
		ssr=fields[1]
		#print "new line"+contig
#	arr=[]
#	largest=0
#	count=1
#	fill =False
		for ssrs in ssr.split(','):
   			#print ssrs.rstrip("\n") 
			ssrs = ssrs.rstrip("\n")
			#print si
			total_count_A += 1
			contig_total_count +=1
			tb = open("tempB", 'r')
			for lines2 in tb:
				fields2 = lines2.split("	")
                		contig2=fields2[0]
                		ssr2=fields2[1]
				for ssrs2 in ssr2.split(','):
					ssrs2=ssrs2.rstrip("\n")
					#print ssrs
					#print contig
					#print ssrs2
					#print contig2
					if contig==contig2 and ssrs==ssrs2:
						#print ssrs
						matched_count_A +=1
						contig_matched_count +=1
	print >> G, i,"	",contig_total_count,"	",contig_matched_count,"	",'{:.1%}'.format(contig_matched_count/contig_total_count)
	ta.close()
	tb.close()
	contig_total_count=0
        contig_matched_count=0
        tb = open("tempB", 'r')
        #print >> G, field[0]," ",field[1].rstrip("\n")
        for lines in tb:
                #total_count_A += 1
                fields = lines.split("	")
                contig=fields[0]
                ssr=fields[1]
                #print "new line"+contig
#       arr=[]
#       largest=0
#       count=1
#       fill =False
                for ssrs in ssr.split(','):
                        #print ssrs.rstrip("\n") 
                        ssrs = ssrs.rstrip("\n")
                        #print si
                        total_count_B += 1
                        contig_total_count +=1
                        ta = open("tempA", 'r')
                        for lines2 in ta:
                                fields2 = lines2.split("	")
                                contig2=fields2[0]
                                ssr2=fields2[1]
                                for ssrs2 in ssr2.split(','):
                                        ssrs2=ssrs2.rstrip("\n")
                                        #print ssrs
                                        #print contig
                                        #print ssrs2
                                        #print contig2
                                        if contig==contig2 and ssrs==ssrs2:
         #                                       print ssrs
                                                matched_count_B +=1
                                                contig_matched_count +=1
	#					print contig_total_count,"	",contig_matched_count
        print >> G2, i,"	",contig_total_count,"	",contig_matched_count,"	",'{:.1%}'.format(contig_matched_count/contig_total_count)
	ta.close()
	tb.close()
print >> G, "Total_Library","  ",total_count_A,"       ",matched_count_A,"     ",'{:.1%}'.format(matched_count_A/total_count_A)
print >> G2, "Total_Library","	",total_count_B,"	",matched_count_B,"	",'{:.1%}'.format(matched_count_B/total_count_B)
os.system("rm tempA")
os.system("rm tempB")
G.close()
G2.close()
