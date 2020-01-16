#!/usr/bin/python
#SamIAm 02-2018
#Written by Joseph Sarro
#This script inputs a sam file and a list of coordinates of SSRs and prints reads that span the SRR 50 bp in each direction.
import sys
import argparse
#Check to see if reads are single end or paired end.  Tye 2 as an argument for paired end and 1 for single end, default is 1.
SP = 1
#if len(sys.argv) >= 2:
#SP = int(sys.argv[1])
#print SP
#files to be opened and written.  F needs to be changed to match your SAM file before running as does datafile below.  
samplename =' '
parser = argparse.ArgumentParser(description='Filters a SAM file with reads that span an SSR at least 50bp in each direction.')
parser.add_argument('sam', type=str,help='Enter a SAM file')
parser.add_argument('-p', dest='pairs',type=str,default='no',help='Enter yes if these are paired reads (Defualt=no)')
args = parser.parse_args()
if args.pairs == 'yes':
	SP=2
print SP
samplename= args.sam
f = open(samplename, 'r')
#f = open('S6b.sam', 'r')
g = open(samplename+'.filter50.sam','w')
count=0
mapcount=0
contigcount=0
hitcount=0
hitrightcount=0
hitleftcount=0
hitnonecount=0
whatcount=0
#Are the reads aired for single end?
if SP ==2:
	#print "hola"
	linecount=1
	readlen=0
	start=0
	for lines in f:
		#count +=1
		fields = lines.split("	")
		sq=fields[0]
		#print sq
		#ignore lines with headers
		if sq != "@SQ" and sq != "@RG" and sq != "@PG":
			linecount+=1
			if linecount % 2 != 0:
				#print "what"
				count +=1
				start=fields[3]
				field1= fields[9]
				readlen=len(field1)
				#print readlen
			else:
				#count +=1
                                endstart=fields[3]
                                field1= fields[9]
                                endreadlen=len(field1)
				end = int(endstart)+endreadlen-1
				contig=fields[2]
				#all mapped reads
				if contig !='*':
					mapcount+=1
					#A file containing adresses of all SSRs, datafile needs to be changed to your file name.  Will mark as a hit if the read spands 50 bp to the left and right of the SRR.  
					#Files containing the length of the sequence, number of bp spanning left of the SSR, and right of the SSR, will be printed to appropriate output files. 
					datafile = file('Cumulative_Bait_Info.txt')
					for line in datafile:
						#print line                     
						if contig in line:
							#print readlen
							#print contig
							contigcount+=1
							hit=line.split("	")
							hitstart=hit[1]
							hitend=hit[2]
							#sub=int(hitstart)-50
							#add=int(hitend)+50
							#print start,  "  ", sub, " ",end, " ",add
							left5span= int(hitstart)-int(start)
							right3span=int(end)-int(hitend)
							if int(start) <= (int(hitstart)-50) and int(end) >= (int(hitend)+50):
								hitcount +=1
								print >> g, lines,
		elif sq == "@SQ" and sq == "@RG" and sq == "@PG":
			print >> g, lines,	
else:
#for all lines in the sam file
	for lines in f:
		#count +=1
		fields = lines.split("	")
		sq=fields[0]
		#print sq
		#ignore lines with headers
		if sq != "@SQ" and sq != "@RG" and sq != "@PG":
			#print "what"
			count +=1
			start=fields[3]
			field1= fields[9]
			readlen=len(field1)
			#print readlen
			end = int(start)+readlen-1
			contig=fields[2]
			#all mapped reads
			if contig !='*':
				mapcount+=1
				#A file containing adresses of all SSRs, datafile needs to be changed to your file name.  Will mark as a hit if the read spands 50 bp to the left and right of the SRR.  
				#Files containing the length of the sequence, number of bp spanning left of the SSR, and right of the SSR, will be printed to appropriate output files. 
				datafile = file('Cumulative_Bait_Info.txt')
				for line in datafile:
					#print line			
					if contig in line:
						#print readlen
						#print contig
						contigcount+=1
						hit=line.split("	")
						hitstart=hit[1]
						hitend=hit[2]
						#sub=int(hitstart)-50
						#add=int(hitend)+50
						#print start,  "  ", sub, " ",end, " ",add
						left5span= int(hitstart)-int(start)
						right3span=int(end)-int(hitend)
						if int(start) <= (int(hitstart)-50) and int(end) >= (int(hitend)+50):
							hitcount +=1
							print >> g, lines,
		else:
			print >> g, lines,
