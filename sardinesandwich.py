#!/usr/bin/python
#SardineSandwich 08-2017
#Written by Joseph Sarro
#This script inputs a sam file and a list of coordinates of SSRs and determines how many mapped reads span the SRR, flaning it at least 50bp in each direction.
import sys
#Check to see if reads are single end or aired end.  Tye 2 as an argument for paired end and 1 for single end, default is 1.
SP = 1
if len(sys.argv) >= 2:
  SP = int(sys.argv[1])
print SP
#files to be opened and written.  F needs to be changed to match your SAM file before running as does datafile below.  
f = open('S6b.sam', 'r')
g = open('fullcov.txt','w')
d= open('left5cov.txt','w')
u= open('right3cov.txt','w')
n= open('nocov.txt','w')
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
								print >> g, readlen, "  ",left5span, "  ",right3span
								#print hitcount
							elif int(start) <= (int(hitstart)-50) and int(end) < (int(hitend)+50):
								hitleftcount +=1
								print >> d, readlen, "  ",left5span, "   ",right3span
							elif int(start) > (int(hitstart)-50) and int(end) >= (int(hitend)+50):
								hitrightcount +=1
								print >> u, readlen, "  ",left5span, "   ",right3span
							elif int(start) > (int(hitstart)-50) and int(end) < (int(hitend)+50):
								hitnonecount +=1
								print >> n, readlen, "  ",left5span, "   ",right3span, "        ",field1
							else:
								whatcount +=1
	
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
							print >> g, readlen, "	",left5span, "	",right3span
							#print hitcount
						elif int(start) <= (int(hitstart)-50) and int(end) < (int(hitend)+50):
							hitleftcount +=1
							print >> d, readlen, "	",left5span, "	",right3span
						elif int(start) > (int(hitstart)-50) and int(end) >= (int(hitend)+50):
							hitrightcount +=1
							print >> u, readlen, "	",left5span, "	 ",right3span
						elif int(start) > (int(hitstart)-50) and int(end) < (int(hitend)+50):
							hitnonecount +=1
							print >> n, readlen, "	",left5span, "	",right3span#, "	",field1 
						else:
							whatcount +=1
				#else:
					#print contig
	#print linesi
#print stats
print "Total reads:", count
print "Mapped reads", mapcount
print "Reads mapped to documented contigs", contigcount
print "Reads flanking both sides of an SSR", hitcount
print "Reads that only flank left", hitleftcount
print "Reads that only flank right", hitrightcount
print "Reads that do not flank either side", hitnonecount
print whatcount
#str = "this is string example....wow!!!";

#print "Length of the string: ", len(str)
