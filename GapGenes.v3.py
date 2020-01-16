#!/usr/bin/python
#GapGenes 2 02-2018
#Written by Joseph Sarro
#This script measures the length of a specified SSR in a read mapped to a contig. It prints the SSR, repeat number, and both flanking ends. A second file contains count information of SSR lengths. 
import os
import argparse
#note: Enter a SAM file in the command line and make sure the bait info file matches the one below
samplename =' '
parser = argparse.ArgumentParser(description='Find SSR lengths and print sequence information.')
parser.add_argument('-sam', required=True, type=str,help='Enter a SAM file')
parser.add_argument('-C', type=str, required=True, help='Enter a coordinates file')
parser.add_argument('-P', type=str, default="unpaired",
                   help='paired or unpaired')
args = parser.parse_args()
samplename= args.sam
f = open(samplename, 'r')
if args.P == "unpaired":
	G= open(samplename+'.GapResults.txt','w')
	#Read Sam file and see if the read aligned to a contig with a known SSR
	for lines in f:
		#count +=1
		fields = lines.split("	")
		sq=fields[0]
		endlen=0
		#print sq
		#ignore lines with headers
		if sq != "@SQ" and sq != "@RG" and sq != "@PG":
			#print "what"
			#start=fields[3]
			field1= fields[9]
			#readlen=len(field1)
			#print readlen
			#end = int(start)+readlen-1
			contig=fields[2]
			flag=fields[1]
			#if int(flag) > 192 and int(flag) < 2112:
			#	print int(flag)
			#all mapped reads
			if contig !='*':# and int(flag) < 192 :
				#A file containing adresses of all SSRs, datafile needs to be changed to your file name.  Will mark as a hit if the read spands 50 bp to the left and right of the SRR.  
				#Get sequence 10bo left and right of the SSR in the contig and search for those sequenes in the read.  If the read is found print the contig and Gap in the read 
				datafile = file(args.C)
				#print contig
				for line in datafile:
					#print line                     
					if contig in line:
						fields2 = line.split("	")
						SSRstart=fields2[1]
						SSRend=fields2[2]
						Seq=fields2[3]
						Leftstart=int(SSRstart)-11
						Rightend=int(SSRend)+10
						endlen=len(field1)
						Left=Seq[int(Leftstart):int(SSRstart)-1] 
						Right=Seq[int(SSRend):int(Rightend)] 
						#Left= [0:int(SSRstart)]
						#Right= Seq[int(SSRend):int(endlen-1)]
						#print Left
						#print "mid"
						#print Right
						#print contig
						#print Right[1:10]
						#print Right[9]
						#print field1.find(Left)
						#print field1.find(Right)
						#search for sequence and mke sure the right sequence is at lest 10 bases away from the start of the left
						FoundLeft=int(field1.find(Left))
						FoundRight=int(field1.find(Right,FoundLeft+10))
						#print FoundLeft
						#print FoundRight
						if FoundRight >9 and FoundLeft >=0:
							#print "BANG"
							Gap=FoundRight-(FoundLeft+10)
							WholeSequence=field1[int(FoundLeft):(int(FoundRight)+10)]
							#print >> G, contig,"	",Gap,"	",WholeSequence
							#string array ends at last number -1
							SSR=field1[int(FoundLeft+10):int(FoundRight)]
							LeftSeq=field1[0:int(FoundLeft+10)]
							RightSeq=field1[int(FoundRight):int(endlen)]
							#print field1[0:2]
							#print field1[int(endlen-1)]
							#print endlen
							print >> G, contig,"	",Gap,"	",LeftSeq,"	",SSR,"	",RightSeq
#####################################
#End SSR search. Begin file formating
#####################################
	#Sort the output file by contig
	#os.system("echo "+samplename+".BAM")
	G.close()
	f.close()
	#sort the receeding output and then generate a new file that counts and prints number of occurences for each gap in a contig
	#This Matrix file can then be joined with other matrix files from other samples using the cutstom join shell script
	os.system("sort "+samplename+".GapResults.txt > "+samplename+".GapResults.sorted.txt")
	import sys
	Print = sys.stdout.write
	S = open(samplename+'.GapResults.sorted.txt', 'r')
	P = open(samplename+'.Matrix.txt','w')
	curnumcount=0
	prevcontig=""
	prevnum=0
	#print header
	print >> P, "Contig","        ",samplename
	for lineS in S:
		curline=lineS.split("	")
		#print curline
		curcontig=curline[0]
		curnum=curline[1]
		#print curnum
		#print prevcontig,curcontig,prevnum,curnum,curnumcount
		if curcontig == prevcontig and int(curnum) == int(prevnum):
			#if the next line is the same contig and SSR, increment the count by 1
			curnumcount+=1
		elif curcontig == prevcontig and int(curnum) != int(prevnum):
			#if the next line is the same contig but a different SSR, print the previsou SSR info and begin storing info for a new one
			print >> P, prevnum,"(",curnumcount,")",
			prevnum=curnum
			curnumcount=1
		else:
			#this is where it starts, at new contig with no SSR info, this line will be deleted in a sed command below.
			#future occurnecs of this result in a new contig
			print >> P, prevnum, "(",curnumcount,")"
			prevcontig=curcontig
			prevnum=curnum
			curnumcount=1
			print >> P, ""
			print >> P, prevcontig,"	",
	#print final SSR info
	print >> P, prevnum, "(",curnumcount,")"
	P.close()
	S.close()
	#sed commands to format the matrix
	#remove spaces from print at beggining of number
	os.system("sed -i -e  's/\s\s(\s/(/g' "+samplename+".Matrix.txt")
	#add comma to multiple entires in same contig
	os.system("sed -i -e  's/\s)\s/),/g' "+samplename+".Matrix.txt")
	#remove space at end of number
	os.system("sed -i -e  's/\s)/)/g' "+samplename+".Matrix.txt")
	#remove first line containing len 0 with count 0
	os.system("sed -i -e  '/^$\|0 ( 0)/d' "+samplename+".Matrix.txt")
	os.system("rm "+samplename+".GapResults.txt")
	#os.system("sed -i -e  's/,\n/\r/g' Matrix.txt")
###########################################################################
#End of single end option. Begin paired end option
###########################################################################
elif args.P == "paired": 
	#print "it's paired"

	G1= open(samplename+'.GapResults.1.txt','w')
	G2= open(samplename+'.GapResults.2.txt','w')
	#Read Sam file and see if the read aligned to a contig with a known SSR
	j = 1
	for lines in f:
		#count +=1
		fields = lines.split("	")
		sq=fields[0]
		endlen=0
		#print sq
		#ignore lines with headers
		if sq != "@SQ" and sq != "@RG" and sq != "@PG":
			#print "what"
			#start=fields[3]
			field1= fields[9]
			#readlen=len(field1)
			#print readlen
			#end = int(start)+readlen-1
			contig=fields[2]
			flag=fields[1]
			#all mapped reads
			if contig !='*':
				#if j % 2 != 0 :
				if (int(flag) < 128) or ((int(flag)> 2111) and (int(flag) < 2176)):
#########################
#Read 1
#########################
					#A file containing adresses of all SSRs, datafile needs to be changed to your file name.  Will mark as a hit if the read spands 50 bp to the left and right of the SRR.  
					#Get sequence 10bo left and right of the SSR in the contig and search for those sequenes in the read.  If the read is found print the contig and Gap in the read 
					datafile = file(args.C)
					#print contig
					for line in datafile:
						#print line                     
						if contig in line:
							fields2 = line.split("	")
							SSRstart=fields2[1]
							SSRend=fields2[2]
							Seq=fields2[3]
							Leftstart=int(SSRstart)-11
							Rightend=int(SSRend)+10
							endlen=len(field1)
							Left=Seq[int(Leftstart):int(SSRstart)-1] 
							Right=Seq[int(SSRend):int(Rightend)] 
							#Left= [0:int(SSRstart)]
							#Right= Seq[int(SSRend):int(endlen-1)]
							#print Left
							#print "mid"
							#print Right
							#print contig
							#print Right[1:10]
							#print Right[9]
							#print field1.find(Left)
							#print field1.find(Right)
							#search for sequence and mke sure the right sequence is at lest 10 bases away from the start of the left
							FoundLeft=int(field1.find(Left))
							FoundRight=int(field1.find(Right,FoundLeft+10))
							#print FoundLeft #####
							#print FoundRight ######
							if FoundRight >9 and FoundLeft >=0:
								#print "BANG"
								Gap=FoundRight-(FoundLeft+10)
								WholeSequence=field1[int(FoundLeft):(int(FoundRight)+10)]
								#print >> G, contig,"	",Gap,"	",WholeSequence
								#string array ends at last number -1
								SSR=field1[int(FoundLeft+10):int(FoundRight)]
								LeftSeq=field1[0:int(FoundLeft+10)]
								RightSeq=field1[int(FoundRight):int(endlen)]
								#print field1[0:2]
								#print field1[int(endlen-1)]
								#print endlen
								print >> G1, contig,"	",Gap,"	",LeftSeq,"	",SSR,"	",RightSeq
				#	j += 1 #increment read1
				
				#elif j%2 == 0 :
				elif (int(flag) < 192) or (int(flag) > 2175):
##########################
#Read 2
##########################
					#A file containing adresses of all SSRs, datafile needs to be changed to your file name.  Will mark as a hit if the read spands 50 bp to the left and right of the SRR.  
					#Get sequence 10bo left and right of the SSR in the contig and search for those sequenes in the read.  If the read is found print the contig and Gap in the read 
					datafile = file(args.C)
					#print contig
					for line in datafile:
						#print line                     
						if contig in line:
							fields2 = line.split("	")
							SSRstart=fields2[1]
							SSRend=fields2[2]
							Seq=fields2[3]
							Leftstart=int(SSRstart)-11
							Rightend=int(SSRend)+10
							endlen=len(field1)
							Left=Seq[int(Leftstart):int(SSRstart)-1] 
							Right=Seq[int(SSRend):int(Rightend)] 
							#Left= [0:int(SSRstart)]
							#Right= Seq[int(SSRend):int(endlen-1)]
							#print Left
							#print "mid"
							#print Right
							#print contig
							#print Right[1:10]
							#print Right[9]
							#print field1.find(Left)
							#print field1.find(Right)
							#search for sequence and mke sure the right sequence is at lest 10 bases away from the start of the left
							FoundLeft=int(field1.find(Left))
							FoundRight=int(field1.find(Right,FoundLeft+10))
							#print FoundLeft
							#print FoundRight
							if FoundRight >9 and FoundLeft >=0:
								#print "BANG"
								Gap=FoundRight-(FoundLeft+10)
								WholeSequence=field1[int(FoundLeft):(int(FoundRight)+10)]
								#print >> G, contig,"	",Gap,"	",WholeSequence
								#string array ends at last number -1
								SSR=field1[int(FoundLeft+10):int(FoundRight)]
								LeftSeq=field1[0:int(FoundLeft+10)]
								RightSeq=field1[int(FoundRight):int(endlen)]
								#print field1[0:2]
								#print field1[int(endlen-1)]
								#print endlen
								print >> G2, contig,"	",Gap,"	",LeftSeq,"	",SSR,"	",RightSeq
				j += 1 #increment read2
#####################################
#End SSR search. Begin file formating
#####################################
	#Sort the output file by contig
	#os.system("echo "+samplename+".BAM")
	G1.close()
	G2.close()
	f.close()
############################
#begin formatting read1
############################
	#sort the receeding output and then generate a new file that counts and prints number of occurences for each gap in a contig
	#This Matrix file can then be joined with other matrix files from other samples using the cutstom join shell script
	os.system("sort "+samplename+".GapResults.1.txt > "+samplename+".GapResults.sorted.1.txt")
	import sys
	Print = sys.stdout.write
	S = open(samplename+'.GapResults.sorted.1.txt', 'r')
	P = open(samplename+'.Matrix.1.txt','w')
	curnumcount=0
	prevcontig=""
	prevnum=0
	#print header
	print >> P, "Contig","        ",samplename
	for lineS in S:
		curline=lineS.split("	")
		#print curline
		curcontig=curline[0]
		curnum=curline[1]
		#print curnum
		#print prevcontig,curcontig,prevnum,curnum,curnumcount
		if curcontig == prevcontig and int(curnum) == int(prevnum):
			#if the next line is the same contig and SSR, increment the count by 1
			curnumcount+=1
		elif curcontig == prevcontig and int(curnum) != int(prevnum):
			#if the next line is the same contig but a different SSR, print the previsou SSR info and begin storing info for a new one
			print >> P, prevnum,"(",curnumcount,")",
			prevnum=curnum
			curnumcount=1
		else:
			#this is where it starts, at new contig with no SSR info, this line will be deleted in a sed command below.
			#future occurnecs of this result in a new contig
			print >> P, prevnum, "(",curnumcount,")"
			prevcontig=curcontig
			prevnum=curnum
			curnumcount=1
			print >> P, ""
			print >> P, prevcontig,"	",
	#print final SSR info
	print >> P, prevnum, "(",curnumcount,")"
	P.close()
	S.close()
	#sed commands to format the matrix
	#remove spaces from print at beggining of number
	os.system("sed -i -e  's/\s\s(\s/(/g' "+samplename+".Matrix.1.txt")
	#add comma to multiple entires in same contig
	os.system("sed -i -e  's/\s)\s/),/g' "+samplename+".Matrix.1.txt")
	#remove space at end of number
	os.system("sed -i -e  's/\s)/)/g' "+samplename+".Matrix.1.txt")
	#remove first line containing len 0 with count 0
	os.system("sed -i -e  '/^$\|0 ( 0)/d' "+samplename+".Matrix.1.txt")
	os.system("rm "+samplename+".GapResults.1.txt")
	#os.system("sed -i -e  's/,\n/\r/g' Matrix.txt")
############################
#begin formatting read2
############################

	#sort the receeding output and then generate a new file that counts and prints number of occurences for each gap in a contig
	#This Matrix file can then be joined with other matrix files from other samples using the cutstom join shell script
	os.system("sort "+samplename+".GapResults.2.txt > "+samplename+".GapResults.sorted.2.txt")
	import sys
	Print = sys.stdout.write
	S = open(samplename+'.GapResults.sorted.2.txt', 'r')
	P = open(samplename+'.Matrix.2.txt','w')
	curnumcount=0
	prevcontig=""
	prevnum=0
	#print header
	print >> P, "Contig","        ",samplename
	for lineS in S:
		curline=lineS.split("	")
		#print curline
		curcontig=curline[0]
		curnum=curline[1]
		#print curnum
		#print prevcontig,curcontig,prevnum,curnum,curnumcount
		if curcontig == prevcontig and int(curnum) == int(prevnum):
			#if the next line is the same contig and SSR, increment the count by 1
			curnumcount+=1
		elif curcontig == prevcontig and int(curnum) != int(prevnum):
			#if the next line is the same contig but a different SSR, print the previsou SSR info and begin storing info for a new one
			print >> P, prevnum,"(",curnumcount,")",
			prevnum=curnum
			curnumcount=1
		else:
			#this is where it starts, at new contig with no SSR info, this line will be deleted in a sed command below.
			#future occurnecs of this result in a new contig
			print >> P, prevnum, "(",curnumcount,")"
			prevcontig=curcontig
			prevnum=curnum
			curnumcount=1
			print >> P, ""
			print >> P, prevcontig,"	",
	#print final SSR info
	print >> P, prevnum, "(",curnumcount,")"
	P.close()
	S.close()
	#sed commands to format the matrix
	#remove spaces from print at beggining of number
	os.system("sed -i -e  's/\s\s(\s/(/g' "+samplename+".Matrix.2.txt")
	#add comma to multiple entires in same contig
	os.system("sed -i -e  's/\s)\s/),/g' "+samplename+".Matrix.2.txt")
	#remove space at end of number
	os.system("sed -i -e  's/\s)/)/g' "+samplename+".Matrix.2.txt")
	#remove first line containing len 0 with count 0
	os.system("sed -i -e  '/^$\|0 ( 0)/d' "+samplename+".Matrix.2.txt")
	os.system("rm "+samplename+".GapResults.2.txt")
	#os.system("sed -i -e  's/,\n/\r/g' Matrix.txt")
else:
	print "imporoper use of -P"
