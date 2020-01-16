#!/usr/bin/python
#FindReadLengths 09-2017
#Written by Joseph Sarro
#This script fings thelegnth of all reads in a SAM file
#files to be opened and written.  F needs to be changed to match your SAM file before running as does datafile below.  
import os
f = open('S1b.sam', 'r')
l= open('ReadLengths.txt','w')
for lines in f:
	#count +=1
	fields = lines.split("	")
	sq=fields[0]
	#print sq
	#ignore lines with headers
	if sq != "@SQ" and sq != "@RG" and sq != "@PG":
		field1= fields[9]
                length = len(field1)
		print >> l, length
