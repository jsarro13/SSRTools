# SSRTools
The following scripts were generated to form a pipeline that takes sequence capture data and generates a data matrix costing of Short Sequence Repeat(SSR) information. In addition to the novel scripts, two previously written scripts were used. FindReadLengths.py was used to verify our sequence capture had worked correctly, and SnipVCF.py to filter downstream SNPs. Both tools can be found at https://git.io/JTYMb.


sardinesandwich.py: This script can be used as a quality control measure after sequence alignment. It will output four files. One file containing reads that flank 50bp or more both upstream and downstream of the SSR, a file containing reads that only flank 50bp or more upstream of the SSR, a file containing reads that flank 50bp or more only downstream of the SSR, and reads that do neither. Files contain three columns, the read length, flanking length to the left of the SSR, and flanking length to the right of the SSR. The script will also print all counts to terminal.

SamIAm.py: This script will print out sequences that flank and SSR at least 50bp in both directions. These sequences can be used for downstream SNP calling 

GapGenes.py: Measures the length of an SSR after sequence alignment. Two output files are generated. The first prints the SSR, repeat number, and both flanking ends. A second file contains a matrix of count information of SSR lengths.

SNPMatrix.py: This script will filter the matrix output of GapGenes.py. It filters with the following criteria: 1) A contig can contain no more than 3 SSRs 2) An SSR must have at least 4 reads. 3) SSR two and three may not be less than 25% the most abundant SSR.

FormatMatrix.py: Formats the output of SNPMatrix.py into a format that can be read by JoinMap.


An example workflow would include trimming data with trimmomatic and mapping with BWA. After mapping, the SAM files can be processed as follows.

python GapGenes.v3.py -sam <SAM file> -C <SSR info file> -P/S
python SnipMatrix.py <Output file from GapGenes>
python Format_Matrix.py
  
  For SNPs, a VCF file must be generated. It is recommended to use a tool like BCFTools. Before calling SNPs you can filter the SAM file to only include reads mapped within 50bp of the SSR. An example of this workflow would be.
  python SamIAm.py <SAM file> -p <yes/no>
  <Call SNPs with your workflow>
    python Format_VCF-Matrix.py
