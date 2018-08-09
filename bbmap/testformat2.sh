#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 23, 2018

Description:  Reads the entire file to find extended information about the format and contents.

Usage:        testformat2.sh <file>


Parameters:

full=t          Process the full file.
speed=f         Print processing time.

printjunk=f     Print headers of junk reads to stdout.
printbarcodes=f Print barcodes to stdout.
printqhist=f    Print quality histogram to stdout.
printihist=f    Print insert size histogram to stdout.

merge=t         Calculate mergability via BBMerge.
sketch=t        (card) Calculate cardinality via BBSketch.
trim=t          Calculate trimmability from quality.


File output parameters (these can be eliminated by setting to null):

junk=junk.txt          Print headers of junk reads to this file.
barcodes=barcodes.txt  Print barcodes to this file.
qhist=qhist.txt        Print quality histogram to this file.
ihist=ihist.txt        print insert size histogram to this file.


Terminology:

Format          File format, e.g. fastq.
Compression     Compression format, e.g. gz.
Interleaved     True if reads are paired in a single file.
MaxLen          Maximum observed read length.
MinLen          Minimum observed read length.
StdevLen        Standard deviation of observed read lengths.
ModeLen         Mode of observed read lengths.
QualOffset      Quality score offset.
NegativeQuals   Number of bases with negative quality scores.

Content         Nucleotides or AminoAcids.
Type            RNA, DNA, or Mixed.
Reads           Number of reads processed.
-JunkReads      Reads with invalid bases or other problems.
-ChastityFail   Reads failing Illumina's chastity filter.
-BadPairNames   Read pairs whose names don't match.

Bases           Number of bases processed.
-Lowercase      Lowercase bases.
-Uppercase      Uppercase bases.
-Non-Letter     Non-letter symbols in bases.
-FullyDefined   A, C, G, T, or U bases.
-No-call        N bases.
-Degenerate     Non-ACGTUN IUPAC symbols.
-Gap            - symbol.
-Invalid        Symbols that are not valid characters for sequence.

GC              GC content: (C+G)/(A+T+U).
Cardinality     Approximate number of unique 31-mers in the file.
Barcodes        Number of observed barcodes.

Mergable        Fraction of read pairs that appear to overlap.
-InsertMean     Average insert size, from merging.
-InsertMode     Insert size mode from, merging.
-AdapterReads   Fraction of reads with adapter sequence, from merging.
-AdapterBases   Fraction of bases that are adapter sequence, from merging.

QErrorRate      Average error rate from quality scores.
-QAvgLog        Logarithmic average quality score.
-QAvgLinear     Linear average quality score.
-TrimmedAtQ5    Fraction of bases trimmed at Q5.
-TrimmedAtQ10   Fraction of bases trimmed at Q10.
-TrimmedAtQ15   Fraction of bases trimmed at Q15.
-TrimmedAtQ20   Fraction of bases trimmed at Q20.

Qhist           Quality score histogram, one line per observed quality bin.
Ihist           Insert size histogram, based on pair merging.
BarcodeList     List of observed barcodes.
JunkList        List of headers of problematic reads.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx2g"
z2="-Xms2g"
EA="-ea"
EOOM=""
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 2000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

testformat() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load pigz
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module load pigz
	fi
	local CMD="java $EA $EOOM $z -cp $CP jgi.TestFormat $@"
#	echo $CMD >&2
	eval $CMD
}

testformat "$@"
