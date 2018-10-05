#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 17, 2018
This script requires at least 17GB RAM.
It is designed for NERSC and uses hard-coded paths.

Description:  Removes all reads that map to the human genome with at least 88% identity after quality trimming.
This is more aggressive than removehuman.sh and uses an unmasked human genome reference.
It removes roughly 99.99% of human 2x150bp reads, but may incur false-positive removals.
NOTE!  This program uses hard-coded paths and will only run on Nersc systems unless you change the path.

Usage:  removehuman.sh in=<input file> outu=<clean output file>

Input may be fasta or fastq, compressed or uncompressed.

Parameters:
threads=auto        (t) Set number of threads to use; default is number of logical processors.
overwrite=t         (ow) Set to false to force the program to abort rather than overwrite an existing file.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
trim=t              Trim read ends to remove bases with quality below minq.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
untrim=t            Undo the trimming after mapping.
minq=4              Trim quality threshold.
ziplevel=2          (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
outm=<file>         File to output the reads that mapped to human.
path=               Set the path to an indexed human genome.

***** All BBMap parameters can be used; run bbmap.sh for more details. *****

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
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
NATIVELIBDIR="$DIR""jni/"

z="-Xmx16000m"
z2="-Xms16000m"
EA="-ea"
EOOM=""
set=0

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

function removehuman() {
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
	local CMD="java $EA $z $z2 -cp $CP align2.BBMap minratio=0.75 maxindel=8 bwr=0.22 bw=26 minhits=1 path=/global/projectb/sandbox/gaag/bbtools/hg19 build=2 pigz unpigz zl=6 qtrim=r trimq=10 untrim idtag usemodulo printunmappedcount ztd=2 maxsites=1 k=14 tipsearch=0 kfilter=25 bloomfilter $@"
	echo $CMD >&2
	eval $CMD
}

removehuman "$@"
