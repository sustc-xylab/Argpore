#!/bin/bash
## ARGpore is designed to identify ARGs and their host populations in nanopore dataset
##Author: Yu XIA - 2017-10-12
##Email: xiay@sustc.edu.cn
##version 1.0
set -e

#### usage info ####
show_help() {
cat << EOF
Usage: ${0##*/} 
version 1.0
arguments:
	-h	display this help 

	-f 	1D.fasta generated by nanopore sequncing as input for argpore

	-s	similarity cutoff [0-1] for filtering lastal results
		default similarity cutoff is 0.6

	-l	alignment length cutoff [0-1] for filtering lastal results
		default alignment length cutoff is 0.9

	-o	prefix for the name of the output files. 
		By default, the output files are created in the same directory as the input fasta as below: 
			input_taxa.tab
			input_arg.tab
			input_arg.w.taxa.tab
		To change the prefix of output files, specify the suffix name using this option

	-t 	number of threads used for lastal, default t=1


output files:
	input_arg.tab	nanopore reads with valid ARGs hits
	input_taxa.tab	nanopore reads with valid taxonomy assignment
	input_arg.w.taxa.tab	ARGs-containing nanopore reads with valid taxonomy assignment

Example usage: bash argpore.sh -f test.fa 
EOF
}

####################
# define arguments
####################
OPTIND=1  # Reset in case getopts has been used previously in the shell.

# initialize your own variables:
N_threads="1"
Input_fa=""
Lencuoff="0.9"
Simcutoff="60"
Output=$Input_fa
nowt=`date +%Y-%m-%d.%H:%M:%S`;
# nowt="2018-05-17.08:52:33"
# the DIR of argpore scirpt
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# echo $DIR


while getopts "t:f:l:s:o:h" opt; do
	case "$opt" in
		h|--help)
			show_help
			exit 0
			;;
		t)
			N_threads=$OPTARG
			;;
		f)
			Input_fa=$OPTARG
			;;
		l)
			Lencuoff=$OPTARG
			;;
		s)
			Simcutoff=`echo "$OPTARG*100"|bc`
			;;
		o)
			Output=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 0
			;;
		'?')
			show_help >&2
			exit 1
			;;
		-?*)
			print 'Warning: Unknown option (ignored) : %s\n' "$1" >&2
			exit 0
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
		*) # default case: if no more options then break out of the loop
			break
			
	esac
done

if [ -z "$Input_fa" ]
then
	echo "No input fasta, -f must be specified"
	exit
fi
if [ -z "$Output" ]
then
	Output=$Input_fa
fi
if [ -z "$Simcutoff" ]
then
	Simcutoff="60"
fi
if [ -z "$Lencuoff" ]
then
	Lencuoff="0.9"
fi
if [ -z "$N_threads" ]
then
	N_threads=1
fi
shift "$((OPTIND-1))"
# echo $Input_fa $Simcutoff $Lencuoff $N_threads $Output $DIR $nowt

# subset the name of the $Input_fa
# $Input_fa : including input.fa path while $Input_fa2 only contain name
myarray=(`echo $Input_fa| tr "/" " "`) 
Input_fa2=${myarray[-1]}



#######
# ARG quantification
# LAST against the SARG-nt
########

echo "start argpore @ `date +"%Y-%m-%d %T"`"
echo "search $Input_fa2 agaisnt SARG-nt with similarity cutoff $Simcutoff and alignment length cutoff $Lencuoff using $N_threads threads"
# mkdir -p ./tmp

lastal -s 2 -T 0 -Q 0 -a 1 -P $N_threads -f BlastTab ${DIR}/ARGs_database_renamed.fnt.subset.lastindex $Input_fa > /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast
echo "finish searching againt SARG-nt database"
echo "parsing SARG-nt results"
grep -v "#" /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast > /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast.modified

ruby ${DIR}/BlastTab.addlen.rb -s -f ${DIR}/ARGs_database_renamed.fnt.subset < /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast.modified > /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast2

ruby ${DIR}/BlastTab.addlen.rb -f $Input_fa < /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast2 > ${Input_fa}_${nowt}_sarg.tab

echo "finish parsing SARG-nt results"

#################
# taxonomy profile
# lastal 1D.fa or 2D.fa agaisnt markers.fasta
#################
echo " "
echo "search $Input_fa2 agaisnt MetaPhlan 2.0 markergene database with similarity cutoff $Simcutoff and alignment length cutoff $Lencuoff using $N_threads threads"

lastal -s 2 -T 0 -Q 0 -a 1 -b 1 -q 2 -P $N_threads -f BlastTab ${DIR}/markers.lastindex $Input_fa > /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast3

echo "finish searching agaisnt markergene database"
echo "parsing markergene results"
grep -v "#" /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast3 > /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast3.modified

ruby ${DIR}/BlastTab.addlen.rb -s -f ${DIR}/markers.fasta < /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast3.modified > /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast4

ruby ${DIR}/BlastTab.addlen.rb -f $Input_fa < /tmp/argpore_${nowt}_${Input_fa2}_tmp.blast4 > ${Input_fa}_${nowt}_marker.tab

echo "finish parsing markergene results"


# STEP Three: summary in nanopore.summary.R
echo " " 
echo "start parsing results in R"
Rscript ${DIR}/argpore.R $Input_fa $Simcutoff $Lencuoff $Output $nowt $DIR --save


echo ""
echo "finish argpore @ `date +"%Y-%m-%d %T"`"




