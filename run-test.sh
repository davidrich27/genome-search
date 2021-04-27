##################################################################################
# NAME:     run_test.sh
# PURPOSE:  Runs naive and optimized nucleotide search and verify produce same result.
##################################################################################

function ECHO_AND_RUN
{
   echo "${@}"
   "${@}"
}

NUM_ARGS=$#
if (( NUM_ARGS == 2 ))
then 
   GENOME_FILE=$1
   SHORTREAD_FILE=$2
else
   echo "run_test.sh: Runs naive and optimized nucleotide search and verify produce same result."
   echo "Usage: <genome_file> <shortread_file>"
   exit 1
fi 

NUC_NAIVE="python nucleotide-search-naive.py"
NUC_OPT="./nucleotide-search"

TEST_NAIVE="out/test-naive.out"
TEST_OPT="out/test-opt.out"

$NUC_NAIVE $GENOME_FILE $SHORTREAD_FILE > $TEST_NAIVE 
$NUC_OPT $GENOME_FILE $SHORTREAD_FILE | grep -v "Took" > $TEST_OPT

TEST=$(diff -q $TEST_NAIVE $TEST_OPT)

if [ -z $TEST ]
then 
   echo "Test passed."
else 
   echo "Test failed."
fi
