
TESTNAME=$1

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TESTDIR=${DIR}/${TESTNAME}_test_output

BASEGENOME=${DIR}/mds42_full.fa
TESTGENOME=${DIR}/mds42_full_${TESTNAME}.fa

if [ ! -e ${TESTGENOME} ]
  then
    echo "${TESTNAME} is not a valid test; ${TESTGENOME} does not exist"
    exit 1
fi


echo "Test output files will be saved to: ${TESTDIR}"

rm -rf $TESTDIR
mkdir $TESTDIR
cd $TESTDIR
$DIR/../bin/progressiveMauveStatic --output=mds42_full_${TESTNAME}.xmfa ${BASEGENOME} ${TESTGENOME}

exit 0
