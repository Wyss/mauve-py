#!/usr/bin/env bash

TESTNAME=$1

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TESTDIR=${DIR}/${TESTNAME}_output

BASESEQ=${DIR}/baseseq.fa
TESTSEQ=${DIR}/test_${TESTNAME}.fa

if [ ! -e ${TESTSEQ} ]; then
  echo -e "\n${TESTNAME} is not a valid test; ${TESTSEQ} does not exist\n"
  echo -e "Valid test names are: "
  TESTRE="test_(.*).fa"
  for f in $(ls ${DIR}); do
    if [[ $f =~ $TESTRE ]]; then
      test_name=${BASH_REMATCH[1]}
      echo -e "\t${test_name}"
    fi
  done
  echo -e
  exit 1
fi

echo "Test output files will be saved to: ${TESTDIR}"

rm -rf $TESTDIR
mkdir $TESTDIR
cd $TESTDIR
$DIR/../bin/progressiveMauveStatic --output=${TESTNAME}.xmfa ${BASESEQ} ${TESTSEQ} 

# Delete all of the .sslist files
rm ${DIR}/*.sslist

exit 0
