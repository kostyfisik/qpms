#!/bin/bash
jobno=$1
# job numbers start from zero, lines from one
freqno=$[jobno + 1]
freq=$(sed "${freqno}q;d" OmegaList)
epfiletmp="EPFile_30ring-xz_$(echo -n $freq | tr . :)"
epfilec="EPFile_30ring-xz_$(echo -n $freq)"
rm ${epfiletmp}.scattered.*
rm ${epfiletmp}.total.*
ln -s EPFile_30ring-xz "$epfiletmp"
/lib64/ld-linux-x86-64.so.2 --library-path scuff-em/lib/ scuff-em/bin/scuff-scatter --Omega $freq --EPFile $epfiletmp <Args
ST=$?
mv "${epfiletmp}.total" out/${epfilec}.total
mv "${epfiletmp}.scattered" out/${epfilec}.scattered
tar cf "out_${jobno}.tar" "out/${epfilec}.total" "out/${epfilec}.scattered"
hostname >out/${freqno}.out
echo $ST $freq >>out/${freqno}.out
rm $epfiletmp "out/${epfilec}.total" "out/${epfilec}.scattered" 
exit $ST
