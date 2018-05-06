for an in 0 1 2 3 ; do
 for am in 0 2 ; do
  for bn in 0 1 2 3 ; do
   for bm in 0 2 ; do
    for bf in 0 1 2 3 ; do
     c99 -o ../../tests/raytests/z/AN${an}M${am}BN${bn}M${bm}F${bf} -DAN${an} -DAM${am} -DBN${bn} -DBM${bm} -DBF${bf} -ggdb -O2 -I .. vswf_translation_test_rays.c ../translations.c ../vswf.c ../gaunt.c ../legendre.c -lgsl -lm -lblas
    done
   done
  done
 done
done

export NORM="129"
for an in 0 1 2 3 ; do
 for am in 0 2 ; do
  for bn in 0 1 2 3 ; do
   for bm in 0 2 ; do
    for bf in 0 1 2 3 ; do
     ./AN${an}M${am}BN${bn}M${bm}F${bf} an${an}m${am}bn${bn}m${bm}f${bf}_${NORM} 1 0 ${NORM}
     sed -e 's/\./,/g' an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.0 > an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.0.tsv
     sed -e 's/\./,/g' an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.1 > an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.1.tsv
     sed -e 's/\./,/g' an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.2 > an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.2.tsv
#    for i in an${an}m${am}bn${bn}m${bm}_130.? ; do
#     sed -e 's/\./,/g' $i >${i}.tsv
#    done
    done
   done
  done
 done
done


export NORM="129"
export smer="125"
export an=1 am=0 bn=1 bm=0 bf=0
     c99 -o ../../tests/raytests/$smer/AN${an}M${am}BN${bn}M${bm}F${bf} -DAN${an} -DAM${am} -DBN${bn} -DBM${bm} -DBF${bf} -ggdb -I .. vswf_translation_test_rays.c ../translations.c ../vswf.c ../gaunt.c ../legendre.c -lgsl -lm -lblas
     cd ../../tests/raytests/$smer
     ./AN${an}M${am}BN${bn}M${bm}F${bf} an${an}m${am}bn${bn}m${bm}f${bf}_${NORM} 1 0 ${NORM}
     sed -e 's/\./,/g' an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.0 > an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.0.tsv
     sed -e 's/\./,/g' an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.1 > an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.1.tsv
     sed -e 's/\./,/g' an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.2 > an${an}m${am}bn${bn}m${bm}f${bf}_${NORM}.2.tsv
     cd -
     