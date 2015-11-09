#!/bin/bash
for i in out_*.tar ; do
    tar xf $i
    rm $i
done

#for i in out/*.scattered out/*.total ; do
#    mv $i $(echo $i | tr : .)
#done

