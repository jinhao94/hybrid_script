less $1 | grep '>' | perl -a -F"\s+" -lne '@F[0]=~/>(.*)_(\d+)/; $o=$1."_".$2; print "$o\t$1"' 
