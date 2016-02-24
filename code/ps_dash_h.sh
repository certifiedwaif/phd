ps u | awk '{
for ( x=1 ; x<=4 ; x++ ) { printf("%s\t",$x) } 
for ( x=5 ; x<=6 ; x++ ) {  if (NR>1) { printf("%13.2fMb\t",hr=$x/1024) }  
else { printf("\t%s\t",$x)  } 
}
for ( x=7 ; x<=10 ; x++ ) { printf("%s\t",$x) } 
for ( x=11 ; x<=NF ; x++ ) { printf("%s ",$x) } 
print "" 
}'