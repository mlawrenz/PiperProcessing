input=$1
sed "/HEADER/d" < $input | sed "s/END/TER/g" > complex-$input

