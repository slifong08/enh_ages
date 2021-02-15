 awk 'BEGIN{FS=OFS=","} {if (a[$1]<$3) {a[$1]=$3; data[$1]=$0}} END{for (i in a) print data[i]}' file
