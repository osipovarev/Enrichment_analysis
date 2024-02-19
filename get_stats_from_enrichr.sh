 #!/usr/bin/env bash


terms=$1
# names of files to extract info from, comma-sparated
infiles=$2
# column values (choose two): 4:adjused p-values; 7:Odds ratio; 8:Combined score
col1=$3
col2=$4

# echo -e "term\tup_pval\tup_odds\tdown_pval\tdown_odds"

while read -r line; 
do
   out_line=$line
   for fi in $(echo -e $infiles | sed 's/,/ /g');
   do
      val_col1=$(awk -F"\t" -v v=$col1 -v l="$line" '$1==l {print $v}' $fi );
      val_col2=$(awk -F"\t" -v v=$col2 -v l="$line" '$1==l {print $v}' $fi );
      out_line+="\t$val_col1\t$val_col2"
   done

   echo -e "$out_line";
done < $terms