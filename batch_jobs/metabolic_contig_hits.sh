for file in `find $(pwd) -name '*.gene' -not -path '*/\.*'`; 
do binID="`basename $file | cut -f1 -d '.'`";
geneCount="`grep -o '>' $file | wc -l`";
echo "${binID} ${geneCount}" >> geneCounts.txt;
done

