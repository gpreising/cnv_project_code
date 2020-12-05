# Create lists for looping
awk '{print $1}' ON* | sort -u | sed '/NA/d' | sed '/PROBE/d' > nr_ON_list.txt
awk '{print $1}' MZ* | sort -u | sed '/NA/d' | sed '/PROBE/d' > nr_MZ_list.txt
ONlist='nr_ON_list.txt' # repeat loop for each probe recognized as ON gain
MZlist='nr_MZ_list.txt' # repeat loop for each probe recognized as MZ gain

#create blank probe category files
touch nr_ON_probe_category_temp.txt
touch nr_MZ_probe_category_temp.txt

# Now for each probe we find each technology that identified it as a gain (bias for NGS) in ON
while read z; do
	grep ${z} ON* | sed 's/:/ /g' | awk '{print $1}' | sort | tr "\n" "_" | sed 's/ON_bias_//g' | sed 's/ON_gains_//g' | sed 's/list.txt_//g' | awk '{print "ON_" $1 "list"}' >> nr_ON_probe_category_temp.txt
done < $ONlist
mv nr_ON_probe_category_temp.txt nr_ON_probe_category.txt
paste nr_ON_probe_category.txt nr_ON_list.txt > nr_ON_tech_overlaps.txt

# Now for each probe we find each technology that identified it as a gain (bias for NGS) in MZ
while read y; do
	grep ${y} MZ* | sed 's/:/ /g' | awk '{print $1}' | sort | tr "\n" "_" | sed 's/MZ_bias_//g' | sed 's/MZ_gains_//g' | sed 's/list.txt_//g' | awk '{print "MZ_" $1 "list"}' >> nr_MZ_probe_category_temp.txt
done < $MZlist
mv nr_MZ_probe_category_temp.txt nr_MZ_probe_category.txt
paste nr_MZ_probe_category.txt nr_MZ_list.txt > nr_MZ_tech_overlaps.txt

#combine lists to create master list
cat nr_ON_tech_overlaps.txt nr_MZ_tech_overlaps.txt > gain_categories.txt
#group gains regardless of species bias
sed 's/ON_//g' gain_categories_redo.txt | sed 's/MZ_//g' > gain_categories2.txt
