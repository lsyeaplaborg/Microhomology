### identify microhomology
scripts_dir='/Volumes/Elements/handover_code/Microhomology/scripts/'
input_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/results/simalirity'
out_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/results/mh_pipline_out'
ref_fa='/Volumes/Elements/handover_code/Microhomology/data/VB18_F3_short.fa'

for file in `ls ${input_dir}`
do 
	echo $file
	python ${scripts_dir}/insertion/MH_insertion_v2.py  ${input_dir}/$file  ${ref_fa}  ${out_dir}
done

