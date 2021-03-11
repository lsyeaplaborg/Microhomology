### count total reads
filt_clone_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/data'
out_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/meta'
scripts_dir='/Volumes/Elements/handover_code/Microhomology/scripts'

python ${scripts_dir}/total_reads.py ${filt_clone_dir} ${out_dir}
