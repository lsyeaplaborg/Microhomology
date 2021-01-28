### identify microhomology
scripts_dir='/Volumes/Elements/code/mh_pipeline/develop/scripts'
input_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/x083_202101/ins_classfication/results/final_res'
out_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/x083_202101/ins_MH/mh_pipeline_out'
ref_fa='/Volumes/Elements/code/mh_pipeline/develop/data/VB18_F3_short.fa'

for file in `ls ${input_dir}`
do 
	echo $file
	python ${scripts_dir}/insertion/MH_insertion_v2.py  ${input_dir}/$file  ${ref_fa}  ${out_dir}
done

### count total reads
filt_clone_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/data/20210125/x083/x083_filt_clones'
out_dir='/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/x083_202101/ins_MH'

python ${scripts_dir}/total_reads.py ${filt_clone_dir} ${out_dir}

### downstream analysis 
pipline_out="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/x083_202101/ins_MH/mh_pipeline_out_all"
meta_dat="/Volumes/Elements/code/mh_pipeline/develop/data/merged_meta.txt"
total_reads_dat="/Volumes/Elements/code/mh_pipeline/develop/data/all_reads.txt"
out_dir="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/x083_202101/ins_MH/downstream_update"
factor_order_mh="WT,Fen1,UNG,53BP1,ATM,H2AX,XLF,PolH,MSH2,Exo1,Pms2,MLH1,Ape2-M,Ape2-F,AID,UNG-MSH2,MSH2-UNG-HE"
min_ins_len=4
category="1,2,3,11,7,12,13"
merge_mh_len=5
prefix="all_ins_g4"

Rscript ${scripts_dir}/insertion/mh_frequncy.R ${pipline_out} ${meta_dat} ${total_reads_dat} ${out_dir} ${factor_order_mh} ${min_ins_len} ${category} ${merge_mh_len} ${prefix}



