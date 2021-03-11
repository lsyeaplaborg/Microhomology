### downstream analysis
scripts_dir='/Volumes/Elements/handover_code/Microhomology/scripts' 
pipline_out="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/results/mh_pipline_out"
meta_dat="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/meta/merged_meta.txt"
total_reads_dat="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/meta/all_reads.txt"
out_dir="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/insertion_mh_20210311/results/mh_downstream"
factor_order_mh="WT,Fen1,UNG,53BP1,ATM,H2AX,XLF,PolH,MSH2,Exo1,Pms2,MLH1,Ape2-M,Ape2-F,AID,UNG-MSH2,MSH2-UNG-HE"
min_ins_len=4
category="1,2,3,11,7,12,13"
merge_mh_len=5
prefix="all_ins_g4"

Rscript ${scripts_dir}/insertion/mh_frequncy.R ${pipline_out} ${meta_dat} ${total_reads_dat} ${out_dir} ${factor_order_mh} ${min_ins_len} ${category} ${merge_mh_len} ${prefix}



