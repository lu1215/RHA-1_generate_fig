#!bin/bash

# usage: bash run_generate_fig.sh

# get cutadapt2.9 (preprocess fastq file)
# pip install pyyaml


# [ True / False ]
preprocess=False
del_meta=False

################
## preprocess ##
################

if [ $preprocess = 'True' ]; then

    # all files need to be analysed
    # user should put these files in input folder and rename it as same as name below

    # [ "RPM" / "miRNA_abu" ]
    norm_method="RPM"


    total_og_files=("OG-02" "OG-06" "OG-10" "OG-14" "OG-04" "OG-08" "OG-12" "OG-16" "OG-03" "OG-07")

    for ((i=0; i<${#total_og_files[@]}; i+=1)); do
        og_file="${total_og_files[i]}"
        echo "preprocessing $og_file"
        cd sRNAanalyst/example/script/
        echo input_file=../../../input/${og_file}.fastq.gz >> 22G_preprocess.conf
        echo base_name=${og_file} >> 22G_preprocess.conf
        if [ $norm_method = "RPM" ]; then
            bash 22G_preprocess_RPM.sh
        else
            bash 22G_preprocess.sh
        fi
        mv ../output/${og_file}.csv ../../../input/
        cd -
    done
fi

#######################
## metagene analysis ##
#######################

# metagene input data name
metagene_groups=("OG-02" "OG-06" "OG-10" "OG-14")
for ((i=0; i<${#metagene_groups[@]}; i+=2)); do
    mv "input/${metagene_groups[i]}.csv" sRNAanalyst/example/data/dataset/
    mv "input/${metagene_groups[i+1]}.csv" sRNAanalyst/example/data/dataset/

    og1="${metagene_groups[i]}"
    og2="${metagene_groups[i+1]}"
    echo "metagene analysis: $og1 and $og2"    

    python code/change_metagene_para.py $i 

    cd sRNAanalyst/example/ || exit 1

    # 執行分析
    sh script/22G_analyze.sh

    # 移動輸出圖片
    if [ "$i" -eq 0 ]; then
        mv output/analyze/fig/Metagene_0.png ../../output/fig_h.png
    else
        mv output/analyze/fig/Metagene_0.png ../../output/fig_f.png
    fi

    cd - || exit 1
done


###########################
## scatter plot analysis ##
###########################

python code/run_scatter.py

######################################
## 22G-level & 22G-RNA distribution ##
######################################


G22_og_files=("OG-04" "OG-08" "OG-12" "OG-16" "OG-03" "OG-07")

for ((i=0; i<${#G22_og_files[@]}; i+=1)); do
    python code/process_22G_bedgraph.py input/${G22_og_files[i]}.csv T26E3.3a.1 1410
done

if [ $del_meta = 'True' ]; then
    rm sRNAanalyst/example/data/dataset/OG*
fi