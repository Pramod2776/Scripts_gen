#!/usr/bin/env bash
#SBATCH -t 00:01:00
#SBATCH	-p shared
#SBATCH --mem=1G
#SBATCH -J Num_to_Sample_ipsc
#SBATCH --job-name=Num_to_Sample_ipsc
#SBATCH --output=Num_to_Sample_ipsc.out
#SBATCH --array=1-36

PERIOD="3m"

FOLDER_NAME=`cat /home/jic307/jiayechen/16p_Project/input_filenames/${PERIOD}_files.txt | awk -v SLURM_ARRAY_TASK_ID=$((${SLURM_ARRAY_TASK_ID} * 2)) 'NR==SLURM_ARRAY_TASK_ID {print}' | sed "s/[_]*S[0-9]*//g" | sed "s/_R[0-9][_]*[0-9]*.fastq.gz//g" | sed "s/_L[0-9]*//g"`

OLD_FOLDER_NAME=${SLURM_ARRAY_TASK_ID}
DIR="/home/jic307/jiayechen/16p_Project/16p_isoform_mapping_RSEM_result/${PERIOD}"

mv "$DIR/$FOLDER_NAME/RSEM_Quant.isoforms.results" "$DIR/$FOLDER_NAME/${FOLDER_NAME}_RSEM_Quant.isoforms.results"
#mv "$DIR/$OLD_FOLDER_NAME" "$DIR/$FOLDER_NAME"

#for FILE in $(ls "$DIR/$FOLDER_NAME" | grep "[0-9]*-[0-9]*-[0-9]*-[0-9]*-[0-9]*-[0-9]*")
#do
#	TEMP=`echo $FILE | sed "s/[0-9]*-[0-9]*-[0-9]*-[0-9]*-[0-9]*-[0-9]*/$FOLDER_NAME/g"`
#	mv "$DIR/$FOLDER_NAME/$FILE" "$DIR/$FOLDER_NAME/$TEMP" 
#done

#STAR_OUTPUT=`ls ${DIR}/${FOLDER_NAME}/*Log.final.out`
#mv $STAR_OUTPUT "$DIR/$FOLDER_NAME/${FOLDER_NAME}.STARLog.final.out"

#SORTED_STAR=`ls ${DIR}/${FOLDER_NAME}/*Aligned.sortedByCoord.out.bam`
#mv $SORTED_STAR "$DIR/$FOLDER_NAME/${FOLDER_NAME}.STARAligned.sortedByCoord.out.bam"

#FULL_PATH=`ls -d $PWD/16p_isoform_mapping_RSEM_result/${PERIOD}/${FOLDER_NAME}`
#echo $FULL_PATH

#echo "${FOLDER_NAME} ${FULL_PATH}" >> /home/jic307/jiayechen/16p_Project/samples.txt	
