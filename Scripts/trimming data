# loop through all files in directory ending with .fastq.gz
for file in t*.fasta; do
  # run fastx_trimmer -z command on each file
#   fastq_to_fasta -i "$file" -o "${file}.fasta"
  fastx_trimmer -i "$file" -o "trim${file}"
done
