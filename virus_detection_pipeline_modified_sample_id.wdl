version 1.0

workflow virus_detection_pipeline {
  input {
    String person_id
    File cram_uri
    File cram_index_uri
    File krakenuniq_db_tar
    File ref_genome
    File ref_genome_index
  }

  call extract_reads {
    input: 
      person_id = person_id,
      cram_uri = cram_uri, 
      cram_index_uri = cram_index_uri, 
      ref_genome = ref_genome, 
      ref_genome_index = ref_genome_index
  }

  Boolean has_ebv = extract_reads.ebv_count > 0
  Boolean has_unmapped = extract_reads.unmapped_count > 0

  if (has_ebv) {
    call convert_to_fastq_ebv {
      input:
        person_id = person_id,
        ebv_raw_reads = extract_reads.ebv_raw_reads
    }

    call qc_fastq_ebv {
      input:
        person_id = person_id,
        ebv_1 = convert_to_fastq_ebv.ebv_1,
        ebv_2 = convert_to_fastq_ebv.ebv_2
    }
  }
   
  if (has_unmapped) {
    call convert_to_fastq_unmapped {
      input:
        person_id = person_id,
        unmapped_reads = extract_reads.unmapped_reads,
        ref_genome = ref_genome
    }
    
    call qc_fastq_unmapped {
      input:
        person_id = person_id,
        unmapped_reads_1 = convert_to_fastq_unmapped.unmapped_reads_1,
        unmapped_reads_2 = convert_to_fastq_unmapped.unmapped_reads_2
    }
  }
  
  if (has_ebv || has_unmapped) {
    call run_krakenuniq {
      input: 
        person_id = person_id,
        reads_ebv1 = qc_fastq_ebv.ebv_1_qced,
        reads_ebv2 = qc_fastq_ebv.ebv_2_qced,
        reads_unmapped1 = qc_fastq_unmapped.unmapped_reads_1_qced,
        reads_unmapped2 = qc_fastq_unmapped.unmapped_reads_2_qced,
        krakenuniq_db = krakenuniq_db_tar
    }
  }
  
  call workflow_complete_flag {
    input:
      person_id = person_id,
      ebv_filtered_reads = extract_reads.ebv_filtered_reads,
      report_ebv = run_krakenuniq.report_ebv,
      report_unmapped = run_krakenuniq.report_unmapped
  }

  output {
    File ebv_filtered_reads = extract_reads.ebv_filtered_reads
    File workflow_complete = workflow_complete_flag.completion_flag
    File? report_ebv = run_krakenuniq.report_ebv
    File? report_unmapped = run_krakenuniq.report_unmapped
  }
}

task extract_reads {
  input {
    String person_id
    File cram_uri
    File cram_index_uri
    File ref_genome
    File ref_genome_index
  }
    
  command <<<
    # 1) Extract EBV reads and QC
    samtools view --reference ~{ref_genome} -h ~{cram_uri} chrEBV > ebv_raw_reads.sam
    samtools view -e "rnext==rname && sclen<=20 && rlen>=120" -F 1024 ebv_raw_reads.sam | cut -f1,4 | gzip > ~{person_id}_ebv_filtered_reads.tsv.gz

    # 2) Extract unmapped read pairs and low quality mapped pairs, with QC to remove garbage
    samtools view -h -C \
      -f 1 -F 256 -F 2048 \
      -e '(flag&4 || flag&8 || mapq==0) && rname!="chrEBV"' \
      --reference ~{ref_genome} \
      ~{cram_uri} \
      -o unmapped_reads.cram
      
    # 3) Count lines
    samtools view -c ebv_raw_reads.sam > ebv_count.txt
    samtools view -c unmapped_reads.cram > unmapped_count.txt
  >>>

  output {
    File ebv_raw_reads = "ebv_raw_reads.sam"
    File ebv_filtered_reads = "~{person_id}_ebv_filtered_reads.tsv.gz"
    File unmapped_reads = "unmapped_reads.cram"
    Int ebv_count = read_int("ebv_count.txt")
    Int unmapped_count = read_int("unmapped_count.txt")
  }

  runtime {
    docker: "mgibio/samtools:v1.21-noble"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 50 SSD"
    preemptible: 10
  }
}

task convert_to_fastq_ebv {
  input {
    String person_id
    File ebv_raw_reads
  }
    
  command <<<
    # Convert EBV SAM to FASTQ
    samtools fastq ~{ebv_raw_reads} \
      -1 ~{person_id}_ebv_1.fastq \
      -2 ~{person_id}_ebv_2.fastq \
      -0 /dev/null -s /dev/null -n
    gzip ~{person_id}_ebv_1.fastq ~{person_id}_ebv_2.fastq
  >>>

  output {
    File ebv_1 = "~{person_id}_ebv_1.fastq.gz"
    File ebv_2 = "~{person_id}_ebv_2.fastq.gz"
  }

  runtime {
    docker: "mgibio/samtools:v1.21-noble"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 10 SSD"
    preemptible: 10
  }
}

task qc_fastq_ebv {
  input {
    String person_id
    File ebv_1
    File ebv_2
  }
  
  command <<<
    num_logical_processors=$(nproc)
    
    fastp \
      -i ~{ebv_1} -I ~{ebv_2} \
      -o ~{person_id}_ebv_1_qced.fastq.gz -O ~{person_id}_ebv_2_qced.fastq.gz \
      --low_complexity_filter \
      --complexity_threshold 30 \
      --disable_adapter_trimming \
      --disable_trim_poly_g \
      --thread $num_logical_processors \
      --json ~{person_id}_fastp_ebv.json \
      --html ~{person_id}_fastp_ebv.html
  >>>
  
  output {
    File ebv_1_qced = "~{person_id}_ebv_1_qced.fastq.gz"
    File ebv_2_qced = "~{person_id}_ebv_2_qced.fastq.gz"
  }
  
  runtime {
    docker: "mgibio/fastp:v0.23.4-noble"
    memory: "4 GB"
    cpu: 2
    disks: "local-disk 10 SSD"
    preemptible: 10
  }
}

task convert_to_fastq_unmapped {
  input {
    String person_id
    File unmapped_reads
    File ref_genome
  }
    
  command <<<
    # Convert unmapped CRAM to FASTQ
    samtools fastq --reference ~{ref_genome} ~{unmapped_reads} \
        -1 ~{person_id}_unmapped_reads_1.fastq \
        -2 ~{person_id}_unmapped_reads_2.fastq \
        -0 /dev/null -s /dev/null -n
    gzip ~{person_id}_unmapped_reads_1.fastq ~{person_id}_unmapped_reads_2.fastq
  >>>

  output {
    File unmapped_reads_1 = "~{person_id}_unmapped_reads_1.fastq.gz"
    File unmapped_reads_2 = "~{person_id}_unmapped_reads_2.fastq.gz"
  }

  runtime {
    docker: "mgibio/samtools:v1.21-noble"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 SSD"
    preemptible: 10
  }
}

task qc_fastq_unmapped {
  input {
    String person_id
    File unmapped_reads_1
    File unmapped_reads_2
  }
  
  command <<<
    num_logical_processors=$(nproc)
  
    fastp \
      -i ~{unmapped_reads_1} -I ~{unmapped_reads_2} \
      -o ~{person_id}_unmapped_reads_1_qced.fastq.gz -O ~{person_id}_unmapped_reads_2_qced.fastq.gz \
      --low_complexity_filter \
      --complexity_threshold 30 \
      --disable_adapter_trimming \
      --disable_trim_poly_g \
      --thread $num_logical_processors \
      --json ~{person_id}_fastp_unmapped.json \
      --html ~{person_id}_fastp_unmapped.html
  >>>
  
  output {
    File unmapped_reads_1_qced = "~{person_id}_unmapped_reads_1_qced.fastq.gz"
    File unmapped_reads_2_qced = "~{person_id}_unmapped_reads_2_qced.fastq.gz"
  }
  
  runtime {
    docker: "mgibio/fastp:v0.23.4-noble"
    memory: "4 GB"
    cpu: 2
    disks: "local-disk 20 SSD"
    preemptible: 10
  }
}

task run_krakenuniq {
  input {
    String person_id
    File? reads_ebv1
    File? reads_ebv2
    File? reads_unmapped1
    File? reads_unmapped2
    File krakenuniq_db
  }

  command <<<
    num_logical_processors=$(nproc)
  
    # The krakenuniq_db is a tar ball file that unpacks to a folder named DBDIR_virus
    pigz -dc ~{krakenuniq_db} | tar xf -
    
    if [ -f "~{reads_ebv1}" ] && [ -f "~{reads_ebv2}" ]; then
      krakenuniq \
        --db DBDIR_virus \
        --preload-size 55G \
        --threads $num_logical_processors \
        --paired \
        --check-names \
        --report-file ~{person_id}_REPORTFILE_ebv.tsv \
        --output ~{person_id}_READCLASSIFICATION_ebv.tsv \
        ~{reads_ebv1} ~{reads_ebv2}
    fi
    
    if [ -f "~{reads_unmapped1}" ] && [ -f "~{reads_unmapped2}" ]; then
      krakenuniq \
        --db DBDIR_virus \
        --preload-size 55G \
        --threads $num_logical_processors \
        --paired \
        --check-names \
        --report-file ~{person_id}_REPORTFILE_unmapped.tsv \
        --output ~{person_id}_READCLASSIFICATION_unmapped.tsv \
        ~{reads_unmapped1} ~{reads_unmapped2}
    fi
  >>>

  output {
    File? report_ebv = "~{person_id}_REPORTFILE_ebv.tsv"
    File? report_unmapped = "~{person_id}_REPORTFILE_unmapped.tsv"
  }

  runtime {
    docker: "monsieurbl/krakenuniq_pigz:v1"
    memory: "64 GB"
    cpu: 8
    disks: "local-disk 100 SSD"
    preemptible: 10
  }
}

task workflow_complete_flag {
  input {
    String person_id
    File ebv_filtered_reads
    File? report_ebv
    File? report_unmapped
  }

  command <<<
    echo "Workflow completed successfully for sample ~{person_id}." > ~{person_id}_workflow_complete.txt
  >>>

  output {
    File completion_flag = "~{person_id}_workflow_complete.txt"
  }

  runtime {
    docker: "google/cloud-sdk:slim"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 1 SSD"
    preemptible: 10
  }
}
