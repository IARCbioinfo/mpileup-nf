#! /usr/bin/env nextflow

params.nsplit = 1

bed = file(params.bedfile)

bam = Channel.fromPath( params.bam_folder+'/*.bam' ).toList()
bai = Channel.fromPath( params.bam_folder+'/*.bam.bai' ).toList()

params.use_file_name = null
sample_names = params.use_file_name ? "FILE" : "BAM"

params.out_folder = "."
params.out_table = null
out_table = params.out_table ? params.out_table : "mpileup_coverage.txt"

fasta_ref = file( params.fasta_ref )
fasta_ref_fai = file( params.fasta_ref+'.fai' )

params.map_qual = 20
params.base_qual = 20 
params.max_DP = 50000

params.help = null

if (params.help) {
    log.info ''
    log.info '-----------------------------------------------------------------'
    log.info 'NEXTFLOW PIPELINE FOR COVERAGE COMPUTATION WITH SAMTOOLS MPILEUP'
    log.info '-----------------------------------------------------------------'
    log.info 'Script for computing coverage at each position in an input bed file, for each sample in an input BAM folder.'
    log.info ''
    log.info ''
    log.info 'Output format looks like the following:'
    log.info '          chr pos K211007X K211008Y K211010Z'
    log.info '          chr1 801957 314 252 287'
    log.info '          chr1 866511 212 128 169'
    log.info ''
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run iarcbioinfo/mpileup-nf [-with-docker] --bedfile bedfile.bed --bam_folder BAM/ --fasta_ref reference.fasta [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder     BAM_DIR                  BAM files directory.'
    log.info '    --fasta_ref      REF_IN_FASTA             Reference genome in fasta format.'
    log.info '    --bedfile        BED FILE                 A BED file for coverage computations.'
    log.info 'Optional arguments: '
    log.info '    --nsplit         INTEGER                  Split the region for calling in nsplit pieces and run in parallel.'
    log.info '    --use_file_name                           Sample names are taken from file names. Default: extracted from the bam file SM tag.'
    log.info '    --out_folder     OUTPUT FOLDER            Output directory, by default is current directory.'
    log.info '    --out_table      OUTPUT TABLE             Output table file, by default "mpileup_coverage.txt".'
    log.info '    --map_qual       INTEGER (PHRED)          Samtools minimum mapping quality.'
    log.info '    --base_qual      INTEGER (PHRED)          Samtools minimum base quality.'
    log.info '    --max_DP         INTEGER                  Samtools maximum coverage before downsampling.'
    log.info ''
    log.info ''
    exit 1
}

process split_bed {

      input:
      file bed

      output:
      file '*_regions' into split_bed mode flatten

      shell:
      '''
      grep -v '^track' !{bed} | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print $1" "$2" "$3}' | bed_cut.r !{params.nsplit}
      '''
}

  process samtools_mpileup {

      tag { region_tag }

      input:
      file split_bed
      file 'BAM/*' from bam
      file 'BAM/*' from bai
      file fasta_ref
      file fasta_ref_fai

      output:
      set val(region_tag), file("${region_tag}.pileup") into pileup

      shell:
      region_tag = split_bed.baseName
      '''
      set -o pipefail
      while read bed_line; do
          samtools mpileup --fasta-ref !{fasta_ref} --region $bed_line --ignore-RG --min-BQ !{params.base_qual} --min-MQ !{params.map_qual} --max-idepth 1000000 --max-depth !{params.max_DP} BAM/*.bam | sed 's/		/	*	*/g' >> !{region_tag}.pileup
      done < !{split_bed}
      '''
  }

  process mpileup2table {

      tag { region_tag }

      input:
      set val(region_tag), file("${region_tag}.pileup") from pileup.filter { tag, file -> !file.isEmpty() }

      output:
      set val(region_tag), file('sample*.txt') into table

      shell:
      '''
      nb_pos=$(wc -l < !{region_tag}.pileup)
      if [ $nb_pos -gt 0 ]; then
          # split and convert pileup file
          pileup2baseindel.pl -i !{region_tag}.pileup
      fi
      '''
  }

  process table2coverage {

      tag { region_tag }

      input:
      set val(region_tag), file('sample*.txt') from table
      file 'BAM/*' from bam
      val sample_names
	
      output:
      file '*_cov.txt' into cov mode flatten

      shell:

      '''
      i=1
      for cur_bam in BAM/*.bam
      do
          if [ "!{sample_names}" == "FILE" ]; then
              # use bam file name as sample name
              bam_file_name=$(basename "${cur_bam%.*}")
              # remove whitespaces from name
              SM="$(echo -e "${bam_file_name}" | tr -d '[[:space:]]')"
          else
              # extract sample name from bam file read group info field
              SM=$(samtools view -H $cur_bam | grep @RG | head -1 | sed "s/.*SM:\\([^	]*\\).*/\\1/" | tr -d '[:space:]')
          fi
          printf "sample$i	$SM\\n" >> names.txt
          i=$((i+1))
      done

      table2coverage.r
      '''
  }

  process merge_table {

      publishDir params.out_folder, mode: 'copy'

      input:
      val out_table
      file all_cov from cov.toList()

      output:
      file "$out_table" into table_merged
	
      shell:
      '''
      head -n1 `ls -1 *cov.txt | head -1` > header.txt
      sed -i '1d' *cov.txt
      cat *cov.txt >> header.txt
      mv header.txt !{out_table}
      '''
  }

