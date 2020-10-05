#!/usr/bin/env nextflow

params.reads = 'data/*.fq.gz'
params.adapt = 'data/adapters.fasta'
params.results = 'results'

files = Channel.fromPath(params.reads)
adapters = file(params.adapt)


process gunzip {
    input:
        file f from files

    output:
        file "*.fq" into unzip_files

    script:
		"""
		gunzip -f "${f}"
		"""

}

process adapter_trimming {
    input:
        val file from unzip_files
        file 'adapters.fasta' from adapters

    output:
        file "*.fastq" into adapt_trimmed

    script:
		"""
		scythe -q sanger -a adapters.fasta -o "${file.baseName}.fastq" $file
		"""

}

process quality_trimming {
    input:
        file fastq from adapt_trimmed

    output:
        file "*" into trimmed

    script:
        """
        sickle se -f $fastq -t sanger -o "${fastq.baseName}" -q 20
        trim-low-abund.py -C 3 -Z 18 -V -M 100e9 "${fastq.baseName}"
        mv "${fastq.baseName}".abundtrim "${fastq.baseName}"
        """
        // khmer, which trim_low_abundance comes from, suggests here khmer.readthedocs.io/en/v2.1.1/user/choosing-table-sizes.html
        // for -M 256G for ~300 Gbp of soil metagenomes. But summit nodes have 4.84 * 24 = 116 GB on normal nodes, and only five 42.7 GB * 48 = 2049 GB nodes in queue smem
}

process sourmash_compute {
    input:
        file trimmed

    output:
        file "*.sig" into sourmash_compute

    script:
        """
        sourmash compute -f $trimmed --scaled 1000 -k 31
        """
        // defualt = -k 31 -n 500
        //  -k 21,31,51 
        // for some downstream uses, may want to use --containment
}

process sourmash_compare {

    input:
        file '*.sig' from sourmash_compute.toList()

    output:
        file "cmp" into sourmash_compare
        file "cmp.labels.txt" into labels
     
    script:
        """
        sourmash compare -k 31 *.sig -o cmp
        """
        //  -p 8 for multi-threading, also add process.$sourmash_compare.cpus = 8 to config
}

process sourmash_plot {
    publishDir params.results, mode: 'copy'

    input:
        file "cmp" from sourmash_compare
        file "cmp.labels.txt" from labels

    output:
        file "cmp.matrix.png" into plots
        file "cmp.csv" into plot_csv

    script:
        """
        sourmash plot cmp --labels --csv cmp.csv
        """
        // if env dosen't have it, this will fix matplotlib errors
        //         mkdir -p $HOME/.config/matplotlib
        //         echo "backend : Agg" > $HOME/.config/matplotlib/matplotlibrc
}

