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
        """
}

process sourmash_compute {
    input:
        file trimmed

    output:
        file "*.sig" into sourmash_compute

    script:
        """
        sourmash compute -f $trimmed
        """
}

process sourmash_compare {

    input:
        file '*.sig' from sourmash_compute.toList()

    output:
        file "cmp" into sourmash_compare
        file "cmp.labels.txt" into labels
     
    script:
        """
        sourmash compare *.sig -o cmp
        """
}

process sourmash_plot {
    publishDir params.results, mode: 'copy'

    input:
        file "cmp" from sourmash_compare
        file "cmp.labels.txt" from labels

    output:
        file "cmp.matrix.png" into plots

    script:
        """
        mkdir -p $HOME/.config/matplotlib
        echo "backend : Agg" > $HOME/.config/matplotlib/matplotlibrc
        sourmash plot cmp --labels
        """
}

process sourmash_tsv {
    publishDir params.results, mode: 'copy'

    input:
        file "cmp" from sourmash_compare
        file "cmp.labels.txt" from labels

    output:
        file "cmp.tsv" into tsv
        file "cmp.labels.txt" into tsv_labels

        """
        #!/usr/bin/env python
        import numpy
        M = numpy.load(open('cmp', 'rb'))
        numpy.savetxt('cmp.tsv', M, delimiter="\t")
        """
}
