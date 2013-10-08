all: analysis.Rmd
	/usr/bin/Rscript -e "knitr::knit2html('analysis.Rmd')"

clean:
	$(RM) -r cache *.fastq *.bam *.fq NCms10010.fastq.xz NCms10010.id *.log tagdust.fa
