all: analysis.Rmd
	/usr/bin/Rscript -e "knitr::knit2html('analysis.Rmd')"

clean:
	$(RM) -r cache *.fastq *.bam *.fq *.fastq.bz2 *.id *.log tagdust.fa
