sdrf: sdrf.html

%.html: %.Rmd
	/usr/bin/Rscript -e "knitr::knit2html('$<')"

clean:
	rm *_sdrf.txt
