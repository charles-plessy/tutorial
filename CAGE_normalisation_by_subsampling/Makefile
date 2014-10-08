all: subsampling.Rmd
	/usr/bin/Rscript -e "knitr::knit2html('subsampling.Rmd')"
	sed -i 's/\.md/\.html/g' subsampling.html

clean:
	$(RM) -r cache
