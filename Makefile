all: README.md
	pandoc README.md > README.html
	sed -i 's/\.md/\.html/g' README.html
