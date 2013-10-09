Tutorials for analysing CAGE and Deep-RACE data.
================================================

Various tutorials on how to analyse
[CAGE](https://en.wikipedia.org/wiki/Cap_analysis_gene_expression) data.

 * Deep-RACE (in preparation)
 * CAGE differential analysis2 (in preparation)
 * [CAGE differential analysis2](./CAGE_differential_analysis2/analysis.md)

These tutorials are designed to be executed on a Linux system's command line
interface (also called _Terminal_or _shell_).  The programs used are assumed to
be installed in advance.  On the [Debian](http://www.debian.org) operating
system, many of them (BWA, SAMtools, BEDTools, ...) are available pre-packaged
and will be installed (altogether with many other programs) by the command
`apt-get install med-bio`.

Other software have to be downloaded and installed by hand.  Place them in the
`bin` directory in your home directory, and set their executable property in
order to use them.  If you had to create the `bin` directory, it will only be
taken into account at your next connection (see
[stackoverflow](http://stackoverflow.com/questions/16366986/adding-bin-directory-in-your-path)
for alternatives).
