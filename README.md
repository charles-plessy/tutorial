Tutorials for analysing CAGE and Deep-RACE data.
================================================

Various tutorials on how to analyse
[CAGE](https://en.wikipedia.org/wiki/Cap_analysis_gene_expression) data.

 * Deep-RACE (in preparation)
 * CAGE differential analysis2 (in preparation)
 * [CAGE differential analysis2](./CAGE_differential_analysis2/analysis.md)

These tutorials are designed to be executed on a Linux system's command line
interface (also called _Terminal_ or _shell_).  I recommend the book _[The Linux
Command Line][]_, by William E. Shotts, Jr, January 2012, [no starch press][]
to people not familiar with entering commands on the keybord.

[The Linux Command Line]: http://linuxcommand.org/tlcl.php "A Complete Introduction"
[no starch press]: http://nostarch.com/tlcl.htm "the finest in geek entertainment"

The programs used are assumed to be installed in advance.  On the
[Debian](http://www.debian.org) operating system, many of them (BWA, SAMtools,
BEDTools, ...) are available pre-packaged and will be installed (altogether
with many other programs) by the command `apt-get install med-bio`.

Other software have to be downloaded and installed by hand.  Place them in the
`bin` directory in your home directory, and set their executable property in
order to use them.  If you had to create the `bin` directory, it will only be
taken into account at your next connection (see
[stackoverflow](http://stackoverflow.com/questions/16366986/adding-bin-directory-in-your-path)
for alternatives).

Here is for example how to download, compile and install the
[tagdust](http://genome.gsc.riken.jp/osc/english/software/src/tagdust.tgz)
software.  By convention, we will download the software in a directory called
`src`.  _Compiling_ means to produce the executable program suitable for your
computer, using the [source code](https://en.wikipedia.org/wiki/Source_code)
that was downloaded.  On Debian systems, the programs necessary for compiling a
program made in the C programming language can be installed through the
`build-essential` package.

```
cd                    # move back to the home directory
mkdir -p src          # create the src directory if it did not exist.
cd src                # enter the src directory
wget http://genome.gsc.riken.jp/osc/english/software/src/tagdust.tgz   # download tagdust
tar xvf tagdust.tgz   # unpack tagdust
cd tagdust            # enter the freshly tagdust directory created by tagdust
make                  # compile the program
cp tagdust $HOME/bin  # copy tagdust to the 'bin' directory in your home directory
```
