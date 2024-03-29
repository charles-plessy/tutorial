Tutorials for analysing CAGE and Deep-RACE data.
================================================

Various tutorials on how to analyse
[CAGE](https://en.wikipedia.org/wiki/Cap_analysis_gene_expression) data.

 * [Deep-RACE](./Deep-RACE1/Deep-RACE1.md) (work in progress)
 * [CAGE differential analysis 1](./CAGE_differential_analysis1/analysis.md) (work in progress)
 * [CAGE differential analysis 2](./CAGE_differential_analysis2/analysis.md)
 * [Simple use of FANTOM5 SDRF files](./FANTOM5_SDRF_files/sdrf.md)
 * [Normalisation of CAGE libraries by sub-sampling](./CAGE_normalisation_by_subsampling/subsampling.md)
 * [Demultiplex nanoCAGE data using TagDust 2](./Demultiplexing/Demultiplexing.md)

These tutorials are designed to be executed on a Linux system's command line
interface (also called _Terminal_ or _shell_).  I recommend the book _[The Linux
Command Line][]_, by William E. Shotts, Jr, January 2012, [no starch press][]
to people not familiar with entering commands on the keyboard.  The [_missing
semseter_](https://missing.csail.mit.edu/) course of MIT looks good as well.

[The Linux Command Line]: https://linuxcommand.org/tlcl.php "A Complete Introduction"
[no starch press]: https://nostarch.com/tlcl.htm "the finest in geek entertainment"

The programs used are assumed to be installed in advance.  On the
[Debian](https://www.debian.org) operating system, many of them (BWA, SAMtools,
BEDTools, ...) are available pre-packaged and will be installed (altogether
with many other programs) by the command `apt-get install med-bio`.

Other software have to be downloaded and installed by hand.  Place them in the
`bin` directory in your home directory, and set their executable property in
order to use them.  If you had to create the `bin` directory, it will only be
taken into account at your next connection (see
[stackoverflow](https://stackoverflow.com/questions/16366986/adding-bin-directory-in-your-path)
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
wget http://genome.gsc.riken.jp/osc/english/software/src/tagdust.tgz   # download TagDust
tar xvf tagdust.tgz   # unpack TagDust
cd tagdust            # enter the freshly tagdust directory created by TagDust
make                  # compile the program
cp tagdust ~/bin      # copy tagdust to the 'bin' directory in your home directory
```

Frequent problems
-----------------

### Command not found.

It is not enough to compile a program.  The command-line interface needs to
find them, and by default it does not search in the current work directory.

A very good explanation is in _[The Linux Command Line][]_'s chapter 24,
section _Script File Location_.  Here is a brief summary.

The standard way to make programs accessible is to add them to one of a set of
pre-defined directories that are collectively called the _PATH_.  For
system-wide installations, the directory is usually `/usr/bin`.  For local
installations by a single user, the directory is usually called `bin`, in the
_home_ directory, also accessible via the shortcut `~/bin`.  If it does not exist,
it can be created like any other directory, but it may  be necessary to log out
and in again in order for the system to recognise this directory in the _PATH_.

In addition, the program needs to have the executable permissions.  These can
be given with the `chmod` command (see _[The Linux Command Line][]_'s chapter
24, section _Executable Permissions_.), or via the file navigator of the
desktop graphical interface.

Lastly, it is possible to run a program that is not in the _PATH_.  For this,
just indicate in which directory it is.  The current directory is always
aliased to `.`, so to run a program called `myscript` that is in the current
directory, type `./myscript`.  (The comment above about executable permissions
still applies).

### What is that sponge ?

`sponge` is a command from the [moreutils](https://joeyh.name/code/moreutils/)
collection, that I use frequently.  On Debian systems, it is easy to install
via the [moreutils](packages.debian.org/moreutils) package.

The goal of `sponge` is to solve the following problem: when one file is read,
piped to a command, and the result is redirected to the file itself, the
contents are not updated as expected, but the file is deleted.  This is because
at the very beginning of the command, the file receiving the redirection is
transformed in an empty file before its contents are even read.  For example,
with a file called `example.fq`:

```
cat example.fq | fastx_trimmer -f 11 > example.fq          # Deletes the file.
cat example.fq | fastx_trimmer -f 11 | sponge example.fq   # Trims the first 10 nucleotides.
```

Without `sponge`, one would need to create a temporary file (which is actually
what `sponge` does in a more proper way behind the scene).

```
cat example.fq | fastx_trimmer -f 11 > example.tmp.fq
mv example.tmp.fq example.fq
```
