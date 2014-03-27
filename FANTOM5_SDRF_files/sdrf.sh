# Return all libraries.

function SDRFlib {
for SDRF in $(find ${F5-.} -name '*sdrf.txt')
do
  grep "$1" $SDRF | cut -f14
done
}

# Returns the description of a library.

function SDRFdesc {
for SDRF in $(find ${F5-.} -name '*sdrf.txt')
do
  grep "$1" $SDRF | cut -f3
done
}
