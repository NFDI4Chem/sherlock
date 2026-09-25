UNAME=`uname`
FILES="Makefile solve m_edit"
if ( echo $UNAME | grep 'CYGWIN' > /dev/null )
then
	DIR='Cygwin'
elif ( echo $UNAME | grep 'Linux' > /dev/null )
then
	DIR='Linux'
elif ( echo $UNAME | grep 'Darwin' > /dev/null )
then
	DIR='MacOSX'
else
	echo "Unknown system $UNAME"
	exit 1
fi
echo $DIR
(cd $DIR ; cp Makefile solve m_edit ..)
chmod u+x solve m_edit
make

