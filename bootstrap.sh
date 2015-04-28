#!/bin/sh
echo -n '1) aclocal...'
aclocal
echo 'done.'
echo -n '2) autoheader...'
autoheader
echo 'done.'
echo -n '3) libtoolize...'
libtoolize
echo 'done.'
echo -n '4) autoconf...'
autoconf
echo 'done.'
echo -n '5) automake...'
automake -a -c
echo 'done.'
