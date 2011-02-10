#!/bin/sh
echo -n '1) aclocal...'
aclocal
echo 'done.'
echo -n '2) autoheader...'
autoheader
echo 'done.'
echo -n '3) autoconf...'
autoconf
echo 'done.'
echo -n '4) automake...'
automake -a -c
echo 'done.'
