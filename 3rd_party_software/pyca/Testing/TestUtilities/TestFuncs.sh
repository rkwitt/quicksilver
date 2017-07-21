
EXIT_PASS=0
EXIT_FAIL=1

function calc {
    echo "scale=4; $*" | bc
}

# If the first arg (probably $?, the output of the last command) is
# nonzero (bad return code), exit with the message given by the second
# arg
function _passtest {
    if [ $1 -ne 0 ]; then
	echo "failing: $2"
	exit $EXIT_FAIL
    fi
}

function _make_and_push_dir {
    if [ ! -d "$1" ]; then
	mkdir "$1"
    fi
    pushd $1
}

# return the return code of command
function _test_execute {
    echo "executing: $*"
    return `$*`
}

# first parameter is command, second is file to pipe output to
function _execute_save_output {
    echo "executing: $1"
    if [ -n $2 ]; then
	$1 &> $2
	_passtest $? "Error executing: $1"
    else
	$1
	_passtest $? "Error executing: $1"
    fi
}

# fail on bad return code
function _execute {
    echo "executing: $*"
    $*
    _passtest $? "Error executing: $*"
}

function _sedsafe_path {
    SEDSAFE=`echo $1 | sed 's/\([\/]\)/\\\\\1/g'`
    echo $SEDSAFE
}

function _print_gnuplot_header {
    PNGFILE=$1
    TITLE=$2
    echo "set terminal png"
    echo "set   autoscale                        # scale axes automatically"
    echo "unset log                              # remove any log-scaling"
    echo "unset label                            # remove any previous labels"
    echo "set xtic auto                          # set xtics automatically"
    echo "set ytic auto                          # set ytics automatically"
    echo "set title \"$TITLE\""
    echo "set xlabel \"Iteration\""
    echo "set ylabel \"Energy\""
    echo "set out \"$PNGFILE\""
}

function _gen_plot {
    DATAFILE=$1
    PNGFILE=$2
    TITLE=$3
    PLOTFILE=plot.gnuplot
    _print_gnuplot_header $PNGFILE $TITLE > $PLOTFILE
    echo "plot    \"$DATAFILE\" using 2:3 title 'total' with linespoints,  \\" >> $PLOTFILE
    echo "	\"$DATAFILE\" using 2:4 title 'image' with linespoints,  \\" >> $PLOTFILE
    echo "	\"$DATAFILE\" using 2:5 title 'vector' with linespoints" >> $PLOTFILE
    echo "set out" >> $PLOTFILE
    _execute gnuplot $PLOTFILE
}

function _gen_plotline {
    DATAFILE=$1
    PLOTCOLS=$2
    LINETITLE=$3
    
    USING=""
    if [ -n $PLOTCOLS ]; then
	USING=" using $PLOTCOLS"
    fi
    TITLE=""
    if [ -n $LINETITLE ]; then
	TITLE=" title '${LINETITLE}'"
    fi
    echo -n "plot    \"$DATAFILE\" $USING $TITLE with linespoints"
}
