#!/usr/bin/env zsh
set -e #-x

function set_env () {
    # @SET_ENV
    # Define $_srcdir and $_bindir
    echo "Where shall we put the source code [/opt/share]: " | tr -d "\n"
    read _srcdir
    eval _srcdir=${_srcdir}
    export _srcdir=${_srcdir:-/opt/share}
    [[ -d "$_srcdir" ]] || mkdir -p $_srcdir

    echo "Where shall we install the programme [/usr/local/bin]: " | tr -d "\n"
    read _bindir
    eval _bindir=${_bindir}
    export _bindir=${_bindir:-/usr/local/bin}
}

function select_blast () {

    INSTALL_DIR=$(pwd)
    _blast_dir=/usr/local
    cd $_blast_dir

    _rmbl_url=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/LATEST/
    _blast_url=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

    _rmbl_idx=$(wget -nv $_rmbl_url -O /dev/fd/1)
    _blast_idx=$(wget -nv $_blast_url -O /dev/fd/1)

    _rmbl_exec=(=$(awk -F'"' '/  [[:digit:]]/ {print $2}' < <(echo $_rmbl_idx)))
    _blast_exec=(=$(awk -F'"' '/  [[:digit:]]/ {print $2}' < <(echo $_blast_idx)))

    # get the ncbi-blast+ binaries
    echo ""
    echo "Choose the approbriate binary distribution: "
    echo ""
    for ((i=2;i<${#_blast_exec}+1;++i)); do
        echo " [$(($i - 1))] ${${_blast_exec[$i]}##*/}"
    done
    echo ""
    read _blast_i
    _blast_exec_i=$_blast_exec[$(($_blast_i + 1))]
    sudo wget -nc $_blast_exec_i
    _blast_exec_i=$(find . -name ncbi-blast-*.tar.gz -print)
    sudo tar zxkf $_blast_exec_i && sudo rm -f $_blast_exec_i

    # get the rmblast binaries
    echo ""
    echo "Choose the approbriate binary distribution: "
    echo ""
    for ((i=2;i<${#_rmbl_exec}+1;++i)); do
        echo " [$(($i - 1))] ${${_rmbl_exec[$i]}##*/}"
    done
    echo ""
    read _rmbl_i
    _rmbl_exec_i=$_rmbl_exec[$(($_rmbl_i + 1))]
    sudo wget -nc $_rmbl_exec_i
    _rmbl_exec_i=$(find . -name rmblast-*.tar.gz -print)
    sudo tar zxkf $_rmbl_exec_i && sudo rm -f $_rmbl_exec_i

    cd $INSTALL_DIR
}

function install_blat () {
    # @BLAT
    # Install Jim Kent's BLAST-like alignment tool

    INSTALL_DIR=$(pwd)
    _blat_url=http://hgwdev.cse.ucsc.edu/~kent/src
    _blat=blatSrc.zip
    _blatdir=$_srcdir/${_blat%%.*}

    cd $_srcdir
    wget -nc $_blat_url/$_blat
    unzip $_blat && rm -f $_blat
    cd $_blatdir

    ## if the default MACHTYPE is long with dashes, e.g. x86_64-pc-linux-gnu,
    ## set MACHTYPE to it's correct short value
    export MACHTYPE=${MACHTYPE:+$(uname -m)}
    [ -z "$MACHTYPE" ] && { echo "MACHTYPE not set"; exit 1 }

    ## make a temporary bin/$MACHTYPE dir in $HOME
    _tmpbin="$HOME/bin/$MACHTYPE"
    [ -d "$_tmpbin" ] || mkdir -p $_tmpbin

    _blatbin="$_blatdir/bin"
    [ -d "$_blatbin" ] || mkdir -p $_blatbin

    _blatlib="${_blatdir}/lib/$MACHTYPE"
    [ -d "$_blatlib" ] || mkdir -p $_blatlib

    make 2>&1 |tee $_blatdir/blat-make.log

    ## This will move the executables into _TMPBIN
    ## Relocate to _BLATBIN and link to _BINDIR
    mv $_tmpbin/* $_blatbin
    for x in $_blatbin/*; do
        sudo ln -sf $x $_bindir/${x##*/}
    done

    ## cleanup $_tmpbin
    rm -rf $_tmpbin
    cd $INSTALL_DIR
}

function install_mercator () {
    # @mercator
    # Install Colin Dewey's Source Code distribution
    # which includes mercator and a great number of other tools

    INSTALL_DIR=$(pwd)
    _cndsrc_url=http://www.biostat.wisc.edu/~cdewey/software
    _cndsrc=cndsrc-2010.10.11.tar.gz
    _cnddir=$_srcdir/cndsrc-2010.10.11

    cd $_srcdir
    wget -nc $_cndsrc_url/$_cndsrc
    tar --no-same-owner --no-same-permissions -xzkf $_cndsrc -C $_srcdir && rm -f $_cndsrc
    cd $_cnddir

    ## test for dependencies
    x=$(which blat)
    if [ ! $? ]; then
        echo ""
        echo "Also going to install the BLAST-like alignment tool - BLAT ..."
        echo ""
        sleep 2
        install_blat
    fi

    ## test version of g++; 4.7 does currently not compile
    _gpp_version=$(g++ --version | awk -F ' +' '/g\+\+/ {print $(NF)}' | awk -F.  '{print $1 "." $2}')
    if [ "$_gpp_version" != 4.6 ]; then
        printf "g++ is version %s\n" $_gpp_version
        echo "g++ version 4.6 required for mercator to compile"
        exit 1
    fi

    make 2>&1 |tee $_cnddir/mercator-make.log
    make install prefix=$_cnddir
    for x in $_cnddir/bin/*; do
        sudo ln -sf $x $_bindir/${x##*/}
    done
    cd $INSTALL_DIR
}

function install_fsa () {
    # @FSA
    # install Fast Statistical Alignment
    # and optionally the dependencies MUMmer, exonerate and mercator

    INSTALL_DIR=$(pwd)
    _fsa_url=http://sourceforge.net/projects/fsa/files/latest/download
    _fsa=fsa_latest.tar.gz

    echo ""
    echo "For long sequence alignments FSA needs MUMmer for 'anchor anealing'"
    echo "Shall we install MUMmer (y/n) [y]: "
    read _with_mummer | tr -d "\n"
    _with_mummer=${_with_mummer:-y}

    echo ""
    echo "FSA can also use the programm exonerate to detect remote homology"
    echo "Shall we install exonerate (y/n) [y]: "
    read _with_exonerate | tr -d "\n"
    _with_exonerate=${_with_exonerate:-y}

    echo ""
    echo "If we want to align genomes with rearrangements we need first"
    echo "construct a homology map with mercator and then run FSA on the"
    echo "homologous segments."
    echo "Shall we install mercator (y/n) [y]: "
    read _with_mercator | tr -d "\n"
    _with_mercator=${_with_mercator:-y}

    ## dependencies mummer, exonerate and mercator
    if [ ! which mummer 2>&1 >/dev/null && "$_with_mummer" = "y" ]; then
        sudo apt-get --yes install mummer
    fi

    if [ ! which exonerate 2>&1 >/dev/null && "$_with_exonerate" = "y" ]; then
        sudo apt-get --yes install exonerate
    fi

    if [ ! which mercator 2>&1 >/dev/null && "$_with_mercator" = "y" ]; then
        install_mercator
    fi

    ## get the source
    cd $_srcdir
    if [ ! -e "$_fsa"  ]; then
        wget $_fsa_url -O $_fsa
    fi
    tar --no-same-owner --no-same-permissions -xzkf $_fsa -C $_srcdir && rm -f $_fsa
    _fsadir=$(find $_srcdir -type d -name 'fsa-*' -print)
    cd $_fsadir

    case $_with_mummer in 
        y) _mummer=--with-mummer ;;
        *) _mummer= ;;
    esac

    case $_with_exonerate in
        y) _exonerate=--with-exonerate ;;
        *) _exonerate= ;;
    esac

    ./configure --prefix=${_bindir%/bin} \
                $_mummer \
                $_exonerate \
                --disable-gui 2>&1 |tee $_fsadir/fsa-config.log

    ## we have to supply -XDignore.symbol.file to javac to compile the gui
    for makefile in ./display/Makefile.*; do
        sed 's/JAVACFLAGS = /&-XDignore.symbol.file /' < ${makefile} > ${makefile}~
        mv -f ${makefile}~ ${makefile}
    done

    make 2>&1 |tee $_fsadir/fsa-make.log
    sudo make install
    cd $INSTALL_DIR
}

function install_glimmer () {
    # @glimmer 
    # install Glimmer3 and ELPH

    INSTALL_DIR=$(pwd)
    _glimmer_url=http://www.cbcb.umd.edu/software/glimmer
    _glimmer=glimmer302a.tar.gz
    _glimmerdir=$_srcdir/glimmer3.02

    _elph_url=ftp://ftp.cbcb.umd.edu/pub/software/elph
    _elph=ELPH-1.0.1.tar.gz
    _elphdir=$_srcdir/ELPH

    echo ""
    echo "Also going to install the dependency ELPH ..."
    echo ""
    sleep 2
    cd $_srcdir

    wget -nc $_glimmer_url/$_glimmer
    tar --no-same-owner --no-same-permissions -zxkf $_glimmer -C $_srcdir && rm -f $_glimmer

    wget -nc $_elph_url/$_elph
    tar --no-same-owner --no-same-permissions -zxkf $_elph -C $_srcdir && rm -f $_elph
	chmod a+rx $_elphdir
	chmod a+rx $_elphdir/sources

    cd $_glimmerdir/SimpleMake
    make 2>&1 |tee ../glimmer-make.log
    for x in $_glimmerdir/bin/*; do
        sudo ln -sf $x $_bindir/${x##*/}
    done

    cd $_elphdir/sources
    make 2>&1 |tee ../elph-make.log
    sudo ln -sf $_elphdir/sources/elph $_bindir/elph
    cd $INSTALL_DIR
}


function install_repeatmasker () {
    # @RepeatMasker
    # Install RepeatMasker and stuff

    INSTALL_DIR=$(pwd)
    echo ""
    echo "Also going to install RMBlast, Tandem Repeat Finder and"
    echo "RepeatMasker Libraries"
    echo ""
    sleep 2

    _trf_url=https://s3-eu-west-1.amazonaws.com/trf64linux/trf407b.linux64
    _trfdir=$_srcdir/trf

    _rm_url=http://www.repeatmasker.org
    _rm=RepeatMasker-open-3-3-0-p1.tar.gz
    _rmdir=$_srcdir/RepeatMasker

    _rm_lib_url=https://s3-eu-west-1.amazonaws.com/repmasklib/repeatmaskerlibraries-20120418.tar.gz
    _rm_lib=repeatmaskerlibraries-20120418.tar.gz

    ## PREREQUISITES
    ## 1) Blast+ and RMBlast
    if ! command -v /usr/local/rmblast/bin/blastn >/dev/null 2>&1; then
        select_blast
        cd /usr/local/
        _rmbldir=$(find . -type d -name rmblast-* -print)
        _blastdir=$(find . -type d -name ncbi-blast-* -print)
        sudo cp -R $_rmbldir/* $_blastdir
        sudo rm -rf $_rmbldir
        sudo mv $_blastdir rmblast
    fi

    ## 2) Tandem Repeat Finder
    if ! command -v trf >/dev/null 2>&1; then
        [ -d "$_trfdir" ] || mkdir -p $_trfdir
        cd $_trfdir
        wget -nc $_trf_url
        chmod +x trf407b.linux64 && sudo ln -sf $_trfdir/trf407b.linux64 $_bindir/trf
    fi

    ## Download Repeatmasker
    cd $_srcdir
    wget -nc $_rm_url/$_rm
    tar zxkf $_rm -C $_srcdir && rm -f $_rm
    cd $_rmdir

    ## Download and unpack the Repeatmasker libraries
    wget -nc $_rm_lib_url
    tar -xvf $_rm_lib #&& rm -f $_rm_lib

    ## Configure RepeatMasker
    _perl=$(which perl)
    for x in RepeatMasker ProcessRepeats DateRepeats; do
        sed '1 s_.*_#!'"${_perl}"'_' < $x > ${x}~ && mv -f ${x}~ $x
        chmod +x $x
    done

    perl configure
    sudo ln -sf $_rmdir/RepeatMasker $_bindir/RepeatMasker
    cd $INSTALL_DIR
}

function install_repeatscout () {
    # @RepeatScout

    INSTALL_DIR=$(pwd)
    echo ""
    echo "Also going to install nseq and trf"
    echo ""
    sleep 2

    _trf_url=https://s3-eu-west-1.amazonaws.com/trf64linux/trf407b.linux64
    _trfdir=$_srcdir/trf

    _nseg_url=ftp://ftp.ncbi.nih.gov/pub/seg/nseg
    _nsegdir=$_srcdir/nseg

    _rs_url=http://bix.ucsd.edu/repeatscout
    _rs=RepeatScout-1.0.5.tar.gz
    _rsdir=$_srcdir/RepeatScout-1

    ## 1) Tandem Repeat Finder
    if ! command -v trf >/dev/null 2>&1; then
        [ -d "$_trfdir" ] || mkdir -p $_trfdir
        cd $_trfdir
        wget -nc $_trf_url
        chmod +x trf407b.linux64 && sudo ln -sf $_trfdir/trf407b.linux64 $_bindir/trf
    fi

    ## 2) Nseg
    if ! command -v nseg >/dev/null 2>&1; then
        [ -d "$_nsegdir" ] || mkdir -p $_nsegdir
        cd $_srcdir
        wget -rnH --cut-dirs=2 $_nseg_url
        cd $_nsegdir
        make
        sudo ln -s $_nsegdir/nmerge $_bindir/nmerge
        sudo ln -s $_nsegdir/nseg $_bindir/nseg
    fi

    ## Download RepeatScout
    cd $_srcdir
    wget -nc $_rs_url/$_rs
    tar zxkf $_rs -C $_srcdir && rm -f $_rs
    cd $_rsdir && make

    for x in ${$(find . -type f -executable)##*/}; do
        sudo ln -s $_rsdir/$x $_bindir/$x
    done

    cd $INSTALL_DIR
}

# @MENU
clear
while true
do
    echo ""
    echo " Install software: "
    echo " [1] BLAT : The BLAST-like alignment tool"
    echo " [2] Mercator : Whole-genome orthology maps (installs also BLAT)"
    echo " [3] FSA : Fast statistical alignment (installs also Mercator and Blat)"
    echo " [4] Glimmer3 : Predict prokaryotic genes (installs also ELPH)"
    echo " [5] RepeatMasker : Mask interspersed repeats and low-complexity regions"
    echo " [6] RepeatScout : Identify repeat family sequences"
    echo " [q] Quit"
    echo ""
    echo " Choose a program for installation: " | tr -d "\n"

    read program

    case "$program" in
        "1" )
            set_env
            install_blat ;;
        "2" )
            set_env
            install_mercator ;;
        "3" )
            set_env
            install_fsa ;;
        "4" )
            set_env
            install_glimmer ;;
        "5" )
            set_env
            install_repeatmasker ;;
        "6" )
            set_env
            install_repeatscout ;;
        q|Q )
            break ;;
         * )
            echo "Illegal choice"
            continue ;;
    esac
done

exit 0
