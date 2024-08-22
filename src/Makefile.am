AM_CPPFLAGS = \
    -I${top_srcdir}/phyml \
    -I${top_srcdir}/seqgen \
    -I${top_srcdir}/levmar \
    -I${top_srcdir}/librna

bin_PROGRAMS = RNAcode

RNAcode_SOURCES = \
    RNAcode.c \
    RNAcode.h \
    rnaz_utils.h \
    rnaz_utils.c \
    utils.h \
    utils.c \
    cmdline.c \
    cmdline.h \
    treeSimulate.c \
    treeSimulate.h \
    treeML.h \
    treeML.c \
    code.c \
    postscript.c \
    code.h \
    extreme_fit.c \
    extreme_fit.h \
    postscript.h \
    score.c \
    score.h \
    misc.c \
    misc.h

LDADD = \
    -L${top_builddir}/levmar -llevmar \
    -L${top_builddir}/phyml -lPHYML \
    -L${top_builddir}/seqgen -lseqgen \
    -L${top_builddir}/librna -lRNA \
    -lm
