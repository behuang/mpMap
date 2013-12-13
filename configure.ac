AC_INIT([gputools], 0.01)
AC_CONFIG_AUX_DIR([tools])

AC_MSG_NOTICE([checking for critical tools used by configure])

AC_ARG_WITH([cuda-home],
	AC_HELP_STRING([--with-cuda-home=cuda_home],
		[specify the root of your cuda install, e.g. /usr/local/cuda]),
	[CUDA_HOME="$withval"])

AC_ARG_WITH([r-include],
	AC_HELP_STRING([--with-r-include=r_include],
		[specify the full path to your R header files, e.g. /usr/share/R/include]),
	[R_INCLUDE="$withval"])

AC_ARG_WITH([r-lib],
	AC_HELP_STRING([--with-r-lib=r-lib],
		[specify the full path to your R shared libraries, e.g. /usr/lib64/R/lib]),
	[R_LIB="$withval"])

AC_MSG_CHECKING([environment variable CUDA_HOME])
if test -z "${CUDA_HOME}"; then
	AC_MSG_RESULT([CUDA_HOME not set; using default /usr/local/cuda])
	CUDA_HOME="/usr/local/cuda"
else
	AC_MSG_RESULT([using CUDA_HOME=${CUDA_HOME}])
fi

HAS_CUDA=false
AC_CHECK_FILE([${CUDA_HOME}/bin/nvcc], [HAS_NVCC="yes"])
if test -z "${HAS_NVCC}"; then
	AC_MSG_WARN([cannot find nvcc; compiling without CUDA support])
	HAS_CUDA=false
else
	AC_MSG_RESULT([nvcc found])
	CUDA_OBJS=rfhaps_gpu.o
	HAS_CUDA=true

	AC_MSG_CHECKING([whether this is the 64 bit linux version of CUDA])
	AC_CHECK_FILE([${CUDA_HOME}/lib64/libcublas.so], [BIT64="yes"])
	if test -z "${BIT64}"; then
		AC_MSG_RESULT([no -- using ${CUDA_HOME}/lib for CUDA libs])
		CUDA_LIB_DIR="/lib"
	else
		AC_MSG_RESULT([yes -- using ${CUDA_HOME}/lib64 for CUDA libs])
		CUDA_LIB_DIR="/lib64"
	fi
fi

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
	AC_MSG_ERROR([could not determine R_HOME])
else
	CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
	if test -z "${R_INCLUDE}"; then
		CUDA_CPPFLAGS=" $("${R_HOME}/bin/R" CMD config --cppflags) "
		PKG_CPPFLAGS="$("${R_HOME}/bin/R" CMD config --cppflags) $("${R_HOME}/bin/Rscript" -e "Rcpp:::CxxFlags()") -openmp -DUSE_OPENMP "
	else
		CUDA_CPPFLAGS=" -I${R_INCLUDE} "
		PKG_CPPFLAGS="-I${R_INCLUDE} $("${R_HOME}/bin/Rscript" -e "Rcpp:::CxxFlags()") -openmp -DUSE_OPENMP "
	fi
	if test -z "${R_LIB}"; then
		R_LIB="`"${R_HOME}/bin/R" CMD config --ldflags` $("${R_HOME}/bin/Rscript" -e "Rcpp:::LdFlags()")"
	else
		R_LIB="-L${R_LIB} -lR $("${R_HOME}/bin/Rscript" -e "Rcpp:::LdFlags()")"
	fi
	AC_MSG_RESULT([using ${R_INCLUDE} for R header files])
	AC_MSG_RESULT([using ${R_LIB} for R shared libraries])
fi
DEVICEOPTS=""

if test -z "${HAS_NVCC}"; then
	:
else
	LIBS=" -L${CUDA_HOME}${CUDA_LIB_DIR} -lcudart "
	PKG_CPPFLAGS+=" -DHAS_CUDA "
fi


AC_MSG_NOTICE([determining host platform])
arch=`uname -m`
osrelease=`uname -r`
osversion=`uname -v`

RPATHFLAG=""
AC_CANONICAL_HOST
case "${host_os}" in
	*darwin*)
		RPATHFLAG=" -rpath -Xlinker ${CUDA_HOME}${CUDA_LIB_DIR} "
		LIBS+=" -lstdc++ "
		AC_PATH_PROG([SW_VERS],[sw_vers],[no],[$PATH])
		if test "$arch" = "i386"; then
			$SW_VERS -productVersion | grep '^10.6' > /dev/null 2>&1
			if test "x$?" = "x0" ; then
				DEVICEOPTS="-m64"
			fi
		elif test "$arch" = "x86_64"; then
			DEVICEOPTS="-m64"
		else
			AC_MSG_WARN([Unrecognised Darwin Architecture])
		fi
	;;
esac

AC_SUBST(CUDA_HOME)
AC_SUBST(CUDA_OBJS)

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(R_LIB)
AC_SUBST(LIBS)

AC_SUBST(CPICFLAGS)
AC_SUBST(DEVICEOPTS)
AC_SUBST(RPATHFLAG)
AC_SUBST(HAS_CUDA)
AC_SUBST(CUDA_CPPFLAGS)

AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT
