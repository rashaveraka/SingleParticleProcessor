set(root_build_options)

#---------------------------------------------------------------------------------------------------
#---ROOT_BUILD_OPTION( name defvalue [description] )
#---------------------------------------------------------------------------------------------------
function(ROOT_BUILD_OPTION opt defvalue)
  if(ARGN)
    set(description ${ARGN})
  else()
    set(description " ")
  endif()
  set(${opt}_defvalue    ${defvalue} PARENT_SCOPE)
  set(${opt}_description ${description} PARENT_SCOPE)
  set(root_build_options  ${root_build_options} ${opt} PARENT_SCOPE )
endfunction()

#---------------------------------------------------------------------------------------------------
#---ROOT_APPLY_OPTIONS()
#---------------------------------------------------------------------------------------------------
function(ROOT_APPLY_OPTIONS)
  foreach(opt ${root_build_options})
     option(${opt} "${${opt}_description}" ${${opt}_defvalue})
  endforeach()  
endfunction()

#---------------------------------------------------------------------------------------------------
#---ROOT_SHOW_OPTIONS([var] )
#---------------------------------------------------------------------------------------------------
function(ROOT_SHOW_OPTIONS)
  set(enabled)
  foreach(opt ${root_build_options})
    if(${opt})
      set(enabled "${enabled} ${opt}")
    endif()
  endforeach()
  if(NOT ARGN)
    message(STATUS "Enabled support for: ${enabled}")
  else()
    set(${ARGN} "${enabled}" PARENT_SCOPE)
  endif()
endfunction()

#---------------------------------------------------------------------------------------------------
#---ROOT_WRITE_OPTIONS(file )
#---------------------------------------------------------------------------------------------------
function(ROOT_WRITE_OPTIONS file)
  file(WRITE ${file} "#---Options enabled for the build of ROOT-----------------------------------------------\n")
  foreach(opt ${root_build_options})
    if(${opt})
      file(APPEND ${file} "set(${opt} ON)\n")
    else()
      file(APPEND ${file} "set(${opt} OFF)\n")
    endif()
  endforeach()
endfunction()


#--------------------------------------------------------------------------------------------------
#---Full list of options with their descriptios and default values
#   The default value can be changed as many times as we wish before calling ROOT_APPLY_OPTIONS()
#--------------------------------------------------------------------------------------------------

ROOT_BUILD_OPTION(afdsmgrd OFF "Dataset manager for PROOF-based analysis facilities")
ROOT_BUILD_OPTION(afs OFF "AFS support, requires AFS libs and objects")
ROOT_BUILD_OPTION(alien ON "AliEn support, requires libgapiUI from ALICE")
ROOT_BUILD_OPTION(asimage ON "Image processing support, requires libAfterImage")
ROOT_BUILD_OPTION(astiff ON "Include tiff support in image processing")
ROOT_BUILD_OPTION(bonjour ON "Bonjour support, requires libdns_sd and/or Avahi")
ROOT_BUILD_OPTION(builtin_afterimage ON "Built included libAfterImage, or use system libAfterImage")
ROOT_BUILD_OPTION(builtin_fftw3 OFF "Built the FFTW3 library internally (downloading tarfile from the Web)")
ROOT_BUILD_OPTION(builtin_ftgl ON "Built included libFTGL, or use system libftgl")
ROOT_BUILD_OPTION(builtin_freetype OFF "Built included libfreetype, or use system libfreetype")
ROOT_BUILD_OPTION(builtin_glew ON "Built included libGLEW, or use system libGLEW")
ROOT_BUILD_OPTION(builtin_openssl OFF "Build OpenSSL internally, or use system OpenSSL")
ROOT_BUILD_OPTION(builtin_pcre OFF "Built included libpcre, or use system libpcre")
ROOT_BUILD_OPTION(builtin_zlib OFF "Built included libz, or use system libz")
ROOT_BUILD_OPTION(builtin_lzma OFF "Built included liblzma, or use system liblzma")
ROOT_BUILD_OPTION(builtin_davix OFF "Built the Davix library internally (downloading tarfile from the Web)")
ROOT_BUILD_OPTION(builtin_gsl OFF "Built the GSL library internally (downloading tarfile from the Web)")
ROOT_BUILD_OPTION(builtin_cfitsio OFF "Built the FITSIO library internally (downloading tarfile from the Web)")
ROOT_BUILD_OPTION(builtin_xrootd OFF "Built the XROOTD internally (downloading tarfile from the Web)")
ROOT_BUILD_OPTION(builtin_llvm ON "Built the LLVM internally")
ROOT_BUILD_OPTION(builtin_tbb OFF "Built the TBB internally")
ROOT_BUILD_OPTION(cxx11 ON "Build using C++11 compatible mode, requires gcc > 4.7.x or clang")
ROOT_BUILD_OPTION(cxx14 OFF "Build using C++14 compatible mode, requires gcc > 4.9.x or clang")
ROOT_BUILD_OPTION(libcxx OFF "Build using libc++, requires cxx11 option (MacOS X only, for the time being)")
ROOT_BUILD_OPTION(castor ON "CASTOR support, requires libshift from CASTOR >= 1.5.2")
ROOT_BUILD_OPTION(ccache OFF "Enable ccache usage for speeding up builds")
ROOT_BUILD_OPTION(chirp ON "Chirp support (Condor remote I/O), requires libchirp_client")
ROOT_BUILD_OPTION(cling ON "Enable new CLING C++ interpreter")
ROOT_BUILD_OPTION(cocoa OFF "Use native Cocoa/Quartz graphics backend (MacOS X only)")
ROOT_BUILD_OPTION(davix ON "DavIx library for HTTP/WEBDAV access")
ROOT_BUILD_OPTION(dcache ON "dCache support, requires libdcap from DESY")
ROOT_BUILD_OPTION(exceptions ON "Turn on compiler exception handling capability")
ROOT_BUILD_OPTION(explicitlink ON "Explicitly link with all dependent libraries")
ROOT_BUILD_OPTION(fftw3 ON "Fast Fourier Transform support, requires libfftw3")
ROOT_BUILD_OPTION(fitsio ON "Read images and data from FITS files, requires cfitsio")
ROOT_BUILD_OPTION(fortran ON "Enable the Fortran components of ROOT")
set(gcctoolchain "" CACHE PATH "Path for the gcctoolchain in case not the system gcc is used to build clang/LLVM")
ROOT_BUILD_OPTION(gviz ON "Graphs visualization support, requires graphviz")
ROOT_BUILD_OPTION(gdml OFF "GDML writer and reader")
ROOT_BUILD_OPTION(geocad OFF "ROOT-CAD Interface")
ROOT_BUILD_OPTION(genvector ON "Build the new libGenVector library")
ROOT_BUILD_OPTION(gfal ON "GFAL support, requires libgfal")
ROOT_BUILD_OPTION(glite ON "gLite support, requires libglite-api-wrapper v.3 from GSI (https://subversion.gsi.de/trac/dgrid/wiki)")
ROOT_BUILD_OPTION(globus OFF "Globus authentication support, requires Globus toolkit")
ROOT_BUILD_OPTION(gnuinstall OFF "Perform installation following the GNU guidelines")
ROOT_BUILD_OPTION(gsl_shared OFF "Enable linking against shared libraries for GSL (default no)")
ROOT_BUILD_OPTION(hdfs ON "HDFS support; requires libhdfs from HDFS >= 0.19.1")
ROOT_BUILD_OPTION(http OFF "HTTP Server support")
ROOT_BUILD_OPTION(jemalloc OFF "Using the jemalloc allocator")
ROOT_BUILD_OPTION(krb5 ON "Kerberos5 support, requires Kerberos libs")
ROOT_BUILD_OPTION(ldap ON "LDAP support, requires (Open)LDAP libs")
ROOT_BUILD_OPTION(mathmore ON "Build the new libMathMore extended math library, requires GSL (vers. >= 1.8)")
ROOT_BUILD_OPTION(memstat ON "A memory statistics utility, helps to detect memory leaks")
ROOT_BUILD_OPTION(minuit2 OFF "Build the new libMinuit2 minimizer library")
ROOT_BUILD_OPTION(monalisa ON "Monalisa monitoring support, requires libapmoncpp")
ROOT_BUILD_OPTION(mt OFF "Multi-threading support")
ROOT_BUILD_OPTION(mysql ON "MySQL support, requires libmysqlclient")
ROOT_BUILD_OPTION(odbc ON "ODBC support, requires libiodbc or libodbc")
ROOT_BUILD_OPTION(opengl ON "OpenGL support, requires libGL and libGLU")
ROOT_BUILD_OPTION(oracle ON "Oracle support, requires libocci")
ROOT_BUILD_OPTION(pch ON)
ROOT_BUILD_OPTION(pgsql ON "PostgreSQL support, requires libpq")
ROOT_BUILD_OPTION(pythia6 ON "Pythia6 EG support, requires libPythia6")
ROOT_BUILD_OPTION(pythia6_nolink OFF "Delayed linking of Pythia6 library")
ROOT_BUILD_OPTION(pythia8 ON "Pythia8 EG support, requires libPythia8")
ROOT_BUILD_OPTION(python ON "Python ROOT bindings, requires python >= 2.2")
ROOT_BUILD_OPTION(qt OFF "Qt graphics backend, requires libqt >= 4.8")
ROOT_BUILD_OPTION(qtgsi OFF "GSI's Qt integration, requires libqt >= 4.8")
ROOT_BUILD_OPTION(roofit OFF "Build the libRooFit advanced fitting package")
ROOT_BUILD_OPTION(root7 OFF "Build the ROOT 7 interface prototype, requires cxx14")
ROOT_BUILD_OPTION(ruby OFF "Ruby ROOT bindings, requires ruby >= 1.8")
ROOT_BUILD_OPTION(r OFF "R ROOT bindings, requires R, Rcpp and RInside")
ROOT_BUILD_OPTION(rfio ON "RFIO support, requires libshift from CASTOR >= 1.5.2")
ROOT_BUILD_OPTION(rpath OFF "Set run-time library load path on executables and shared libraries (at installation area)")
ROOT_BUILD_OPTION(sapdb ON "MaxDB/SapDB support, requires libsqlod and libsqlrte")
ROOT_BUILD_OPTION(shadowpw ON "Shadow password support")
ROOT_BUILD_OPTION(shared ON "Use shared 3rd party libraries if possible")
ROOT_BUILD_OPTION(soversion OFF "Set version number in sonames (recommended)")
ROOT_BUILD_OPTION(sqlite ON "SQLite support, requires libsqlite3")
ROOT_BUILD_OPTION(srp ON "SRP support, requires SRP source tree")
ROOT_BUILD_OPTION(ssl ON "SSL encryption support, requires openssl")
ROOT_BUILD_OPTION(tbb OFF "TBB multi-threading support, requires TBB")
ROOT_BUILD_OPTION(table OFF "Build libTable contrib library")
ROOT_BUILD_OPTION(tcmalloc OFF "Using the tcmalloc allocator")
ROOT_BUILD_OPTION(thread ON "Using thread library (cannot be disabled)")
ROOT_BUILD_OPTION(tmva ON "Build TMVA multi variate analysis library")
ROOT_BUILD_OPTION(unuran OFF "UNURAN - package for generating non-uniform random numbers")
ROOT_BUILD_OPTION(vc OFF "Vc adds a few new types for portable and intuitive SIMD programming")
ROOT_BUILD_OPTION(vdt ON "VDT adds a set of fast and vectorisable mathematical functions")
ROOT_BUILD_OPTION(winrtdebug OFF "Link against the Windows debug runtime library")
ROOT_BUILD_OPTION(xft ON "Xft support (X11 antialiased fonts)")
ROOT_BUILD_OPTION(xml ON "XML parser interface")
ROOT_BUILD_OPTION(x11 ON "X11 support")
ROOT_BUILD_OPTION(xrootd ON "Build xrootd file server and its client (if supported)")

option(fail-on-missing "Fail the configure step if a required external package is missing" OFF)
option(minimal "Do not automatically search for support libraries" OFF)
option(gminimal "Do not automatically search for support libraries, but include X11" OFF)
option(all "Enable all optional components" OFF)
option(testing "Enable testing with CTest" OFF)
option(roottest "Include roottest, if roottest exists in root or if it is a sibling directory." OFF)

#--- Minor chnages in defaults due to platform--------------------------------------------------
if(WIN32)
  set(x11_defvalue OFF)
  set(memstat_defvalue OFF)
  set(davix_defvalue OFF)
elseif(APPLE)
  set(x11_defvalue OFF)
  set(cocoa_defvalue ON)
  set(davix_defvalue OFF)
endif()

#--- The 'all' option swithes ON major options---------------------------------------------------
if(all)
 set(gdml_defvalue ON)
 set(http_defvalue ON)
 set(qt_defvalue ON)
 set(qtgsi_defvalue ON)
 set(roofit_defvalue ON)
 set(minuit2_defvalue ON)
 set(r_defvalue ON)
 set(root7_defvalue ON)
 set(table_defvalue ON)
 set(unuran_defvalue ON)
 if(NOT APPLE)
   set(vc_defvalue ON)
 endif()
endif()

#---VC does not support yet Arm and PPC processors----------------------------------------------
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "aarch64" OR CMAKE_SYSTEM_PROCESSOR STREQUAL "ppc64le")
   message(STATUS "A system not supported by Vc, ${CMAKE_SYSTEM_PROCESSOR}, was detected. Disabling Vc by default.")
   set(vc_defvalue OFF)
endif()

#---Options depending of CMake Generator-------------------------------------------------------
if( CMAKE_GENERATOR STREQUAL Ninja)
   set(fortran_defvalue OFF)
endif()

#---Apply minimal or gminimal------------------------------------------------------------------
foreach(opt ${root_build_options})
  if(NOT opt MATCHES "thread|cxx11|cling|builtin_llvm|builtin_ftgl|explicitlink")
    if(minimal)
      set(${opt}_defvalue OFF)
    elseif(gminimal AND NOT opt MATCHES "x11|cocoa")
      set(${opt}_defvalue OFF)
    endif()
  endif()
endforeach()

#---Define at moment the options with the selected default values-----------------------------
ROOT_APPLY_OPTIONS()

#---Avoid creating dependencies to 'non-standard' header files -------------------------------
include_regular_expression("^[^.]+$|[.]h$|[.]icc$|[.]hxx$|[.]hpp$")

#---Add Installation Variables------------------------------------------------------------------
include(RootInstallDirs)

#---RPATH options-------------------------------------------------------------------------------
#  When building, don't use the install RPATH already (but later on when installing)
set(CMAKE_SKIP_BUILD_RPATH FALSE)         # don't skip the full RPATH for the build tree
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) # use always the build RPATH for the build tree
set(CMAKE_MACOSX_RPATH TRUE)              # use RPATH for MacOSX
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) # point to directories outside the build tree to the install RPATH

# Check whether to add RPATH to the installation (the build tree always has the RPATH enabled)
if(rpath OR gnuinstall)
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR}) # install LIBDIR
  set(CMAKE_SKIP_INSTALL_RPATH FALSE)          # don't skip the full RPATH for the install tree
elseif(APPLE)
  set(CMAKE_INSTALL_NAME_DIR "@rpath")
  set(CMAKE_INSTALL_RPATH "@loader_path/../lib")    # self relative LIBDIR
  set(CMAKE_SKIP_INSTALL_RPATH FALSE)          # don't skip the full RPATH for the install tree
else()
  set(CMAKE_SKIP_INSTALL_RPATH TRUE)           # skip the full RPATH for the install tree
endif()


