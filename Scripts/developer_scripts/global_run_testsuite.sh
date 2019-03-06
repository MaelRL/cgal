#!/bin/bash

#!/bin/bash
#usage : script [-c -l -n -s -k] testsuite_dir

##########################
####   LAUNCH CTEST   ####
##########################

export SHOW_PROGRESS=""
export KEEP_TESTS=""
export DO_NOT_UPLOAD=""
export DO_NOT_TEST=""
export SCP="scp"
export WGET="wget"
export WGET_OPTS="--no-check-certificate --no-verbose"
export CURL="curl"
export CURL_OPTS="-k --remote-name --silent --location-trusted"
export CGAL_URL="https://cgal.geometryfactory.com/CGAL/Members/Releases"
export UPLOAD_RESULT_DESTINATION="cgaltest@cgaltest.geometryfactory.com:incoming"
export LATEST_LOCATION="${CGAL_URL}/LATEST"
export TAR="tar"
export GZIP="gzip"
export GUNZIP="gunzip"
export COMPRESSOR="${GZIP}"
export CONSOLE_OUTPUT="y"
export CGAL_ROOT=`pwd`
export USE_TARGZ="n"
export USE_TARBZ="n"
export CGAL_RELEASE=""
export LOGS_DIR=""
export LOCK_FILE=""
export LIST_TEST_PACKAGES=""
export ACTUAL_LOGFILE=""
export CGAL_DIR=""
USE_LATEST_UNZIPPED=""
WITH_DOCKER="n"

# ----------------------------------------------------------------------------------------
# Downloads the file "LATEST" whose contents indicates which release to test
# ----------------------------------------------------------------------------------------
download_latest()
{
  if [ -r "LATEST" ]; then
    rm -rf LATEST
  fi
    log "${ACTUAL_LOGFILE}" "getting LATEST"
    if [ -n "${USE_CURL}" ]; then
      ${CURL} ${CURL_OPTS} "${LATEST_LOCATION}" >> "${ACTUAL_LOGFILE}" 2>&1
    else
      ${WGET} ${WGET_OPTS} "${LATEST_LOCATION}" >> "${ACTUAL_LOGFILE}" 2>&1
    fi
  if [ ! -f "LATEST" ]; then
    error "COULD NOT DOWNLOAD LATEST!"
  fi
}

# ----------------------------------------------------------------------------------------
# Exits the testsuite if the latest release has been already tested.
# This is tested by comparing files LATEST and RELEASE_NR, where
# RELEASE_NR is a copy of the previous LATEST.
# ----------------------------------------------------------------------------------------
abort_if_latest_already_tested()
{
  if [ -r "RELEASE_NR" ]; then
    cmp LATEST RELEASE_NR >> "${ACTUAL_LOGFILE}"
    if [ ! ${?} != 0 ]; then
      log "${ACTUAL_LOGFILE}" "This release has already been tested."
      rm -f "$LOCK_FILE";
      exit 1;
    fi
  fi
}



# ----------------------------------------------------------------------------------------
# get CGAL
# ----------------------------------------------------------------------------------------
get_cgal()
{
  if [ -z "$CGAL_LOCATION" ]; then
    for i in `cat LATEST`
    do
      CGAL_LOCATION="${CGAL_URL}/${i}";
      CGAL_ZIPFILE="${i}";
    done
  else
    CGAL_ZIPFILE=`echo "$CGAL_LOCATION" | sed 's|.*/||'`
  fi

  CGAL_RELEASE_ID=`echo $CGAL_ZIPFILE | sed "s/.tar.gz//"`
  if [ ! "${CGAL_RELEASE_ID}" = "${CGAL_ZIPFILE}" ]; then
    USE_TARGZ="y"
  else
    CGAL_RELEASE_ID=`echo $CGAL_ZIPFILE | sed "s/.tar.bz2//"`
    if [ ! "${CGAL_RELEASE_ID}" = "${CGAL_ZIPFILE}" ]; then
      USE_TARBZ="y"
    fi
  fi

  log "${ACTUAL_LOGFILE}" "CGAL_ZIPFILE    = ${CGAL_ZIPFILE}"
  log "${ACTUAL_LOGFILE}" "CGAL_RELEASE_ID = ${CGAL_RELEASE_ID}"

  log "${ACTUAL_LOGFILE}" "getting CGAL"
  rm -f "${CGAL_ZIPFILE}"
  if [ -n "${USE_CURL}" ]; then
      ${CURL} ${CURL_OPTS} "${CGAL_LOCATION}" >> "${ACTUAL_LOGFILE}" 2>&1
  else
      ${WGET} ${WGET_OPTS} "${CGAL_LOCATION}" >> "${ACTUAL_LOGFILE}" 2>&1
  fi
  if [ ${?} != 0 ]; then
    error "Could not get CGAL"
  fi
  log_done "${ACTUAL_LOGFILE}"
}


# ----------------------------------------------------------------------------------------
# Unzips and untars the downloaded CGAL release
# ----------------------------------------------------------------------------------------
unzip_cgal()
{
  cd "${CGAL_ROOT}"

  log "${ACTUAL_LOGFILE}" "unzipping CGAL"
  if [ "${USE_TARGZ}" = "y" ]; then
    DECOMPRESSOR="${GUNZIP}"
    log_done "${ACTUAL_LOGFILE}"
  fi

  if [ "${USE_TARBZ}" = "y" ]; then
    DECOMPRESSOR="bunzip2"
  fi

  log "${ACTUAL_LOGFILE}" "untarring CGAL"
  ${DECOMPRESSOR} -c "${CGAL_ZIPFILE}" | ${TAR} xf - >> "${ACTUAL_LOGFILE}" 2>&1
  if [ ${?} != 0 ]; then
    error "Could not untar CGAL"
  fi

  # check, if CGAL_DIR exists
  if [ -d "${CGAL_ROOT}/${CGAL_RELEASE_ID}" ]; then
    # Reset CGAL-I symlink
    log "${ACTUAL_LOGFILE}" "Resetting CGAL-I symlink to ${CGAL_ROOT}/${CGAL_RELEASE_ID}"
    rm -f CGAL-I
    ln -s "${CGAL_ROOT}/${CGAL_RELEASE_ID}" CGAL-I
    # Reset CGAL-3.x-I symlink
    CGAL_RELEASE=`echo "${CGAL_RELEASE_ID}" | sed 's/I\([^-]*\)-.*/I\1/'`
    log "${ACTUAL_LOGFILE}" "Resetting ${CGAL_RELEASE} symlink to ${CGAL_ROOT}/${CGAL_RELEASE_ID}"
    rm -f "${CGAL_RELEASE}"
    ln -s "${CGAL_ROOT}/${CGAL_RELEASE_ID}" "${CGAL_RELEASE}"
  else
    error "directory ${CGAL_ROOT}/${CGAL_RELEASE_ID} does not exist"
  fi

  log_done "${ACTUAL_LOGFILE}"
}



# Parse command line arguments
for arg in "$@"
do
  case "$arg" in 
    "-c")
      echo "Using latest unzipped release instead of getting a new one from the server"
      USE_LATEST_UNZIPPED="y"
      ;;
    "-l")
      echo "Not uploading results to dashboard"
      DO_NOT_UPLOAD="y"
      ;;
    "-n")
       echo "No testsuite will be launched. Installation only."
      DO_NOT_TEST="y"
      ;;
    "-s")
      echo "Showing progress."
      SHOW_PROGRESS="y"
      ;;
    "-k")
      echo "Compiled test/ directory will be kept."
      KEEP_TESTS="y"
      ;;
    "-d")
      echo "The testsuite will be launched through docker"
      WITH_DOCKER="y"
      ;;
    *)
      CGAL_LOCATION=$arg
  esac
done

# Load settings
if [ -f "$CGAL_ROOT/.autocgalrc" ]; then
  . "$CGAL_ROOT/.autocgalrc"
else
  echo "CONFIGURATION FILE  .autocgalrc NOT FOUND" >&2;
  exit 1
fi


CGAL_DIR=`readlink "${CGAL_ROOT}/CGAL-I"`
source "${CGAL_DIR}/developer_scripts/log.sh"
LOGS_DIR="${CGAL_ROOT}/AUTOTEST_LOGS"
LOCK_FILE="${CGAL_ROOT}/autotest_cgal_with_cmake.lock"
LIST_TEST_PACKAGES="${CGAL_ROOT}/list_test_packages"

# Setup logfile
ACTUAL_LOGFILE="${CGAL_ROOT}/`basename ${0}`.log"
rm -f "${ACTUAL_LOGFILE}"

echo "Running `basename ${0}` "'$Revision$' >> "${ACTUAL_LOGFILE}" 

cd "$CGAL_ROOT"

# Starts the process

if [ -z "${USE_LATEST_UNZIPPED}" ]; then
  if [ -z "$CGAL_LOCATION" ]; then
    download_latest
    abort_if_latest_already_tested
  fi
  get_cgal
  unzip_cgal
fi

if [ "${WITH_DOCKER}" = "y" ]; then
  #launch docker container
  echo "export SHOW_PROGRESS=$SHOW_PROGRESS"> env.sh
  echo "export KEEP_TESTS=$KEEP_TESTS">> env.sh
  echo "export DO_NOT_UPLOAD=$DO_NOT_UPLOAD">> env.sh
  echo "export DO_NOT_TEST=$DO_NOT_TEST">> env.sh
  echo "export SCP=$SCP">> env.sh
  echo "export WGET=$WGET">> env.sh
  echo "export WGET_OPTS=\"$WGET_OPTS\"">> env.sh
  echo "export CURL=$CURL">> env.sh
  echo "export CURL_OPTS=\"$CURL_OPTS\"">> env.sh
  echo "export CGAL_URL=$CGAL_URL">> env.sh
  echo "export UPLOAD_RESULT_DESTINATION=$UPLOAD_RESULT_DESTINATION">> env.sh
  echo "export LATEST_LOCATION=$LATEST_LOCATION">> env.sh
  echo "export TAR=$TAR">> env.sh
  echo "export GZIP=$GZIP">> env.sh
  echo "export GUNZIP=$GUNZIP">> env.sh
  echo "export COMPRESSOR=$COMPRESSOR">> env.sh
  echo "export CONSOLE_OUTPUT=$CONSOLE_OUTPUT">> env.sh
  echo "export CGAL_ROOT=/cgal_root">> env.sh
  echo "export USE_TARGZ=$USE_TARGZ">> env.sh
  echo "export USE_TARBZ=$USE_TARBZ">> env.sh
  echo "export CGAL_RELEASE=$CGAL_RELEASE">> env.sh
  echo "export CGAL_DIR=/cgal_root/\${CGAL_LAST}">>env.sh
  echo "export LIST_TEST_PACKAGES=/cgal_root/list_test_packages">> env.sh
  echo "export LOGS_DIR=/cgal_root/AUTOTEST_LOGS">> env.sh
  echo "export LOCK_FILE=/cgal_root/autotest_cgal_with_cmake.lock">> env.sh
  echo "export ACTUAL_LOGFILE=/cgal_root/\`basename \${0}\`.log">> env.sh
  docker pull cgal/testsuite-docker:debian-for-arm
  docker run --rm  -t -e CGAL_LAST=${CGAL_RELEASE_ID} -v ${CGAL_ROOT}/ssh:/tmp_ssh -v ${DEPS_DIR}:/deps -v ${CGAL_ROOT}/:/cgal_root cgal/testsuite-docker:debian-for-arm
else
  bash ${CGAL_DIR}/developer_scripts/run_testsuite_with_cmake
fi
