# -----------------------------------------------------------------------------
# Insert PyCA version/commit info into
# Code/Python/Common/PyCAVersion.py.in (creating PyCAVersion.py in binary dir)
# PyCA_SOURCE_DIR and PyCA_BINARY_DIR must be passed in to this script

# File Locations
SET(VERSION_FILE_IN ${PyCA_SOURCE_DIR}/Code/Python/Common/PyCAVersion.py.in)
SET(VERSION_FILE_OUT_DIR ${PyCA_BINARY_DIR}/Python/Common)
SET(VERSION_FILE_OUT ${VERSION_FILE_OUT_DIR}/PyCAVersion.py)

# Get the PyCA version as a string
EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} --git-dir ${PyCA_SOURCE_DIR}/.git --work-tree ${PyCA_SOURCE_DIR} describe --dirty
  OUTPUT_VARIABLE PYCA_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get the PyCA commit
EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} --git-dir ${PyCA_SOURCE_DIR}/.git --work-tree ${PyCA_SOURCE_DIR} rev-parse HEAD
  OUTPUT_VARIABLE PYCA_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get the current branch
EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} --git-dir ${PyCA_SOURCE_DIR}/.git --work-tree ${PyCA_SOURCE_DIR} rev-parse --abbrev-ref HEAD
  OUTPUT_VARIABLE PYCA_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE)
# make directory
FILE(MAKE_DIRECTORY ${VERSION_FILE_OUT_DIR})
# insert values into file creating PyCAVersion.py
CONFIGURE_FILE(
  ${VERSION_FILE_IN}
  ${VERSION_FILE_OUT}
  )
