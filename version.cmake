# based on https://www.mattkeeter.com/blog/2018-01-06-versioning/

execute_process(COMMAND git describe --tags --dirty --always
                OUTPUT_VARIABLE GIT_REV
                ERROR_QUIET)

# Check whether we got any revision (which isn't
# always the case, e.g. when someone downloaded a zip
# file from Github instead of a checkout)
if ("${GIT_REV}" STREQUAL "")
    set(GIT_REV "N/A")
    set(GIT_BRANCH "N/A")
else()
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        OUTPUT_VARIABLE GIT_BRANCH)

    string(STRIP "${GIT_REV}" GIT_REV)
    string(STRIP "${GIT_BRANCH}" GIT_BRANCH)
endif()

set(VERSION "const char* GIT_REV=\"${GIT_REV}\";
const char* GIT_BRANCH=\"${GIT_BRANCH}\";")

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/qpms_version.c)
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/qpms_version.c VERSION_)
else()
    set(VERSION_ "")
endif()

if (NOT "${VERSION}" STREQUAL "${VERSION_}")
    file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/qpms_version.c "${VERSION}")
endif()

