cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH "../cmake-modules/")
ENABLE_TESTING()

FIND_PACKAGE(RLibs)

include_directories (.)
add_definitions(-DSTANDALONE)

link_libraries(${R_LIBRARY})

add_executable (sMQM 
MQMdata.cpp        MQMmapQTL.h     standalone.h       MQMsupport.cpp
MQMdata.h          MQMmixture.cpp  MQMregression.cpp  MQMsupport.h
MQMinterfaces.cpp  MQMmixture.h    MQMregression.h    
MQMinterfaces.h    MQMprob.cpp     MQMscan.cpp        
MQMmapQTL.cpp      MQMprob.h       MQMscan.h
util.h util.c      MQMmain.cpp
)


