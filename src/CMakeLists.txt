cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH "../cmake-modules/")
ENABLE_TESTING()

FIND_PACKAGE(RLibs)

include_directories (.)

link_libraries(${R_LIBRARY})

add_executable (sMQM 
MQMdata.cpp        MQMmapQTL.h     MQMreDefine.h      MQMsupport.cpp
MQMdata.h          MQMmixture.cpp  MQMRegression.cpp  MQMsupport.h
MQMinterfaces.cpp  MQMmixture.h    MQMRegression.h    MQMutil.cpp
MQMinterfaces.h    MQMprob.cpp     MQMscan.cpp        MQMutil.h
MQMmapQTL.cpp      MQMprob.h       MQMscan.h

)


