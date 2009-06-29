include_directories (.)

# add_definitions(`sdl-config --cflags`)
# link_libraries(SGE `sdl-config --libs` log4cxx dl)

add_executable (sMQM 
MQMdata.cpp        MQMmapQTL.h     MQMreDefine.h      MQMsupport.cpp
MQMdata.h          MQMmixture.cpp  MQMRegression.cpp  MQMsupport.h
MQMinterfaces.cpp  MQMmixture.h    MQMRegression.h    MQMutil.cpp
MQMinterfaces.h    MQMprob.cpp     MQMscan.cpp        MQMutil.h
MQMmapQTL.cpp      MQMprob.h       MQMScan.h

)

# add_custom_command(TARGET reversi PRE_BUILD COMMAND ./bn.sh)
