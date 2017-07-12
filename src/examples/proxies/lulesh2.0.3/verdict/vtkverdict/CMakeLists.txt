

configure_file(verdict.h.in ${CMAKE_CURRENT_BINARY_DIR}/verdict.h )
set( verdict_SRCS
  V_EdgeMetric.cpp
  V_GaussIntegration.cpp
  V_HexMetric.cpp
  V_KnifeMetric.cpp
  V_PyramidMetric.cpp
  V_QuadMetric.cpp
  V_TetMetric.cpp
  V_TriMetric.cpp
  VerdictVector.cpp
  V_WedgeMetric.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/verdict.h
  )


add_library(verdict ${verdict_SRCS} )
set_target_properties( ${verdict_LIBRARY} PROPERTIES DEFINE_SYMBOL verdict_EXPORTS)
