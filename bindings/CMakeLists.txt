cmake_minimum_required(VERSION 3.14)

project(
  mxtopotoolbox
  VERSION 0.0.1
  LANGUAGES C
)

include(FindMatlab)

include(FetchContent)
FetchContent_Declare(
  topotoolbox
  GIT_REPOSITORY https://github.com/TopoToolbox/libtopotoolbox.git
  GIT_TAG 69b4cee41208feb843333c732fd055ee5937cae9
)
FetchContent_MakeAvailable(topotoolbox)

matlab_add_mex(
  NAME tt_has_topotoolbox
  MODULE
  SRC tt_has_topotoolbox.c
  R2018a
  LINK_TO topotoolbox
)

matlab_add_mex(
  NAME tt_fillsinks
  MODULE
  SRC tt_fillsinks.c
  R2018a
  LINK_TO topotoolbox
)

matlab_add_mex(
  NAME tt_identifyflats
  MODULE
  SRC tt_identifyflats.c
  R2018a
  LINK_TO topotoolbox
)

matlab_add_mex(
  NAME tt_gwdt_computecosts
  MODULE
  SRC tt_gwdt_computecosts.c
  R2018a
  LINK_TO topotoolbox
)

install(
  TARGETS tt_has_topotoolbox tt_fillsinks tt_identifyflats tt_gwdt_computecosts
  DESTINATION "."
  COMPONENT Runtime
)
