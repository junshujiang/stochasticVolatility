cmake_minimum_required(VERSION 3.16)
project(cloglike LANGUAGES C)

# ===== 手动配置 =====
# 生成的动态库名（不带前缀/后缀）
set(LIB_NAME gaussiantest)

# 头文件目录列表：把需要的 include 目录都写在这里
# 例：INLA 的 include、gmrflib 目录等
set(INCLUDE_DIRS
  /home/jianj0c/anaconda3/envs/my_r_env/lib/R/library/INLA/include         # 含 cgeneric.h
 /home/jianj0c/project/STP/notebook_V/tmp/r-inla/gmrflib         # 含 interpol.h（如 pc-priors.c 需要）
)

# 源文件列表：把需要编译进 .so 的所有 .c 全部列出来
# 例：你的 cloglike 源、pc-priors.c、gmrflib/interpol.c
set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/gaussiantest.c
 /home/jianj0c/project/STP/notebook_V/tmp/r-inla/inlaprog/src/pc-priors.c
 /home/jianj0c/project/STP/notebook_V/tmp/r-inla/gmrflib/interpol.c
  # /more/sources.c ...
)

# 可选：额外链接库（按需填写，比如需要 GSL 时）
# set(LINK_LIBS gsl gslcblas)

add_library(${LIB_NAME} SHARED ${SOURCE_FILES})

target_include_directories(${LIB_NAME} PRIVATE ${INCLUDE_DIRS})
target_compile_options(${LIB_NAME} PRIVATE -Wall -O2 -fPIC)
target_link_libraries(${LIB_NAME} m)

# 产物名不加 'lib' 前缀（得到 <basename>.so）
set_target_properties(${LIB_NAME} PROPERTIES PREFIX "")

# 打印信息（便于核对）
message(STATUS "Library name : ${LIB_NAME}")
message(STATUS "Include dirs : ${INCLUDE_DIRS}")
message(STATUS "Source files : ${SOURCE_FILES}")