# 声明CMake的最低版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名和版本
project(MyProject VERSION 1.0)

# 指定C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# # 添加包含目录（假设头文件在 "include" 文件夹内）
# include_directories(include)

# 添加可执行文件
set(SOURCE_FILES main.cpp basefunc.cpp functions.cpp in_output.cpp PKUMC_init.cpp statistic.cpp stoppingPower.cpp)
add_executable(MyProject ${SOURCE_FILES})

# # 链接库（假设库文件在 "lib" 文件夹内）
# link_directories(lib)
# target_link_libraries(MyProject someLibrary)

# 复制到构建目录（假设它与CMakeLists.txt在同一目录下）
set(FILES_TO_COPY
    parameters.txt
    h.88
    scoef.88
    # ... 其他文件
)
foreach(FILE ${FILES_TO_COPY})
    configure_file(${FILE} ${FILE} COPYONLY)
endforeach()