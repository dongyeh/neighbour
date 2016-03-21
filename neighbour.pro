TEMPLATE = app
QMAKE_CXXFLAGS += -fopenmp
CONFIG += console
CONFIG -= qt

LIBS += -lgomp -lpthread \
        -L/opt/hdf5/1.8.16/lib -lhdf5 \
        -L/opt/H5hut/lib -lH5hut \

INCLUDEPATH += ./include \
               /usr/include \
               /opt/hdf5/1.8.16/include \
               /opt/H5hut/include \

SOURCES += main.cpp \
    src/search.cpp \
    src/geninput.cpp \
    src/neighbour_search.cpp \
    src/direct.cpp \
    src/link_list_algorithm.cpp

HEADERS += \
    include/search.h \
    include/geninput.h \
    include/global.h \
    include/timer.hpp \
    include/timsort.hpp \
    include/neighbour_search.h \
    include/direct.h \
    include/link_list_algorithm.h

TARGET = bin/neighbour
