TEMPLATE = app
QMAKE_CXXFLAGS += -fopenmp
CONFIG += console
CONFIG -= qt

LIBS += -lgomp -lpthread -lhdf5 \
        -L/opt/H5Part/1.6.6/lib -lH5Part \

INCLUDEPATH += ./include \
               /usr/local/include \
               /opt/H5Part/1.6.6/include \

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
