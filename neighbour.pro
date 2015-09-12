TEMPLATE = app
QMAKE_CXXFLAGS += -fopenmp
CONFIG += console
CONFIG -= qt

LIBS += -lgomp -lpthread \
        -L/opt/H5Part/1.6.6/lib -lH5Part \
        -L/opt/hdf5/1.8.12/lib -lhdf5 \

INCLUDEPATH += /usr/local/include \
               /opt/H5Part/1.6.6/include \
               /opt/hdf5/1.8.12/include \
               /opt/boost/1.55.0/include

SOURCES += main.cpp \
    search.cpp \
    geninput.cpp \
    neighbour_search.cpp \
    direct.cpp \
    link_list_algorithm.cpp

HEADERS += \
    search.h \
    geninput.h \
    global.h \
    timer.hpp \
    timsort.hpp \
    neighbour_search.h \
    direct.h \
    link_list_algorithm.h

TARGET = neighbour
