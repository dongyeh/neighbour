TEMPLATE = app
QMAKE_CXXFLAGS += -fopenmp
CONFIG += console debug
CONFIG -= qt

LIBS += -lgomp -lpthread\

INCLUDEPATH += /export/home/dwang/opt/boost/1.51.0/include

SOURCES += main.cpp \
    search.cpp \
    geninput.cpp \
    neighbour_search.cpp \
    direct.cpp

HEADERS += \
    search.h \
    geninput.h \
    global.h \
    timer.hpp \
    timsort.hpp \
    neighbour_search.h \
    direct.h

TARGET = neighbour
