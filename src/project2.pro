TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    functions.cpp

HEADERS += \
    functions.h
