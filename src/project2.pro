TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    functions.cpp \
    unity.cpp \

HEADERS += \
    functions.h \
    unity.h
