TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp
LIBS += -L/usr/local/lib/armadillo-4.400.1 -larmadillo
INCLUDEPATH += /usr/include
DEPENDPATH += /usr/include