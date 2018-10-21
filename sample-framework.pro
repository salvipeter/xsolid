# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG *= c++14 qt opengl
QT += gui widgets opengl xml

TRANSFINITE = /home/salvi/project/transfinite
INCLUDEPATH += $${TRANSFINITE}/src/geom $${TRANSFINITE}/src/transfinite
LIBS += -L$${TRANSFINITE}/release/geom -L$${TRANSFINITE}/release/transfinite

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp

LIBS *= -lQGLViewer -L/usr/lib/OpenMesh -lOpenMeshCore -ltransfinite -lgeom -lgle -lGL -lGLU

RESOURCES = sample-framework.qrc
