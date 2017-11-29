#-------------------------------------------------
#
# Project created by QtCreator 2016-11-05T18:41:36
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = swimmer_3d
CONFIG   += console
CONFIG   -= app_bundle
CONFIG += c++11


TEMPLATE = app


SOURCES += main.cpp \
    spring.cpp \
    mass.cpp \
    vector3d.cpp \
    environment.cpp \
    angularspring.cpp

HEADERS += \
    spring.h \
    mass.h \
    vector3d.h \
    environment.h \
    angularspring.h
