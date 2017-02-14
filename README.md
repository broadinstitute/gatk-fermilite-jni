# gatk-fermilite-jni
JNI code for fermi lite.

This project builds dynamic libraries with a JNI API.
It allows Java code to call Heng Li's fermi lite assembler.

The makefile in src/main/c will build an appropriate library for Mac OSX or x86_64 Linux.
Pre-compiled dynamic libraries for these OS's exist in src/main/resources.
