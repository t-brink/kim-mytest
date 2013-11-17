#  Copyright (c) 2013 Tobias Brink
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# load all basic KIM make configuration
include ../Makefile.KIM_Config

# set test specific details
TEST_NAME := mytest

LOCALOBJ = core.o elastic.o dsl.o

LOCALCLEAN =

# APPEND to compiler option flag lists
#FFLAGS   +=
#CCFLAGS  +=
CXXFLAGS += -std=c++11 -Inlopt/include
LDFLAGS  += -Wl,-rpath,nlopt/lib -Lnlopt/lib
LDDYNAMICLIB += -lnlopt -lm

# load remaining KIM make configuration
include $(KIM_DIR)/MAKE_SYSTEM/Makefile.Test
