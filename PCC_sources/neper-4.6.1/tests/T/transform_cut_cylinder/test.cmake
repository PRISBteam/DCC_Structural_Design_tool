# This file is part of the Neper software package.
# Copyright(C) 2003-2022, Romain Quey.
# See the COPYING file in the top-level directory.

set (test_command -T -n 20 -transform "cut(cylinder(1.5,1.4,1.3,0,1,2,1.2))" -o test)
include(../../test.cmake)
