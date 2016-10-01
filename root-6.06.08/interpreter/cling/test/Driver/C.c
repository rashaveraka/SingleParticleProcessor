//------------------------------------------------------------------------------
// CLING - the C++ LLVM-based InterpreterG :)
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

// RUN: cat %s | %cling -x c | FileCheck %s

// Validate cling C mode.

int printf(const char*,...);
printf("CHECK 123\n"); // CHECK: CHECK 123

// fix value printing!
