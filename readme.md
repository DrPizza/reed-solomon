A fairly straightforward port of Backblaze's [Reed Solomon library](https://www.backblaze.com/blog/reed-solomon/) available [here](https://github.com/Backblaze/JavaReedSolomon), ported to C++.

Includes SSE3/SSSE3 optimizations per [Screaming Fast galois Field Arithmetic](http://www.snia.org/sites/default/files2/SDC2013/presentations/NewThinking/EthanMiller_Screaming_Fast_Galois_Field%20Arithmetic_SIMD%20Instructions.pdf), and simple parallel encoding/verification using the Intel's TBB.

Only tested in Visual Studio 2015. Shouldn't be too hard to port it to other places.

The code is almost all in headers, as I find maintaining separate header/implementation pairs tedious beyond belief. Maybe one day I'll split it.

MIT licensed, as per the Backblaze original. Note that TBB is GPL licensed with the linking exemption. As such, it's OK to include and link it into projects with other licenses (both open and closed), with no source distribution requirement. But modified versions of TBB must be distributed with source.