// Compatibility shims for static linking against older glibc targets.
// These symbols are referenced by bundled static dependencies during
// fully-static builds with the conda cross toolchain.

#if defined(__linux__)

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>
#include <time.h>

extern "C" ssize_t getrandom(void* buf, size_t buflen, unsigned int flags)
{
#ifdef SYS_getrandom
   return syscall(SYS_getrandom, buf, buflen, flags);
#else
   (void)buf;
   (void)buflen;
   (void)flags;
   errno = ENOSYS;
   return -1;
#endif
}

extern "C" int clock_gettime(clockid_t clk_id, struct timespec* tp)
{
#ifdef SYS_clock_gettime
   return (int)syscall(SYS_clock_gettime, clk_id, tp);
#else
   (void)clk_id;
   (void)tp;
   errno = ENOSYS;
   return -1;
#endif
}

#endif
