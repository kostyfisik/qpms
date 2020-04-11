#ifndef _QPMS_OSHACKS_H
#define _QPMS_OSHACKS_H
#include <unistd.h>

#ifdef _SC_NPROCESSORS_ONLN
static inline long get_ncpus(void) {
	return sysconf(_SC_NPROCESSORS_ONLN);
}
#elif (0)
#include <sys/types.h>
#include <sys/sysctl.h>
static inline long get_ncpus(void) {
	int32_t ncpu;
	size_t len = sizeof(ncpu);
	sysctlbyname("hw.physicalcpu", &ncpu, &len, NULL, 0);
	return ncpu;
}
#else
static inline long get_ncpus(void) { return -1; }
#endif
	

#endif // _QPMS_OSHACKS_H
