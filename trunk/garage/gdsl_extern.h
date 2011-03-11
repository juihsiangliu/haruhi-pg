#ifndef GDSL_EXTERN_H
#define GDSK_EXTERN_H

#include <stdlib.h>
#include <pthread.h>
#include <gdsl_queue.h>
#include <gdsl_types.h>

// this is a thread safe function... i use a mutex to lock this function ...
gdsl_queue_t gdsl_queue_copy(const gdsl_queue_t const src);

#endif
