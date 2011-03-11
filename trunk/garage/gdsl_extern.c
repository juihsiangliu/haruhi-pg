#include "gdsl_extern.h"



static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER ;



static int gdsl_queue_copy_element(const gdsl_element_t E, gdsl_location_t LOCATION, void *USER_DATA)
{
	gdsl_queue_t dest_queue = (gdsl_queue_t) USER_DATA;
	gdsl_element_t data_in_queue = E; 	
	gdsl_queue_insert(dest_queue,data_in_queue);
}




gdsl_queue_t gdsl_queue_copy(const gdsl_queue_t const src)
{
	gdsl_queue_t dest = malloc(sizeof(struct _gdsl_queue));
//	gdsl_queue_t dest = gdsl_queue_alloc("copy of src",gdsl_alloc_func_t,gdsl_free_func_t);
//	gdsl_queue_map_forward(src,gdsl_queue_copy_element,dest);
	return dest;
}
