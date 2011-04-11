#include "vector.h"

static int get_element_size(struct Vector *this)
{
	return this->_element_size;
}




static int get_size(struct Vector *this)
{
	return this->_size;
}




static int get_max_size(struct Vector *this)
{
	return this->_max_size;
}



static void enlarge(struct Vector *this,const int new_max)
{
	void *buf = malloc(this->_element_size*new_max);
	memset(buf,0,this->_element_size*new_max);
	memcpy(buf,this->_data,this->_element_size*this->_max_size);
	free(this->_data);
	this->_data = buf;
	this->_max_size = new_max;
}




static void push_back(struct Vector *this,void *data)
{
	if(this->_size == this->_max_size) enlarge(this,this->_max_size*2);
	memcpy(this->_data+this->_element_size*this->_size, data, this->_element_size);
	this->_size++;
}




static void *get(struct Vector *this,const int idx)
{
	assert(idx<this->_size);
	return this->_data+this->_element_size*idx;
}


static void set(struct Vector *this,const int idx,void *data)
{
	if(idx >= this->_max_size) enlarge(this,2*idx);
	memcpy(this->_data+this->_element_size*idx,data,this->_element_size);
	
	if(idx >= this->_size) this->_size = idx+1;
}



static void erase(struct Vector *this,const int idx)
{
	assert(idx < this->_size);
	memset(this->_data+this->_element_size*idx,0,this->_element_size);
}




static void swap(struct Vector *this,const int idx_a,const int idx_b)
{
	assert(idx_a < this->_size);
	assert(idx_b < this->_size);

	void *tmp = malloc(this->_element_size);
	memcpy(tmp,this->get(this,idx_a),this->_element_size);  // tmp = a
	memcpy(this->_data+this->_element_size*idx_a,this->_data+this->_element_size*idx_b,this->_element_size); // a = b
	memcpy(this->_data+this->_element_size*idx_b,tmp,this->_element_size); // b = tmp
	free(tmp);
}


static void *front(struct Vector *this)
{
	return this->_data;
}



static void *back(struct Vector *this)
{
	return this->_data+this->_element_size*(this->_size-1);
}





struct Vector *create_size_Vector(const int element_size,const int vector_size)
{
	struct Vector *ptr = malloc(sizeof(struct Vector));
	
	// element
	ptr->_element_size = element_size;
	ptr->_size = 0;
	ptr->_max_size = vector_size;
	ptr->_data = malloc(ptr->_element_size*ptr->_max_size);
	memset(ptr->_data,0,ptr->_element_size*ptr->_max_size);
	// function
	ptr->get_element_size = get_element_size;
	ptr->get_size = get_size;
	ptr->get_max_size = get_max_size;
	ptr->push_back = push_back;
	ptr->get = get;
	ptr->set = set;
	ptr->erase = erase;
	ptr->swap = swap;
	ptr->front = front;
	ptr->back = back;

	return ptr;
	
}


struct Vector *create_Vector(const int element_size)
{
	return create_size_Vector(element_size,VECTOR_MAX_SIZE);
}


void free_Vector(struct Vector *this)
{
	free(this->_data);
	free(this);
}










