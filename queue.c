typedef struct queue_t
{
   struct queue_t *prev ;  // aponta para o elemento anterior na fila
   struct queue_t *next ;  // aponta para o elemento seguinte na fila
} queue_t ;


void queue_append (queue_t **queue, queue_t *elem){
	if(queue==NULL)
	{
		printf("Queue does not exist\n");
	}
	else if(elem==NULL)
	{
		printf("Element does not exist\n");
	}
	else if(elem->prev != NULL || elem->next != NULL)
	{
		printf("Element already in a queue\n");
	}
	else{

		if(*queue==NULL)//se a fila estiver vazia
		{ 
			//printf("Inserting first element\n");
			*queue = elem;
			elem->next = elem;
			elem->prev = elem;
			//printf("Inserted first element\n");
		}

		else{
			//printf("Inserting another element\n");
			queue_t* aux = *queue;//(*queue) corresponde ao primeiro elemento
			elem->next = aux;	   //elem->next aponta para o primeiro
			elem->prev = aux->prev;//elem->prev aponta para o último
			aux->prev->next = elem;//ultimo.next aponta para elem
			aux->prev = elem;	   //primeiro.prev aponta para o último
			//printf("Inserted another element\n");
		}
	}
	//printf("Entered queue_append\n");	
}

queue_t *queue_remove (queue_t **queue, queue_t *elem){
	if(queue == NULL)
	{
		printf("Queue does not exist\n");
		return NULL;
	}
	else if(elem == NULL)
	{
		printf("Element does not exist\n");
		return NULL;
	}
	else if(*queue == NULL)
	{
		printf("Queue is empty\n");
		return NULL;
	}
	else
	{
		
		queue_t* aux; 
		if(elem->next == elem && elem->prev == elem && *queue == elem)//se houver apenas um elemento
		{	
			*queue = NULL;
			elem->next = NULL;
			elem->prev = NULL;
			return elem;
		}
		//printf("elem      :%p\nqueue     :%p\n*queue    :%p\n&(*queue) :%p\n&queue    :%p\n&(**queue):%p\n",elem,queue,(*queue),&(*queue),&queue,&(**queue));
		if(elem == *queue)//se o elemento a remover for o primeiro da fila
		{	
			elem->next->prev = (*queue)->prev;//segundoElem.prev = ultimoElem
			(*queue)->prev->next = elem->next;
			*queue = elem->next;
			elem->next = NULL;
			elem->prev = NULL;
			return elem;
		}

		aux = *queue;
		while(aux->next!=elem)
		{	
			aux = aux->next;
			if(aux == *queue)//se tiver percorrido um ciclo
			{
				printf("Element not in queue\n");
				return NULL;
			}

		}
		aux->next = elem->next;
		elem->next->prev = aux;
		elem->next = NULL;
		elem->prev = NULL;
		return elem;
	}
}

void queue_print (char *name, queue_t *queue, void print_elem (void*) ){
	queue_t* aux = queue;
	if(queue==NULL){
		//printf("Queue is empty\n");
		printf("%s: []\n",name);
		return;
	}
	printf("%s:",name);
	printf(" [");
	print_elem(aux);
	while(aux->next != queue){
		aux = aux->next;
		printf(" ");
		print_elem(aux);	
	}
	printf("]\n");
}

int queue_size (queue_t *queue){
	queue_t* aux = queue;
	if(queue==NULL){
		printf("Queue is empty\n");
		return 0;
	}
	int count = 1;
	while(aux->next != queue){
		aux = aux->next;
		count++;
	}
	return count;
}

// int main(){
// 	queue_t myQueue;
// 	myQueue.name = 'a';
// 	myQueue.prev = NULL;
// 	myQueue.next = NULL;
// 	myQueue.count = 0;
// 	queue_t* ptrToMyQueue = &myQueue;
// 	queue_t anotherElem;
// 	//printf("%p\n",anotherElem.prev);
// 	queue_t anotherElem2;
// 	queue_t anotherElem3;
// 	anotherElem.prev = NULL; 
// 	anotherElem2.prev = NULL; 
// 	anotherElem3.prev = NULL;

// 	//printf("%p\n",anotherElem.prev);

// 	queue_append(&ptrToMyQueue,&anotherElem);
// 	queue_append(&ptrToMyQueue,&anotherElem2);
// 	queue_append(&ptrToMyQueue,&anotherElem3);

// 	queue_t *currElem = &myQueue;
// 	while(currElem->next!=NULL){
// 		//printf("(%c) count: %d\n",currElem->name, currElem->count);
// 		if(currElem->count != 0){
// 			printf("prev elem name: %c, curr elem name:%c\n",currElem->prev->name, currElem->name);
// 		}
		
// 		currElem = currElem->next;
// 	}
// 	//printf("(%c) count: %d\n",currElem->name, currElem->count);
// 	printf("prev elem name: %c, curr elem name:%c\n",currElem->prev->name, currElem->name);

// 	printf("Total elements: %d\n", queue_size(&myQueue));

// 	queue_t* removedElem = queue_remove(&ptrToMyQueue, &anotherElem2);
// 	printf("removed elem name: %c\n",removedElem->name);

// 	printf("Total elements: %d\n", queue_size(&myQueue));

// 	currElem = &myQueue;
// 	while(currElem->next!=NULL){
// 		//printf("(%c) count: %d\n",currElem->name, currElem->count);
// 		if(currElem->count != 0){
// 			printf("prev elem name: %c, curr elem name:%c\n",currElem->prev->name, currElem->name);
// 		}
		
// 		currElem = currElem->next;
// 	}
// 	//printf("(%c) count: %d\n",currElem->name, currElem->count);
// 	printf("prev elem name: %c, curr elem name:%c\n",currElem->prev->name, currElem->name); 

// }