#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<float.h>
#include"queue.c"
#include<string.h>

//Desenvolvido por
//Arthur Miguel Seehagen de Oliveira RA:1793853
//Felipe Manikowski	                 RA:

//Estrutura Node:
//armazena o nó de um grafo
//utilizada para:
//	- implementar o grafo por listas de adjacências
//	- armazenar os pontos de entrada em um vetor do tipo Node
//	- armazenar os pontos de entrada em uma heap de mínimo
//onde é utilizada:
//	- função prim()
typedef struct Node{
	struct Node* prev;
	struct Node* next;
	int x;
	int y;
	int id;
	int parentId;
	float weight;
	int visited;
}Node;

//Estrutura Graph:
//grafo representado por matriz de adjacências
typedef struct Graph{
	int size;
	float* adjMatrix;
	Node** adjList;
}Graph;

Node* createNode(int x,int y,int id){
	Node* newNode = (Node*) malloc(sizeof(Node));
	newNode->x = x;
	newNode->y = y;
	newNode->next = NULL;
	newNode->prev = NULL;
	newNode->id = id;
	newNode->visited = 0;
	newNode->weight = -1;
	newNode->parentId = -1;
	return newNode;
}

float calcDist(int x1, int y1, int x2, int y2)
{
	return(sqrt(pow(x1-x2,2) + pow(y1-y2,2)));
}

float calcDistNodes(Node n1, Node n2){
	return(sqrt((n1.x-n2.x)*(n1.x-n2.x) + (n1.y-n2.y)*(n1.y-n2.y)));
}

// void print_elem (void *ptr)
// {
//    Node* elem = ptr;

//    if (!elem)
//       return ;

//    elem->prev ? printf ("%d,%d", elem->prev->x, elem->prev->y) : printf ("*") ;
//    printf ("<%d> <%d>", elem->x, elem->y);
//    elem->next ? printf ("%d,%d", elem->next->x, elem->prev->y) : printf ("*") ;
// }

void print_elem(void *ptr)
{
   Node* elem = ptr;

   if (!elem)
      return ;
   printf ("(%d,%d)", elem->x, elem->y);
}


void printAdjList(Graph g, int n){
	int i;
	for(i=0;i<n;i++)
	{
		printf("Lista[%d]",i);
		queue_print (" ", (queue_t*) g.adjList[i], print_elem);
	}
}

void printArray(Node* array,int size){
	int i;
	for(i=0;i<size;i++){
		printf("%f, ",array[i].weight);
	}
	printf("\n");
}

//mapa de índices da heap
int* heapIndexMap;
//função que mantém a propriedade de uma heap de mínimo
void minHeapify(Node* A,int i, int size){
	int l = 2*i+1;
	int r = 2*i+2;
	int smallest;
	if(l<=size-1 && A[l].weight<A[i].weight){
		smallest = l;
	}
	else{
		smallest = i;
	}
	if(r<=size-1 && A[r].weight<A[smallest].weight){
		smallest = r;
	}
	if(smallest != i){
		heapIndexMap[A[i].id] = smallest;
		heapIndexMap[A[smallest].id] = i;
		Node aux = A[i];
		A[i] = A[smallest];
		A[smallest] = aux;
		minHeapify(A,smallest,size);
	}
}

//função que cria uma heap de mínimo
void buildMinHeap(Node* array,int size){
	int i;
	for(i=floor(size/2)-1;i>=0;i--){
		minHeapify(array,i,size);
	}
}

//função que extrai o menor elemento de uma heap de mínimo
Node heapExtractMin(Node* A, int heapsize){
	// int heapsize = sizeof(A)/sizeof(float);
	// printf("heapsize: %d\n",heapsize);
	if(heapsize<1){
		printf("Heap Underflow\n");
		exit(0);
	}
	else{
		Node min = A[0];
		heapsize--;
		heapIndexMap[A[heapsize].id] = 0;
		memcpy(&A[0],&A[heapsize],sizeof(Node));
		//A[0] = A[heapsize];
		A = (Node*)realloc(A,(heapsize)*sizeof(Node));
		minHeapify(A,0,heapsize);
		return min;
	}
}

//função que decrementa o valor de um elemento na heap de mínimo
void heapDecreaseKey(Node* A,int i,float key){
	if(key > A[i].weight){
		printf("New key is greater than current key\n");
		return;
	}
	A[i].weight = key;
	Node aux;
	while(i>0 && A[(int)((i-1)/2)].weight>A[i].weight){
		int parent = floor((i-1)/2);
		heapIndexMap[A[parent].id] = i;
		heapIndexMap[A[i].id] = parent;
		aux = A[i];
		A[i] = A[parent];
		A[parent] = aux;
		i = parent;
	}
}

//função que insere um elemento na heap de mínimo
void minHeapInsert(Node* A,int heapsize,float key){
	A = (Node*)realloc(A,(heapsize+1)*sizeof(Node));
	A[heapsize].weight = FLT_MAX;
	heapDecreaseKey(A,heapsize,key);
}



//função que inicializa o grafo por matriz de ajacências
//recebe um vetor do tipo Node como parâmetro
Graph initializeGraphMatrix(Node* vertex, int size){
	int i,j,x,y;
	Graph graph;
	graph.adjMatrix = (float*) malloc(size*size*sizeof(float));

	//preenche a matrix triangular superior que representa o grafo
	//um dado elemento Aij da matriz corresponde à distância euclidiana entre os nós i e j
	for(i=1;i<size;i++)
	{	
		x = vertex[i-1].x;
		y = vertex[i-1].y;
		for(j=i;j<size;j++)
		{	
			graph.adjMatrix[i*size+j] = calcDist(x,y,vertex[j].x,vertex[j].y);
		}
	}
	return graph;
}
//função que inicializa o grafo por listas de adjacências
//recebe um vetor do tipo Node como parâmetro
Graph* initializeGraphAdjList(Node* vertex, int size){
	int i,j;
	Graph* graph = (Graph*)malloc(sizeof(Graph));
	graph->adjList = (Node**)malloc(size*sizeof(Node*));
	for(i=0;i<size;i++)
	{
		graph->adjList[i]=NULL;
		for(j=0;j<size;j++)
		{
			if(j==i)
			{
				continue;
			}

			Node* currentNode = createNode(vertex[j].x,vertex[j].y,vertex[j].id);
			// Node* elem = &currentNode;
			// printf("(%d, %d)\n",elem->x,elem->y);
			//printf(" i: %d, j: %d\n",i,j);

			queue_append((queue_t**)&graph->adjList[i], (queue_t*)currentNode);
			//printAdjList(graph,size);
		}
	}
	return graph;
}

//função que insere um nó em um grafo
void insertNode(Graph* g, Node* n){
	g->size+=1;
	printf("g size:%d\n",g->size);
	n->id = g->size;
	g->adjList = (Node**)realloc(g->adjList,(g->size)*sizeof(Node*));
	if(g->adjList==NULL){
		printf("error\n");
	}
	printf("after realloc\n");
	// Node* newNode = (Node*)malloc(sizeof(Node));
	// printf("after malloc\n");
	// *newNode = n;
	queue_append((queue_t**)&g->adjList[g->size-1], (queue_t*)n);
}

//função que adiciona um nó a uma lista de adjacências
void insertAdjList(Graph* g, Node* n, int i){
	queue_append((queue_t**)&g->adjList[i], (queue_t*)n);
}


//função que adiciona uma aresta ao grafo
void insertEdge(Graph* g,Node* n1,Node* n2){
	Node* newNode1 = malloc(sizeof(Node));
	Node* newNode2 = malloc(sizeof(Node));
	newNode1 = n1;
	newNode2 = n2;
	queue_append((queue_t**)&g->adjList[n1->id], (queue_t*)newNode2);
	queue_append((queue_t**)&g->adjList[n2->id], (queue_t*)newNode1);
}

// void initializeGraphAdjList2(Graph* graph, Node* vertex, int size){
// 	int i,j;
// 	printf("b\n");
// 	graph->adjList = malloc(size*sizeof(Node*));
// 	for(i=0;i<size;i++)
// 	{
// 		for(j=0;j<size;j++)
// 		{
// 			if(j==i)
// 			{
// 				continue;
// 			}
// 			Node currentNode;
// 			currentNode.x = vertex[j].x;
// 			currentNode.y = vertex[j].y;
// 			currentNode.prev = NULL;
// 			currentNode.next = NULL;
// 			Node* elem = &currentNode;
// 			printf("(%d, %d)\n",elem->x,elem->y);
// 			printf(" i: %d, j: %d\n",i,j);
// 			queue_append((queue_t**)&graph->adjList[i], (queue_t*)&currentNode);
// 			printf("(%d, %d)\n",graph->adjList[i]->x,graph->adjList[i]->y);
// 			printf("(%d, %d)\n",graph->adjList[i]->next->x,graph->adjList[i]->next->y);
// 			printAdjList(*graph,size);
// 		}
// 	}
// }

void printGraphMatrix(Graph g, int n){
	int i,j,k;
	for(i=1;i<n;i++)
	{
		k = 1;
		while(i-k>0)
		{
			printf("       ");
			k++;
		}
		for(j=i;j<n;j++)
		{
			printf("%.3f, ",g.adjMatrix[i*n+j]);
		}
		printf("\n");
	}
}

//fila que armazena o ciclo
Node* cycle = NULL;
//variável que armazena o custo total do ciclo
float cycleCost = 0;

void DFSVisit(Graph graph, Node* points, int u){
	points[u].visited = 1;
	//printf("u:%d\n",u);
	Node* aux = graph.adjList[u];
	if(aux!=NULL){
		if(aux->visited == 0){
			//printf("2 visiting: (%d,%d)\n",aux->x,aux->y);
			points[aux->id].parentId = u;
			points[aux->id].visited = 1;
			Node* v = createNode(aux->x,aux->y,aux->id);
			queue_append((queue_t**)&cycle, (queue_t*)v);
			cycleCost+=calcDistNodes(*(v->prev),*v);
			//printf("2 appended: (%d,%d)\n",v->x,v->y);
			//queue_print ("Cycle: ", (queue_t*)cycle, print_elem);
			DFSVisit(graph,points,v->id);
		}
		while(aux->next!= graph.adjList[u]){
			aux = aux->next;
			if(aux->visited == 0){
				//printf("3 visiting: (%d,%d)\n",aux->x,aux->y);
				points[aux->id].parentId = u;
				points[aux->id].visited = 1;
				Node* v = createNode(aux->x,aux->y,aux->id);
				queue_append((queue_t**)&cycle, (queue_t*)v);
				cycleCost+=calcDistNodes(*(v->prev),*v);
				//printf("3 appended: (%d,%d)\n",v->x,v->y);
				//queue_print ("Cycle: ", (queue_t*)cycle, print_elem);
				DFSVisit(graph,points,v->id);
			}
		}
	}

}

//função que realiza a busca em profundidade e preenche uma lista de nós
void DFS(Graph graph, Node* points, int size){
	int i;
	for(i=0;i<size;i++){
		points[i].visited =  0;
	}
	for(i=0;i<size;i++){
		if(points[i].visited == 0){
			//printf("1 visiting: (%d,%d)\n",points[i].x,points[i].y);
			Node* v = createNode(points[i].x,points[i].y,points[i].id);
			queue_append((queue_t**)&cycle, (queue_t*)v);
			//printf("1 appended: (%d,%d)\n",v->x,v->y);
			//queue_print ("Cycle: ", (queue_t*)cycle, print_elem);
			DFSVisit(graph,points,i);
		}
	}
	Node* v = createNode(points[0].x,points[0].y,points[0].id);
	queue_append((queue_t**)&cycle, (queue_t*)v);
	cycleCost+=calcDistNodes(*(v->prev),*v);
}


//função que gera a árvore geradora mínima
Graph prim(Graph graph, Node* points, int size){
	int i,k;
	Node u;
	Graph tree;
	tree.size = size;
	tree.adjList = (Node**)malloc(size*sizeof(Node*));
	Node *heap = (Node*)malloc(size*sizeof(Node));
	for(i=1;i<size;i++){
		heapIndexMap[i] = points[i].id;
		heap[i] = points[i];
		heap[i].weight = FLT_MAX;
		points[i].weight = FLT_MAX;
		points[i].visited = 0;
		tree.adjList[i]=NULL;
	}
	tree.adjList[0]=NULL;
	heapIndexMap[0] = points[0].id;
	heap[0] = points[0];
	heap[0].weight = 0;
	points[0].weight = 0;
	buildMinHeap(heap,size);
	for(i=size;i>0;i--)
	{
		//printf("i:%d, size:%d\n",i,size);
		//printf("Heap Before Extract\n");
		//printArray(heap,i);
		u = heapExtractMin(heap,i);
		points[u.id].visited = 1;
		//printf("Heap After Extract\n");
		//printArray(heap,i-1);
		//printf("Extracted:\n(%d,%d), id:%d\n",u.x,u.y,u.id);
		if(i<size)
		{
			Node* v = createNode(u.x,u.y,u.id);
			v->weight = u.weight;
			insertAdjList(&tree,v,u.parentId);
			//printf("\nCurrent Graph\n");
			//printAdjList(tree,size);
			//printf("(%d,%d)\n",tree.adjList[0]->x,tree.adjList[0]->y);
		}

		
		
		Node* aux = graph.adjList[u.id];
		//printf("aux: (%d,%d)\n",graph.adjList[u.id]->x,graph.adjList[u.id]->y);
		
		
		k=0;
		//percorre a lista de adjacência do nó "u"
		while(k<size-1)
		{
			//printf("current aux: ");
			//print_elem(aux);
			//printf("\n");
			//printf("aux id: %d, visited:%d\n",aux->id,points[aux->id].visited);
			float w = calcDistNodes(u,*aux);
			if(points[aux->id].visited==0 && w<points[aux->id].weight)
			{
				//printf("true\n");
				//printf("heap[%d].weight:%f, w:%f\n",j,heap[j].weight,w);
				heap[heapIndexMap[aux->id]].parentId = u.id;
				heapDecreaseKey(heap,heapIndexMap[aux->id],w);
				points[aux->id].weight = w;
				
			}
			aux = aux->next;
			k++;	
		}	
	}
	return tree;
}





int main(){
	clock_t start,end,total;
	FILE *input;
	input = fopen("input.txt", "r");
	int n;
	fscanf(input, "%d", &n);
	Node* points = (Node*) malloc(n*sizeof(Node));
	int i = 0;

	while (EOF != fscanf(input, "%d %d", &points[i].x, &points[i].y))
    {
    	points[i].next = NULL;
		points[i].prev = NULL;
		points[i].id = i;
        i++;
    }

    //free(points);
    //Graph graph = initializeGraphMatrix(points,n);
    heapIndexMap = (int*)malloc(n*sizeof(int));
    Graph* graph2 = initializeGraphAdjList(points,n);
    //printGraphMatrix(graph,n);
    //printAdjList(graph2,n);

    //float* heap = (float*)malloc(10*sizeof(float));
    // for(i=0;i<n;i++){
    // 	points[i].weight=n-i;
    // }
    // printArray(points,n);
    // buildMinHeap(points,n);
    // printArray(points,n);

    // i=n;
    // while(i>=0){
    // 	printf("current size: %d\n", i);
    // 	printArray(points,i);
    // 	Node min = heapExtractMin(points,i);
    // 	printf("min: %f\n",min.weight);
    // 	// sleep(5);
    // 	i--;
    // }

    // heapDecreaseKey(points,5,0);
    // printArray(points,n);
    // heapDecreaseKey(points,1,1);
    // printArray(points,n);
    // minHeapInsert(points,n,11.0);
    // printArray(points,n+1);

    // Node n1;
    // n1.x = 10;
    // n1.y = 10;
    // n1.prev = NULL;
    // n1.next = NULL;
    // Node n2;
    // n2.x = 11;
    // n2.y = 11;
    // n2.prev = NULL;
    // n2.next = NULL;

    // insertNode(&graph2,&n1);
    // insertNode(&graph2,&n2);
    // printAdjList(graph2,graph2.size);
    // insertEdge(&graph2,&n1,&n2);
    //insertAdjList(&graph2,&n1,5);
    //printAdjList(graph2,graph2.size);

    start = clock();
    Graph mst = prim(*graph2,points,n);
    end = clock();
    for(i=0;i<n;i++){
    	Node* aux = graph2->adjList[i];
    	while(aux!=NULL){
			Node* elem = (Node*)queue_remove ((queue_t**) &aux, (queue_t*)aux);
			//printf("removed: (%d,%d)\n",elem->x,elem->y);
			free(elem);
    	}
    }
    total=end-start;
    //free(graph2);
    //printAdjList(mst,n);
    start = clock();
    DFS(mst,points,n);
    end = clock();
    total+=end-start;
    printf("%.6f %.6f\n",(float)(total)/CLOCKS_PER_SEC,cycleCost);
    //queue_print ("Cycle: ", (queue_t*) cycle, print_elem);
    //printf("Cycle cost: %f\n",cycleCost);

    fclose(input);

    FILE *output = fopen("cycle.txt", "w");
    Node* aux = cycle;
    fprintf(output, "%d %d\n",aux->x,aux->y);
    while(aux->next!=cycle){
    	aux = aux->next;
    	fprintf(output, "%d %d\n",aux->x,aux->y);
    }
    fclose(output);

    output = fopen("tree.txt", "w");
    for(i=0;i<n;i++)
    {
    	Node* aux = mst.adjList[i];
    	if(aux!=NULL)
    	{
    		fprintf(output, "%d %d\n",points[i].x,points[i].y);
    		fprintf(output, "%d %d\n\n",aux->x,aux->y);
    		while(aux->next!=mst.adjList[i]){
    			aux = aux->next;
    			fprintf(output, "%d %d\n",points[i].x,points[i].y);
    			fprintf(output, "%d %d\n\n",aux->x,aux->y);
    		}
    	}	
    }
    fclose(output);
	


}