#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif


void checklist(double *lst,int n){//checks for error in list
	int i;
	for(i=0;i<n-1;i++){
		if(lst[i]>lst[i+1])
			printf("ERROR at %d\n",i);
	}
	printf("%d elements checked\n",i+1);
}

void printlist(double *lst,int n){//prints list
	for(int i=0;i<n;i++){
		printf("%f\n",lst[i]);
	}
}

void shellsort(double *lst,int lstsize){ //shellsort algorithm 
	int i,j,p,gap,n;
	double temp;
	#ifdef _OPENMP
	int threads = omp_get_max_threads();//sets amount of threads to use
	#endif

	//Ciura 2001(fastest)
	int gaps[8] = {701, 301, 132, 57, 23, 10, 4, 1}; //gap sizes
	
	int gapsize = sizeof(gaps)/sizeof(int);
	for(gap = 0;gap<gapsize;gap++){
		n = gaps[gap]; //gap size
		#pragma omp parallel for shared(lst,n,lstsize,threads) private(i,j,p,temp) num_threads(threads) schedule(dynamic,threads*2)
		for(i = 0; i<lstsize-n;i++){
			p = i;
			for(j = 0;j<=p;j+=1){
				if(lst[p]>lst[p+n]){ //replaces elements that are n apart if the earlier element is bigger
					temp = lst[p];
					lst[p] = lst[p+n];
					lst[p+n] = temp;
				}
				p-=n;
			}
		}
	}
	
	//checklist(lst,lstsize);
	

}

int main(){
	srand(time(NULL));//makes list elements random every time the program runs
	int size = 200000;

	double *list_to_sort;
	list_to_sort = (double *) malloc(size*sizeof(double));//malloc to be able to allocate more memory
	
	
	for(int i=0;i<size;i++)//creates list and its elements 
		list_to_sort[i] = rand()%100;

	shellsort(list_to_sort,size);
	free(list_to_sort);
	return 0;	
}

//Other implementations that were tested

/* 
//First iteration using gapsize n/2 for every loop

n = lstsize;
for(n/=2;n>0;n/=2){
	#pragma omp parallel for shared(lst,n,lstsize,threads) private(i,j,p,temp) num_threads(threads) schedule(dynamic,threads*2)
	for(i = 0; i<lstsize-n;i++){
		p = i;
		for(j = 0;j<=i;j+=n){
			if(lst[p]>lst[p+n]){
				temp = lst[p];
				lst[p] = lst[p+n];
				lst[p+n] = temp;
			}
			p-=n;
		}
	}
}
*/

/*
//Hibbard 1963(slower than the two others)

int size= 8,k=size,gaps[size],pos=0;
int value = pow(2,k)-1;
while(value>=1){
	gaps[pos] = value;
	k--;
	pos++;
	value = pow(2,k)-1;
}
*/

/*
//Tokuda 1992(close to Ciura)

int size= 8,k=size,gaps[size],pos=0;
double value = (1.0/5.0)*(9*pow((9.0/4.0),k-1)-4);
while(value>1){
	gaps[pos] = (int)value+1;
	printf("%d ",gaps[pos]);
	k--;
	pos++;
	value = (1.0/5.0)*(9*pow((9.0/4.0),k-1)-4);
}
gaps[pos] = 1;
*/

/*
//Ciura 2001(fastest)

int gaps[8] = {701, 301, 132, 57, 23, 10, 4, 1};
	
	int gapsize = sizeof(gaps)/sizeof(int);
	for(gap = 0;gap<gapsize;gap++){
		n = gaps[gap];
		#pragma omp parallel for shared(lst,n,lstsize,threads) private(i,j,p,temp) num_threads(threads) schedule(dynamic,threads*2)
		for(i = 0; i<lstsize-n;i++){
			p = i;
			for(j = 0;j<=p;j+=1){
				if(lst[p]>lst[p+n]){
					temp = lst[p];
					lst[p] = lst[p+n];
					lst[p+n] = temp;
				}
				p-=n;
			}
		}
	}
*/