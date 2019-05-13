/*
This program takes a integers and sort them, using different sorting algorithms:
- Insertion sort
- Shell sort
- Selection sort
- Merge sort
- Heap sort
- Quick sort

Kristina Håkansson, BTH 2019-04
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <vector>

template<typename T>
void print(T data[], int n)
{
	for(int i = 0; i < n; i++)
		std::cout << data[i] << " ";
	std::cout << std::endl;
}

//INSERTION SORT
//Bästa: O(n), medel: O(n2)
//”Princip”: ta i tur och ordning ett element i taget och placera in...
//...elementet på rätt plats bland redan sorterade element
//(Kan vara bra när värdena är nästan sorterade från början, eller vid låga n.)

template<typename T>
void insertionSort(T data[], int n)
{
	for (int i = 1; i < n; i++)
	{
		T elementToInsert = data[i];					//Spara undan talet som ska insertas
		int j = i - 1;									//Börja jämföra med talet innan
		while (j >= 0 && elementToInsert < data[j])		//Vandra neråt, till sista eller det undansparade värdet inte är större längre
		{
			data[j + 1] = data[j];						//Skriv över värde längre till höger med det jämförda värdet.
			j--;										//Gå ett steg ner
		}
		data[j + 1] = elementToInsert;					//Inserta det undansparade värdet
	}
}

//SHELL SORT (inkrementsekvens 3, 1)
//O(n3/2)
//”Princip”: utför ”insertionsort” men på följder där det finns en...
//...bestämd distans mellan elementen, utförs upprepande med kortare...
//...och kortare distans för att slutligen ha distansen 1.
//Körtiden beror på inkrementsekvenserna (distanssekvenserna).
template<typename T>
void shellSort(T data[], int n)
{
	// Gap = 3, Gap = 1
	for (int gap = 3; gap > 0; gap -= 2)
	{
		for (int i = gap; i < n; i += gap)
		{
			T temp = data[i];
			int j = i - gap;
			while (j >= 0 && temp < data[j])
			{
				data[j + gap] = data[j];
				j -= gap;
			}
			data[j + gap] = temp;
		}
	}
}

//SELECTION SORT
//Bästa: O(n2), Medel: O(n2)
//”Princip”: leta upp det minsta elementet bland de osorterade och...
//...placera det efter redan sorterade element tills alla element sorterats
template<typename T>
void selectionSort(T data[], int n)
{
	for (int i = 0; i < n; i++)
	{
		int indexOfSmallest = i;
		for (int k = i + 1; k < n; k++)
		{
			if (data[k] < data[indexOfSmallest])
				indexOfSmallest = k;
		}
		std::swap(data[i], data[indexOfSmallest]);
	}
}

//MERGE SORT
//Devide and conquer
//Recursive
//”Princip”: Dela upp arrayen (halvera) till dess att det endast finns...
//...”små” sorterade listor och fläta ihop dessa till större sorterade listor.
//Använder extra array under sorteringen, dvs inte en in-place-algoritm
//O(nlogn)
template<typename T>
void merge(T data[], int firstIndex, int lastIndex);

template<typename T>
void mergeSort(T data[], int firstIndex, int lastIndex);

template<typename T>
void mergeSort(T data[], int n) { mergeSort(data, 0, n - 1); }

template<typename T>
void mergeSort(T data[], int firstIndex, int lastIndex)
{
	if (firstIndex < lastIndex)
	{
		int midIndex = (firstIndex + lastIndex) / 2;
		mergeSort(data, firstIndex, midIndex);
		mergeSort(data, midIndex + 1, lastIndex);
		merge(data, firstIndex, lastIndex);
	}
}

template<typename T>
void merge(T data[], int firstIndex, int lastIndex)
{
	int midIndex = (firstIndex + lastIndex) / 2;
	int nL = midIndex - firstIndex + 1;				//number of elements in left vector
	int nR = lastIndex - midIndex;					//number of elements in right vector

	//create temp arrays and copy data into them
	std::vector<T> L(nL), R(nR);
	for (int i = 0; i < nL; i++)			//first array
		L[i] = data[firstIndex + i];
	for (int j = 0; j < nR; j++)			//second array
		R[j] = data[midIndex + 1 + j];

	//Start index
	int i = 0;			//first array
	int j = 0;			//second array
	int k = firstIndex;	//merged array

	//Compare element in L[] with element in R[], and insert smallest in data[] (merged array).
	while (i < nL && j < nR)
	{
		if (L[i] < R[j])
		{
			data[k] = L[i];
			i++;
		}
		else
		{
			data[k] = R[j];
			j++;
		}
		k++;
	}

	//Copy the remaining elements of L[], if there are any.
	while (i < nL)
	{
		data[k] = L[i];
		i++;
		k++;
	}

	//Copy the remaining elements of R[], if there are any.
	while (j < nR)
	{
		data[k] = R[j];
		j++;
		k++;
	}
}

//HEAP SORT
//In-place sortering (använder inte någon extra datastruktur,...
//...däremot ok att använda enstaka variabler).
//”Princip”:
//1. Konstruera en (max-)heap för arrayen med de n element som ska sorteras (roten på index 0)
//2. Ta bort det största från den del av arrayen som representera en heap och placera elementet i arrayen...
//...direkt efter heapen totalt n gånger, efter varje borttaging måste heapifering (heapify) genomföras
//I båda fallen utförs perkolering nedåt(perculateDown)

template<typename T>
void perculateDown(T data[], int n, int root);

template<typename T>
void heapSort(T data[], int n)
{
	//1. Build max heap by rearranging array
	for (int i = n / 2 - 1; i >= 0; i--)		//Perculate down every element starting with...
		perculateDown(data, n, i);							//...parent of last element.

	//2. One by one extract an element from heap 
	for (int i = n - 1; i >= 0; i--)
	{
		std::swap(data[0], data[i]);			// Swap root and end
		perculateDown(data, i, 0);				// perculate down on the reduced heap 
	}
}

template<typename T>
void perculateDown(T data[], int n, int root)
{
	int largest = root; 
	int left = 2 * root + 1; 
	int right = 2 * root + 2; 

	// If left child exists and is larger than root: that becomes new largest
	if (left < n && data[left] > data[largest])
		largest = left;

	// If right child exists and is larger than largest so far: that becomes new largest 
	if (right < n && data[right] > data[largest])
		largest = right;

	// If largest is not root 
	if (largest != root)
	{
		std::swap(data[root], data[largest]);

		// Recursively heapify the affected sub-tree 
		perculateDown(data, n, largest);
	}
}


//QUICK SORT
//Rekursiv 
//Divide-and-conquer
//"Princip” : Utför grovsortering genom att ”välja” ett element kallat pivotelement...
//...och placera alla mindre än detta till vänster och alla större än(eller lika med) till höger.
//...Detta upprepas på vänster och höger grovsorterade delar.
//Medelfall O(nlogn)
//Värsta fall: O(n2)
//Ofta används en enklare sorteringsalgoritm (ex insertionsort) för grovsorterade ”delar” som är av storlek <= 10.

template<typename T>
int partition(T data[], int firstIndx, int lastIndex);

template<typename T>
void quickSort(T data[], int firstIndex, int lastIndex);

template<typename T>
void quickSort(T data[], int n) { quickSort(data, 0, n - 1); }

template<typename T>
void quickSort(T data[], int firstIndex, int lastIndex)
{
	if (firstIndex < lastIndex)
	{
		int pivotIndex = partition(data, firstIndex, lastIndex);
		quickSort(data, firstIndex, pivotIndex);
		quickSort(data, pivotIndex + 1, lastIndex);
	}
}

template<typename T>
int partition(T data[], int firstIndex, int lastIndex)
{
	T pivotValue = data[firstIndex];
	int pivotPosition = firstIndex;

	for (int pos = firstIndex + 1; pos <= lastIndex; pos++)
	{
		if (data[pos] <= pivotValue)
		{
			std::swap(data[pivotPosition + 1], data[pos]);
			std::swap(data[pivotPosition], data[pivotPosition + 1]);
			pivotPosition++;
		}
	}
	return pivotPosition;
}


#endif
