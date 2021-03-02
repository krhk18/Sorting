/*
This program takes a integers and sort them, using different sorting algorithms:
- Bubble sort
- Insertion sort
- Shell sort
- Selection sort
- Merge sort
- Heap sort
- Quick sort

Kristina H�kansson, BTH 2019-04
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <vector>
#include <iostream>

template<typename T>
void print(T data[], int n)
{
	for(int i = 0; i < n; i++)
		std::cout << data[i] << " ";
	std::cout << std::endl;
}

//BUBBLE SORT
template<typename T>
void bubbleSort(T data[], int n)
{
    int temp;
    bool madeAswap;

    do {
        madeAswap = false;
        for(int count = 0; count < (n - 1); count++) {
            if(data[count] > data[count + 1]) {
                temp = data[count];
                data[count] = data[count + 1];
                data[count + 1] = temp;
                madeAswap = true;
            }
        }
    } while(madeAswap);
}

//INSERTION SORT
//B�sta: O(n), medel: O(n2)
//�Princip�: ta i tur och ordning ett element i taget och placera in...
//...elementet p� r�tt plats bland redan sorterade element
//(Kan vara bra n�r v�rdena �r n�stan sorterade fr�n b�rjan, eller vid l�ga n.)

template<typename T>
void insertionSort(T data[], int n)
{
	for (int i = 1; i < n; i++)
	{
		T elementToInsert = data[i];					//Spara undan talet som ska insertas
		int j = i - 1;									//B�rja j�mf�ra med talet innan
		while (j >= 0 && elementToInsert < data[j])		//Vandra ner�t, till sista eller det undansparade v�rdet inte �r st�rre l�ngre
		{
			data[j + 1] = data[j];						//Skriv �ver v�rde l�ngre till h�ger med det j�mf�rda v�rdet.
			j--;										//G� ett steg ner
		}
		data[j + 1] = elementToInsert;					//Inserta det undansparade v�rdet
	}
}

//SHELL SORT (inkrementsekvens 3, 1)
//O(n3/2)
//�Princip�: utf�r �insertionsort� men p� f�ljder d�r det finns en...
//...best�md distans mellan elementen, utf�rs upprepande med kortare...
//...och kortare distans f�r att slutligen ha distansen 1.
//K�rtiden beror p� inkrementsekvenserna (distanssekvenserna).
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
//B�sta: O(n2), Medel: O(n2)
//�Princip�: leta upp det minsta elementet bland de osorterade och...
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
//�Princip�: Dela upp arrayen (halvera) till dess att det endast finns...
//...�sm� sorterade listor och fl�ta ihop dessa till st�rre sorterade listor.
//Anv�nder extra array under sorteringen, dvs inte en in-place-algoritm
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
//In-place sortering (anv�nder inte n�gon extra datastruktur,...
//...d�remot ok att anv�nda enstaka variabler).
//�Princip�:
//1. Konstruera en (max-)heap f�r arrayen med de n element som ska sorteras (roten p� index 0)
//2. Ta bort det st�rsta fr�n den del av arrayen som representera en heap och placera elementet i arrayen...
//...direkt efter heapen totalt n g�nger, efter varje borttaging m�ste heapifering (heapify) genomf�ras
//I b�da fallen utf�rs perkolering ned�t(perculateDown)

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
//"Princip� : Utf�r grovsortering genom att �v�lja� ett element kallat pivotelement...
//...och placera alla mindre �n detta till v�nster och alla st�rre �n(eller lika med) till h�ger.
//...Detta upprepas p� v�nster och h�ger grovsorterade delar.
//Medelfall O(nlogn)
//V�rsta fall: O(n2)
//Ofta anv�nds en enklare sorteringsalgoritm (ex insertionsort) f�r grovsorterade �delar� som �r av storlek <= 10.

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
